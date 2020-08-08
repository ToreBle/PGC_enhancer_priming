library("ggplot2")
library("tidyverse")
library("rGREAT")
library("GenomicRanges")
library("ggpubr")

dir <- "PATH/PGC_enhancer_priming"

#Load PGCLC genees
PGCLC_genes <- read.table(paste0(dir,"/scRNAseq/PGCLC_genes.txt"))
colnames(PGCLC_genes) <- c("ensembl_id", "gene")

#Load ATAC and H3K27ac defined peaks to link them with PGCLC genes using GREAT
all_peaks <- read.table(paste0(dir,"/PGCLC_enhancer/PGCLC_enhancer_ATAC.bed"))
great_genes_Enhancer <- submitGreatJob(all_peaks, species = "mm10", rule="twoClosest", adv_twoDistance = 100, version = "4")

pdf("GREAT_PGCLC_enhancer_distance.pdf", height = 4, width = 5)
assignments <- plotRegionGeneAssociationGraphs(great_genes_Enhancer, type=2)
dev.off()

assignments <- as.data.frame(na.omit(assignments))
linked_enhancer <- inner_join(assignments, PGCLC_genes, by="gene")
linked_enhancer <- subset(linked_enhancer, linked_enhancer$distTSS>3500 | linked_enhancer$distTSS< (-3500))

#H3K27ac quantification with deeptools
H3K27ac_data <- read.table(paste0(dir,"/PGCLC_enhancer/GREAT_PGCLC_enhancer.tab"), skip=3)
colnames(H3K27ac_data) <- c(rep("H3K27ac_ESC", 2), rep("H3K27ac_EpiLC", 2), rep("H3K27ac_EpiSC", 2), rep("H3K27ac_d2PGCLC", 2), rep("H3K27ac_d6PGCLC", 2))
H3K27ac_data <- sapply(unique(names(H3K27ac_data)), function(col) rowMeans(H3K27ac_data[names(H3K27ac_data) == col]))
H3K27ac_data <- log2(H3K27ac_data+1)

#Clean-up of peaks with low H3K27ac
H3K27ac_data <- cbind(linked_enhancer[,1:3], H3K27ac_data)
H3K27ac_data$add <- H3K27ac_data$H3K27ac_d2PGCLC+H3K27ac_data$H3K27ac_d6PGCLC
ggplot(H3K27ac_data, aes(x=H3K27ac_EpiLC, y=H3K27ac_d2PGCLC, color=add))+geom_point() + theme_classic()+
  theme(axis.title=element_text(size=24), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 16), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))
H3K27ac_enhancer <- subset(H3K27ac_data, H3K27ac_data$add>4.42)
H3K27ac_enhancer$add <- NULL

#Linkage of enhancers to the assigned genes after clean-up
colnames(H3K27ac_enhancer) <- c("chr", "start", "end", "ESC", "EpiLC", "EpiSC", "d2PGCLC", "d6PGCLC")
colnames(linked_enhancer) <- c("chr", "start", "end", "width", "strand", "gene", "distTSS", "ens_gene")
H3K27ac_enhancer <- makeGRangesFromDataFrame(H3K27ac_enhancer, keep.extra.columns=T, ignore.strand=T)
linked_enhancer <- makeGRangesFromDataFrame(as.data.frame(linked_enhancer), keep.extra.columns=T, ignore.strand=F)
Enhancer_linked_genes <- subsetByOverlaps(linked_enhancer, H3K27ac_enhancer)
write.table(Enhancer_linked_genes, row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(dir,"/PGCLC_enhancer/PGCLC_enhancer_assigned.bed"))


#Distance of the PGCLC enhancers to their assigned target genes
dist <- data.frame(Enhancer_linked_genes)
ggplot(dist, aes(distTSS/1000))+geom_histogram(binwidth = 5, fill="darkblue")+ theme_classic()+theme(legend.position="", axis.text.x= element_text(color = "black", size=18, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=14), axis.line = element_line(colour = "white"), axis.title=element_text(size=18), plot.title = element_text(size=16))+labs(x="Distance to PGCLC gene (kb)", y="Absoulte PGCLC Enhnacer")

#Counting of the Enhancers per PGCLC gene
count_enhancers <- data.frame()
for(genes in PGCLC_genes$gene){
  tmp <- data.frame(gene=genes, Enhancer=length(subset(Enhancer_linked_genes, Enhancer_linked_genes$gene==genes)))
  count_enhancers <- rbind(count_enhancers, tmp)
}

pdf(paste("PGCLC_enhancer_density.pdf", sep=""), height = 4, width = 5)
ggplot(count_enhancers, aes(x=(Enhancer))) + geom_density(color="white", fill="darkblue", size=2)+ theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=18, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 16), legend.position="", legend.text=element_text(size=16))+labs(x="Enhancer / Gene", y="Density", title = paste("Enhancer: ",length(Enhancer_linked_genes), ", PGCLC genes with enhancers: ", length(unique(Enhancer_linked_genes$gene)),sep=""))
dev.off()
