library(biomaRt)
library(edgeR)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(Rtsne)

dir <- "PATH/PGC_enhancer_priming"

#scRNAseq data from https://doi.org/10.1038/s41586-019-1825-8
scRNA_counts <- read.table(paste0(dir,"/Heterogeneity/GSE121650_rna_counts.tsv"), header=T, row.names = "ensembl_id")


#RPKM normalization
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "http://dec2016.archive.ensembl.org")
s=getBM(attributes=c("ensembl_gene_id", "start_position","end_position"), mart=ensembl)
s$length=s$end_position - s$start_position
detach("package:biomaRt", unload=TRUE)

genes <- s %>% dplyr::filter(ensembl_gene_id %in% rownames(scRNA_counts))
genes <- genes[order(rownames(scRNA_counts), genes$ensembl_gene_id),]
sc_rpkm <- rpkm(scRNA_counts, genes$length)


#Select only assigned epiblast cells from meta data: https://github.com/rargelaguet/scnmt_gastrulation
meta <- read.table(paste0(dir,"/Heterogeneity/sample_metadata.txt"), header=T)
meta_epi <- subset(meta, meta$lineage10x=="Epiblast" & meta$stage!="E7.5")
meta_epi$sample <- as.character(meta_epi$sample)


#tSNE plot of the scRNAseq data
set.seed(1989)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(sc_rpkm))
tsne_prep <- data.frame(sc_rpkm) %>% dplyr::select(epi$sample)
d <- dist(t(tsne_prep))
r <- Rtsne(d, is.matrix=T)
tnse_plot <- data.frame(x = r$Y[,1], y = r$Y[,2], State = epi$stage)
pdf("tSNE_Epiblast.pdf", height = 5, width = 5.5)
ggplot(tnse_plot, aes(y, x, colour = State))+geom_point(size=2) +  scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"), "Stage") + 
  labs(x = "Component 2: t-SNE", y = "Component 1: t-SNE")+
  theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()

#Transcriptional Noise
#-Transfer a correlation matrix into a data.frame
flattenCorrelationMatrix <- function(cmatrix, pvalue) {
  ut <- upper.tri(cmatrix)
  data.frame(
    row = rownames(cmatrix)[row(cmatrix)[ut]],
    column = rownames(cmatrix)[col(cmatrix)[ut]],
    cor  =(cmatrix)[ut],
    p = pvalue[ut]
  )
}


#Calculation of the transcriptional noise for each stage by:
#Selection of the stage
#Calculation the variance for all genes of a stage and select the top 500 (top_n)
#Spearman correlation of all highly varibale genes
#Transformation of the Spearman correlation into the transcriptional noise by sqrt((1-cor)/2)

States <- c("E4.5", "E5.5", "E6.5")
noise <- data.frame()

for(i in States){
  print(i)
  all <- data.frame(sc_rpkm)  %>% dplyr::select((subset(epi, epi$stage==i)$sample))
  all$var <- apply(all,1,var)
  all <- top_n(all, 500, var) %>% dplyr::select(-var)
  corr <- rcorr(as.matrix(all), type = c("spearman"))
  c <- flattenCorrelationMatrix(corr$r, corr$P) %>% dplyr::mutate(State=i, cor=sqrt((1-cor)/2))
  noise <- rbind(noise, c)
}

pdf("Transcriptional_epiblast.pdf", height=4, width = 4)
noise <- subset(noise, noise$p==0)
ggplot(noise, aes(State,cor, fill=State))+geom_violin()+geom_boxplot(width=0.3, outlier.shape = NA, fill="white")+
  labs(x = "", y = "Transcriptional noise")+
  theme_classic()+scale_fill_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=14, color = "black"), axis.text.x = element_text(size = 18), axis.line = element_blank(), legend.position="")+
  stat_compare_means(ref.group = "E5.5", label = "p.signif", label.y = 0.58, size=14, symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("*", "*C","*B","*A", "ns")))+ylim(0.25,0.65)
dev.off()

#Expression of PGCLC genes within single-cell epiblast
PGCLC_genes <- read.table(paste0(dir,"/scRNAseq/PGCLC_genes.txt"))
colnames(PGCLC_genes) <- c("ensembl_id", "gene")

scPGC_genes <- as.data.frame(sc_rpkm) %>% dplyr::filter(rownames(sc_rpkm) %in% PGCLC_genes$ensembl_id)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(scPGC_genes))
scPGC <- scPGC_genes %>% dplyr::select(epi$sample)
save(scPGC, file=paste0(dir,"/Heterogeneity/scPGC.Rdata"))

genes_expression <- data.frame()
cell_expression <- data.frame()
States <- c("E4.5", "E5.5", "E6.5")
for(i in States){
  print(i)
  PGCLC_expression <- scPGC %>% dplyr::select((subset(epi, epi$stage==i)$sample))  
  subset_genes <- data.frame(mean = rowMeans(PGCLC_expression), Stage = i)
  subset_cell <- data.frame(mean = rowMeans(t(PGCLC_expression)), Stage = i)
  genes_expression <- rbind(genes_expression, subset_genes)
  cell_expression <- rbind(cell_expression, subset_cell)
}

genes_expression$log <- log2(genes_expression$mean+1)

#Removing a single outlier in E4.5
cell_expression <- subset(cell_expression, cell_expression$mean<40)

color_dynamics <- c("#0101DF", "#DF0101", "#FF8000")
pdf("Expression_dynamics_epiblast.pdf", height=4, width = 4)
ggboxplot(genes_expression, x = "Stage", y = "log", color = "Stage", add="jitter", outlier.shape=NA)+labs(x =" ", y="log2 Expression per gene")+scale_color_manual(values=color_dynamics)+theme_classic()+
  theme(legend.position="", axis.text.x= element_text(color = "black", size=18), axis.text.y= element_text(color = "black", size=14), axis.line = element_line(colour = "white"), axis.title=element_text(size=18), plot.title = element_text(size=16))
ggboxplot(cell_expression, x = "Stage", y = "mean", color = "Stage", add="jitter", outlier.shape=NA)+labs(x =" ", y="PGCLC genes - mean rpkm")+scale_color_manual(values=color_dynamics)+theme_classic()+
  theme(legend.position="", axis.text.x= element_text(color = "black", size=18), axis.text.y= element_text(color = "black", size=18), axis.line = element_line(colour = "white"), axis.title=element_text(size=18), plot.title = element_text(size=16))
dev.off()


#tSNE plot for PGCLC genes only
set.seed(1989)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(scPGC_genes))
tsne_prep <- data.frame(scPGC_genes) %>% dplyr::select(epi$sample)
d <- dist(t(tsne_prep))
r <- Rtsne(d, is.matrix=T)
tnse_plot <- data.frame(x = r$Y[,1], y = r$Y[,2], State = epi$stage)
pdf("tSNE_PGCLC_genes_Epiblast.pdf", height = 5, width = 5.5)
ggplot(tnse_plot, aes(y, x, colour = State))+geom_point(size=2) +  scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"), "Stage") + 
  labs(x = "Component 2: t-SNE", y = "Component 1: t-SNE")+
  theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()