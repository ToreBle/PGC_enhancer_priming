library("ggplot2")
library("Seurat")
library("ggpubr")

dir <- "PATH/PGC_enhancer_priming"

#Load complete UMI count matrix and define cell stages
Expression <- read.csv(paste0(dir,"/scRNAseq/UMI_counts.csv"), row.names = 1)

ESC <- Expression[,c(grep(".1", names(Expression)))]
SL <- Expression[,c(grep(".7", names(Expression)))]
d1_EpiLC <- Expression[,c(grep(".3", names(Expression)))]
d2_EpiLC <- Expression[,c(grep(".2", names(Expression)))]
d3_EpiLC <- Expression[,c(grep(".8", names(Expression)))]
EpiSC <- Expression[,c(grep(".5", names(Expression)))]


#Define d2/d4 PGCLC cluster and the EB cluster (non-PGCLC cluster of d2/d4 EB)
d2_EB_all <- Expression[,c(grep(".4", names(Expression)))]
d4_EB_all <- Expression[,c(grep(".6", names(Expression)))]

#Define d2 PGCLC
d2_PGCLC_cluster <- read.csv(paste0(dir,"/scRNAseq/d2_PGCLC_cluster.csv")) %>% dplyr::select(-row.mean, -X)
d2_PGCLC <- Expression %>% dplyr::select(colnames(d2_PGCLC_cluster))
d2_EB <- d2_EB_all %>% dplyr::select(-one_of(names(d2_PGCLC)))

#Define d4 PGCLC
d4_PGCLC_cluster <- read.csv(paste0(dir,"/scRNAseq/d4_PGCLC_cluster.csv")) %>% dplyr::select(-row.mean, -X)
d4_PGCLC <- Expression %>% dplyr::select(colnames(d4_PGCLC_cluster))
d4_EB <- d4_EB_all %>% dplyr::select(-one_of(names(d4_PGCLC)))


#Creating the metadata for Seurat
ES <- data.frame(id = colnames(ESC), Stage = "others", Sample= "ESC")
SL <- data.frame(id = colnames(SL), Stage = "others", Sample= "SL")
d1 <- data.frame(id = colnames(d1_EpiLC), Stage = "others", Sample= "d1 EpiLC")
d2 <- data.frame(id = colnames(d2_EpiLC), Stage = "negative", Sample= "d2 EpiLC")
d3 <- data.frame(id = colnames(d3_EpiLC), Stage = "others", Sample= "d3 EpiLC")
SC <- data.frame(id = colnames(EpiSC), Stage = "others", Sample= "EpiSC")

PGC_negative1 <- data.frame(id = colnames(d2_EB), Stage = "negative", Sample= "d2 EB")
PGCLC1 <- data.frame(id = colnames(d2_PGCLC), Stage = "PGCLC", Sample= "d2 EB")
PGC_negative2 <- data.frame(id = colnames(d4_EB), Stage = "negative", Sample= "d4 EB")
PGCLC2 <- data.frame(id = colnames(d4_PGCLC), Stage = "PGCLC", Sample= "d4 EB")
  
meta_seurat_PGCLC <- rbind(PGC_negative1, PGC_negative2, PGCLC1, PGCLC2, ES, SL, d1, d2, d3, SC)
rownames(meta_seurat_PGCLC) <- meta_seurat_PGCLC$id
save(meta_seurat_PGCLC, file=paste0(dir,"/scRNAseq/meta_scRNAseq_PGCLC.Rdata"))


#Seurat object to select upregulated PGCLC genes (vs. remaining EB and d2 EpiLC)
Seurat <- CreateSeuratObject(Expression, project = "PGCLC", meta.data = meta_seurat_PGCLC)
Seurat <- NormalizeData(Seurat)

upregulated <- FindMarkers(Seurat, group.by = Seurat$Stage, ident.1 = "PGCLC", ident.2 = "negative", test.use = "negbinom", only.pos = TRUE)
sub <- subset(upregulated, upregulated$p_val_adj<0.001)
sub <- subset(sub, sub$pct.2<0.42)


#Create table of all PGCLC genes
library("biomaRt")
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "http://mar2016.archive.ensembl.org")
PGCLC_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id", values=rownames(sub), mart=ensembl)
detach("package:biomaRt", unload=TRUE)

PGCLC_genes <- dplyr::rename(PGCLC_genes, ens_gene = ensembl_gene_id, gene = external_gene_name)
write.table(PGCLC_genes, row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(paste0(dir,"/scRNAseq/PGCLC_genes.txt"), sep=""))

#Expression of all PGCLC genes as mean expression per cell (Figure 1e) or mean expression in all cells of a stage (Supplementary Figure 1c)
List <- list(ESC=ESC,d1_EpiLC=d1_EpiLC, d2_EpiLC=d2_EpiLC, d3_EpiLC=d3_EpiLC, EpiSC=EpiSC, d2_EB=d2_EB, d2_PGCLC=d2_PGCLC, d4_EB=d4_EB, d4_PGCLC=d4_PGCLC)
rownames(PGCLC_genes) <- PGCLC_genes$ens_gene

genes_expression <- data.frame()
cell_expression <- data.frame()

for(i in 1:length(List)){
  subset <- List[[i]][rownames(PGCLC_genes),]
  subset_genes <- data.frame(mean = rowMeans(subset), State = names(List[i]))
  subset_cell <- data.frame(mean = rowMeans(t(subset)), State = names(List[i]))
  genes_expression <- rbind(genes_expression, subset_genes)
  cell_expression <- rbind(cell_expression, subset_cell)
}

color_dynamics <- c("#0101DF", "#8904B1", "#DF0101", "#B40404", "#FF8000", "darkgray", "#01DF01", "gray", "#088A08")
labels <- c("ESC" = "ESC", "d1_EpiLC" = "d1 EpiLC", "d2_EpiLC" = "d2 EpiLC", "d3_EpiLC" = "d3 EpiLC", "EpiSC" = "EpiSC", "d2_EB" = "d2 EB", "d2_PGCLC" = "d2 PGCLC", "d4_EB" = "d4 EB", "d4_PGCLC" = "d4 PGCLC")

pdf(paste("Expression_dynamics_PGCLC.pdf", sep=""), height=5, width = 8)
genes_expression$log <- log2(genes_expression$mean+1) 
ggboxplot(genes_expression, x = "State", y = "log", color = "State", add="jitter", outlier.shape=NA)+labs(x =" ", y=expression(paste("log"[2],"(Expression per gene + 1)", sep="")))+scale_color_manual(values=color_dynamics)+theme_classic()+
  theme(legend.position="", axis.text.x= element_text(color = "black", size=18, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=14), axis.line = element_line(colour = "white"), axis.title=element_text(size=18), plot.title = element_text(size=16))+
  scale_x_discrete(labels=labels)
ggboxplot(cell_expression, x = "State", y = "mean", color = "State", add="jitter", outlier.shape=NA)+labs(x =" ", y="Expression per cell")+scale_color_manual(values=color_dynamics)+theme_classic()+
  theme(legend.position="", axis.text.x= element_text(color = "black", size=21, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=18), axis.line = element_line(colour = "white"), axis.title=element_text(size=21), plot.title = element_text(size=16))+
  scale_x_discrete(labels=labels)
dev.off()