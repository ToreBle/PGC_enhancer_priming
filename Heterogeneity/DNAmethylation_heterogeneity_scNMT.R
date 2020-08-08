library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggpubr)
library(scales)

dir <- "PATH/PGC_enhancer_priming"

#Preparation of the meta data from https://github.com/rargelaguet/scnmt_gastrulation and public GEO ids:
meta <- read.table(paste0(dir,"/Heterogeneity/sample_metadata.txt"), header=TRUE)
meta_epi <- subset(meta, meta$lineage10x=="Epiblast" & meta$stage!="E7.5")

GEO_id <- read.table(paste0(dir,"/Heterogeneity/GEO_id_met.txt"), header=TRUE)
meta_epi <- inner_join(meta_epi, GEO_id, by="sample")
rownames(meta_epi) <- meta_epi$id
meta_epi <- meta_epi %>% dplyr::select(sample, stage, id)

load(paste0(dir,"/Heterogeneity/DNAmethylation_global.Rdata"))
DNAmethylation_global <- meta_data

load(paste0(dir,"/Heterogeneity/DNAmethylation_PGCLC_enhancer.Rdata"))
mCpG_comparision <- merge(meta_data, DNAmethylation_global, by=0, all=T)
rownames(mCpG_comparision) <- mCpG_comparision$Row.names
mCpG_comparision$Row.names <- NULL

meta_data <- merge(mCpG_comparision, meta_epi, by=0)
rownames(meta_data) <- meta_data$Row.names
meta_data$Row.names <- NULL


#Correlation of CpG methylation at enhancers vs genome-wide and 
#QC for the coverage at PGCLC enhancers
pdf("mCpG_genome_wide_and_cov_correlation.pdf", height = 4, width = 5)
p1 <- ggscatter(meta_data, x = "mCpG_local", y = "mCpG_global", xlab = "PGCLC Enhancer - mCpG/CpG", ylab = "Genome wide - mCpG/CpG", na.rm=T, color = "stage")+
  scale_color_manual(values=c("#0101DF", "#DF0101", "#FF8000"))+theme_classic()+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))+
  stat_cor(method = "spearman", size=5)+guides(colour = guide_legend(override.aes = list(size=6)))

p2 <- ggscatter(meta_data, x = "mCpG_local", y = "coverage_local",
                xlab = "PGCLC Enhancer - mCpG/CpG", ylab = "PGCLC enhancer - CpG coverage", na.rm=T, color = "stage")+
  scale_color_manual(values=c("#0101DF", "#DF0101", "#FF8000"))+theme_classic()+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))+
  stat_cor(method = "spearman", size=5, label.y=3.60)+guides(colour = guide_legend(override.aes = list(size=6)))
print(p1)
print(p2+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))))
dev.off()

#Preparation for the ordering of the dissimility matrix
# by increasing CpG methylation at PGCLC enhancers in each stage
E4_naive <- meta_data %>%
  filter(stage=="E4.5") %>%
  arrange(mCpG_local)
rownames(E4_naive) <- E4_naive$id

E5_formative <- meta_data %>%
  filter(stage=="E5.5") %>%
  arrange(mCpG_local)
rownames(E5_formative) <- E5_formative$id

E6_primed <- meta_data %>%
  filter(stage=="E6.5") %>%
  arrange(mCpG_local)
rownames(E6_primed) <- E6_primed$id

meta_diss <- rbind(E4_naive, E5_formative, E6_primed)

#Plotting of the prepared Dissimilarity matrix (from Heterogeneity/Dissimilarity_matrices.R)
load(paste0(dir,"/Heterogeneity/diss_scM_Epiblast_PGCLC.Rdata"))
diss_ordered <- diss_ordered[match(rownames(meta_diss), rownames(diss_ordered)), ]
diss_ordered <- diss_ordered[,match(rownames(meta_diss), colnames(diss_ordered))]

#Preparing the meta data beside the main matrix
meta_diss <- meta_diss %>% dplyr::select(mCpG_local, stage)
anno <- HeatmapAnnotation(df = meta_diss, 
                          col = list(stage = c("E4.5" = "#0101DF", "E5.5" = "#DF0101", "E6.5"="#FF8000"),
                                     mCpG_local = colorRamp2(c(0, 0.6, 0.7), c("blue", "#FFFFFF", "orange"))),
                          
                          annotation_legend_param = list(stage = list(nrow=1, title = "Stage", title_position = "topcenter", title_gp = gpar(fontsize = 14),
                                                                      labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"), space = unit(3, "mm")),
                                                         mCpG_local = list(direction = "horizontal", title = "mCpG/CpG", title_position = "leftcenter", title_gp = gpar(fontsize = 14), 
                                                                           labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"))),
                          show_annotation_name = FALSE, show_legend =c(TRUE, FALSE))

pdf("Dissimilarity_matrix_PGCLC_enhancer_Epiblast.pdf", height = 5.5, width = 6)
enhancer_diss <- Heatmap(diss_ordered, show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F, top_annotation = anno, name="Dissimilarity", col = colorRamp2(c(0.2,0.4,0.6), c("#0101DF", "white", "#DF0101")), heatmap_legend_param = list(title = "mCpG heterogeneity", title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12), title_position = "leftcenter-rot", legend_height = unit(4.5, "cm")))#use_raster = T to reduce the image complexity
draw(enhancer_diss, annotation_legend_side = "top")#merge_legend = TRUE)
dev.off()

#Epigenetic heterogeneity of the PGCLC enhancers
#Transfer a correlation matrix into a data.frame for the comparision of CpG methylation heterogeneity
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

Heterogeneity_PGCLC_enhancer <- flattenCorrMatrix(diss_ordered)
E4 <- Heterogeneity_PGCLC_enhancer %>% filter(row %in% E4_naive$id & column %in% E4_naive$id) %>% mutate(Stage = "E4.5")
E5 <- Heterogeneity_PGCLC_enhancer %>% filter(row %in% E5_formative$id & column %in% E5_formative$id) %>% mutate(Stage = "E5.5")
E6 <- Heterogeneity_PGCLC_enhancer %>% filter(row %in% E6_primed$id & column %in% E6_primed$id) %>% mutate(Stage = "E6.5")
epigenetic_score_enhancer <- rbind(E4, E5, E6)

pdf("Epigenetic heterogeneity PGCLC enhancer.pdf", height = 4, width = 5)
ggplot(epigenetic_score_enhancer, aes(x=Stage, y=cor, fill=Stage))+geom_violin()+geom_boxplot(width=0.3, fill="white")+theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(0.5, "cm"))+
  labs(x="", y="Epigenetic heterogeneity")+scale_fill_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+
  stat_compare_means(ref.group = "E5.5", label = "p.signif", size=14, symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("*", "*C","*B","*A", "ns")), label.y = 0.85)+scale_y_continuous(labels = scales::percent, limits = c(0,0.95))
dev.off()


#Genome-wide CpG methylation heterogeneity
load(paste0(dir,"/Heterogeneity/diss_scM_Epiblast_genome.Rdata"))
Stages <- c("E4.5", "E5.5", "E6.5")
meta_data_genome_diss <- data.frame()
for(Stage in Stages){
  genome <- meta_data %>%
    filter(stage==Stage) %>%
    arrange(mCpG_global) %>%
    dplyr::select(mCpG_global, stage, id)
  rownames(genome) <- genome$id
  genome$id <- NULL
  meta_data_genome_diss <- rbind(meta_data_genome_diss, genome)
}

diss_ordered <- diss_ordered[match(rownames(meta_data_genome_diss), rownames(diss_ordered)), ]
diss_ordered <- diss_ordered[,match(rownames(meta_data_genome_diss), colnames(diss_ordered))]

anno <- HeatmapAnnotation(df = meta_data_genome_diss, 
                          col = list(stage = c("E4.5" = "#0101DF", "E5.5" = "#DF0101", "E6.5"="#FF8000"),
                                     mCpG_global = colorRamp2(c(0, 0.2, 0.7, 0.8), c("darkblue", "blue", "#FFFFFF", "orange"))),
                          
                          annotation_legend_param = list(stage = list(nrow=1, title = "Stage", title_position = "topcenter", title_gp = gpar(fontsize = 14),
                                                                      labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"), space = unit(3, "mm")),
                                                         mCpG_global = list(direction = "horizontal", title = "mCpG/CpG", title_position = "leftcenter", title_gp = gpar(fontsize = 14), 
                                                                            labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"))),
                          show_annotation_name = FALSE, show_legend =c(TRUE, FALSE))

pdf("Dissimilarity_matrix_global.pdf", height = 5.2, width = 6.2)
genome_diss <- Heatmap(diss_ordered, show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F, top_annotation = anno, name="Dissimilarity", col = colorRamp2(c(0,0.2,0.5,0.6), c("darkblue", "#0101DF", "white", "#DF0101")), heatmap_legend_param = list(title = "mCpG heterogeneity", title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12), title_position = "leftcenter-rot", legend_height = unit(4.5, "cm")))#use_raster = T to reduce the image complexity
draw(genome_diss, annotation_legend_side = "top")#merge_legend = TRUE)
dev.off()

#Genom-wide heterogeneity score
Heterogeneity_genome <- flattenCorrMatrix(diss_ordered)
E4 <- Heterogeneity_genome  %>% filter(row %in% E4_naive$id & column %in% E4_naive$id) %>% mutate(Stage = "E4.5")
E5 <- Heterogeneity_genome %>% filter(row %in% E5_formative$id & column %in% E5_formative$id) %>% mutate(Stage = "E5.5")
E6 <- Heterogeneity_genome %>% filter(row %in% E6_primed$id & column %in% E6_primed$id) %>% mutate(Stage = "E6.5")
heterogeneity_score_genome <- rbind(E4, E5, E6)

pdf("Epigenetic heterogeneity genome.pdf", height = 4, width = 4)
ggplot(heterogeneity_score_genome, aes(x=Stage, y=cor, fill=Stage))+geom_violin()+geom_boxplot(width=0.3, fill="white")+theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(0.5, "cm"))+
  labs(x="", y="Epigenetic heterogeneity")+scale_fill_manual(values = c("#0101DF", "#DF0101", "#FF8000", "green"))+
  stat_compare_means(ref.group = "E5.5", label = "p.signif", size=14, symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("*", "*C","*B","*A", "ns")), label.y = 0.4)+scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.45))
dev.off()

#Comparsion of different mCpG Heterogeneities within PGCLC, EpiLC enhancers or genome-wide
Enhancer <- c("PGCLC", "genome", "EpiLC")
Heterogeneity <- data.frame()
for(e in Enhancer){
  load(paste0(dir, "/Heterogeneity/diss_scM_Epiblast_", e, ".Rdata"))
  diss_ordered <- diss_ordered[match(rownames(meta_data_genome_diss), rownames(diss_ordered)), ]
  diss_ordered <- diss_ordered[,match(rownames(meta_data_genome_diss), colnames(diss_ordered))]
  tmp <- flattenCorrMatrix(diss_ordered)
  E4 <- tmp %>% filter(row %in% E4_naive$id & column %in% E4_naive$id) %>% mutate(Stage = "E4.5", Enhancer= e)
  E5 <- tmp %>% filter(row %in% E5_formative$id & column %in% E5_formative$id) %>% mutate(Stage = "E5.5", Enhancer= e)
  E6 <- tmp %>% filter(row %in% E6_primed$id & column %in% E6_primed$id) %>% mutate(Stage = "E6.5", Enhancer= e)
  all <- rbind(E4, E5, E6)
  Heterogeneity <- rbind(Heterogeneity, all)
}

Heterogeneity$Enhancer <- factor(Heterogeneity$Enhancer, levels = c("genome", "EpiLC", "PGCLC"))

pdf("Epigenetic heterogeneity Comparision.pdf", height = 5, width = 6)
ggline(Heterogeneity, x="Stage", y="cor", color="Enhancer", add="mean_ci", error.plot="errorbar", size=2)+theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=16), legend.key.height=unit(2, "cm"), legend.position="right", legend.key.size = unit(1, "cm"))+
  labs(x="", y="Epigenetic heterogeneity")+scale_color_manual(labels= c("genome-\nwide", "EpiLC\nenhancer", "PGCLC\nenhancer"), values = c("darkgray", "#FF8000", "#DF0101"))+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggplot(Heterogeneity, aes(x=Stage, y=cor, fill=Enhancer))+geom_boxplot()+theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))+
  labs(x="", y="Epigenetic heterogeneity")+scale_fill_manual(values = c("#0101DF", "#DF0101", "#FF8000", "green"))+
  scale_y_continuous(labels = scales::percent, limits = c(0,0.8))
dev.off()

#Group each stage in low and low high CpG methylation depending on
#whether the mean methylation of all PGCLC enhancers is above or beneath
#the average for that stage
Stages <- c("E4.5", "E5.5", "E6.5")
meta_data_enhancer <- data.frame()
for(Stage in Stages){
  enhancer <- meta_data %>%
    filter(stage==Stage) %>%
    arrange(mCpG_local)
  #enhancer <- enhancer %>% mutate(Group = c(rep(paste(Stage,"low", sep="_"), nrow(subset(enhancer, enhancer$mCpG_local<mean(enhancer$mCpG_local)))), rep(paste(Stage,"high", sep="_"), nrow(subset(enhancer, enhancer$mCpG_local>mean(enhancer$mCpG_local)))))) %>%
  enhancer <- enhancer %>% mutate(Group = c(rep("low", nrow(subset(enhancer, enhancer$mCpG_local<mean(enhancer$mCpG_local)))), rep("high", nrow(subset(enhancer, enhancer$mCpG_local>mean(enhancer$mCpG_local)))))) %>% dplyr::select(mCpG_local, stage, sample, Group, id)
  rownames(enhancer) <- enhancer$id
  enhancer$id <- NULL
  meta_data_enhancer <- rbind(meta_data_enhancer, enhancer)
}


#Classificiation of the data in respect to low or high methylation:
#for CpG methylation, Chromatin accessiblity at PGCLC enhancers
#and the expression of PGCLC

#Starting with DNA methylation
scM <- meta_data_enhancer
scM$Group <- factor(scM$Group, levels = c("low", "high"))

pdf("mCpG PGCLC enhancer over time.pdf", height = 4.5, width = 5)
ggplot(scM, aes(x=Group, y=mCpG_local, color=stage, group=Group))+facet_grid(.~stage)+
  geom_boxplot(outlier.shape = NA)+geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1))+theme_classic()+
  theme(strip.text.x = element_text(size=18), strip.background = element_rect(colour="white"))+annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(1, "cm"))+
  labs(x="CpG methylation", y="DNA methylation - mCpG")+scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
dev.off()

#Continuing with scN (Chromatin accessibility - mGpC)
scN_header <- read.table(paste0(dir, "/Heterogeneity/scM_scN_acc_1kb.tab"), nrows = 1, header = FALSE)
scN_header$V1 <- NULL
scN_tab <- read.table(paste0(dir, "/Heterogeneity/scM_scN_acc_1kb.tab"), skip=3)
colnames(scN_tab) <- apply(scN_header, 1, function(x) substring(x,1,10))
data <- sapply(unique(colnames(scN_tab)), function(i) rowMeans(scN_tab[,colnames(scN_tab) == i], na.rm=T))

Stages <- c("E4.5", "E5.5", "E6.5")
scN <- data.frame()
for(i in Stages){
  tmp <- data.frame()
  for(tide in c("low", "high")){
    sub <- subset(meta_data_enhancer, meta_data_enhancer$Group==tide && meta_data_enhancer$stage==i)
    tmp1 <- data[,match(rownames(sub), colnames(data))]
    tmp1 <- data.frame(GpC=colMeans(tmp1, na.rm = T), Stage=sub$stage, group=sub$Group)
    tmp <- rbind(tmp, tmp1)}
  scN <- rbind(scN, tmp)}

#Removing outliers (that are probably technical because only cells with a CpG coverage has been selected)
scN <- subset(scN, scN$GpC>0.21 & scN$GpC<0.6)
scN$group <- factor(scN$group, levels = c("low", "high"))

pdf("Acessibility of PGCLC enhancer over time.pdf", height = 4.5, width = 5.5)
ggplot(scN, aes(x=group, y=GpC, color=Stage, group=group))+facet_grid(.~Stage)+
  geom_boxplot(outlier.shape = NA)+geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1))+theme_classic()+
  theme(strip.text.x = element_text(size=18), strip.background = element_rect(colour="white"))+annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(1, "cm"))+
  labs(x="CpG methylation", y="Chromatin accessibility - mGpC")+scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+stat_compare_means(aes(label = sprintf("p = %5.3f", round(as.numeric(..p.format..), digit=3))), size=5, ref.group = "low", label.y = 0.54)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
dev.off()

#Final completing the analysis with the expression analysis of the PGCLC genes
load(paste0(dir,"/Heterogeneity/scPGC.Rdata")) #From scRNAseq_analysis
Stages <- c("E4.5", "E5.5", "E6.5")
PGCLC_genes <- data.frame()
for(i in Stages){
  tmp <- data.frame()
  for(tide in c("low", "high")){
    sub <- subset(meta_data_enhancer, meta_data_enhancer$Group==tide && meta_data_enhancer$stage==i)
    tmp1 <- scPGC[,match(sub$sample, colnames(scPGC))]
    tmp1 <- data.frame(Expression=colMeans(tmp1, na.rm = T), Stage=sub$stage, group=sub$Group)
    tmp <- rbind(tmp, tmp1)}
  PGCLC_genes <- rbind(PGCLC_genes, tmp)}

#Removing outliers (that are probably technical because only cells with a CpG coverage has been selected)
PGCLC_genes <- subset(PGCLC_genes, PGCLC_genes$Expression>0 & PGCLC_genes$Expression<15)
PGCLC_genes$group <- factor(PGCLC_genes$group, levels = c("low", "high"))

pdf("Expression of PGCLC enhancer over time.pdf", height = 4.5, width = 5)
ggplot(PGCLC_genes, aes(x=group, y=Expression, color=Stage, group=group))+facet_grid(.~Stage)+
  geom_boxplot(outlier.shape = NA)+geom_jitter(position = position_jitterdodge(jitter.width = 1, dodge.width = 1))+theme_classic()+
  theme(strip.text.x = element_text(size=18), strip.background = element_rect(colour="white"))+annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(1, "cm"))+
  labs(x="CpG methylation", y="Gene expression - rpkm")+scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+stat_compare_means(label = "p.format", size=5, ref.group = "low", label.y = 14.5)
dev.off()
