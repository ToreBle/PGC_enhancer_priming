library(ggpubr)
library(tidyverse)
library(patchwork)
library(coin)
library(rstatix)

dir <- "PATH/PGC_enhancer_priming"

load(paste0(dir,"/Comparision/Epigenetic_PGCLC_enhancer.Rdata"))
Enhancer <- Epigenetic_PGCLC_enhancer

#Tile plot to summarize the ChIP-seq data
Active <- data.frame(mod=c(rep("H3K4me1",3), rep("H3K4me2",3), rep("H3K4me3",3)), RPGC=c(mean(Enhancer$H3K4me1_ESC), mean(Enhancer$H3K4me1_EpiLC), mean(Enhancer$H3K4me1_EpiSC), mean(Enhancer$H3K4me2_ESC), mean(Enhancer$H3K4me2_EpiLC), mean(Enhancer$H3K4me2_EpiSC), mean(Enhancer$H3K4me3_ESC), mean(Enhancer$H3K4me3_EpiLC), mean(Enhancer$H3K4me3_EpiSC)), State=c(rep(c("ESC","EpiLC","EpiSC"),3)))
Repressive <- data.frame(mod=c(rep("H3K9me2",3), rep("H3K9me3",3), rep("H3K27me2", 3), rep("H3K27me3", 3)), RPGC=c(mean(Enhancer$H3K9me2_ESC), mean(Enhancer$H3K9me2_EpiLC), mean(Enhancer$H3K9me2_EpiSC), mean(Enhancer$H3K9me3_ESC), mean(Enhancer$H3K9me3_EpiLC), mean(Enhancer$H3K9me3_EpiSC), mean(Enhancer$H3K27me2_ESC), mean(Enhancer$H3K27me2_EpiLC), mean(Enhancer$H3K27me2_EpiSC),mean(Enhancer$H3K27me3_ESC), mean(Enhancer$H3K27me3_EpiLC), mean(Enhancer$H3K27me3_EpiSC)), State=c(rep(c("ESC","EpiLC","EpiSC"),4)))

Active$State <- factor(Active$State, levels = c("EpiSC", "EpiLC", "ESC"))
Repressive$State <- factor(Repressive$State, levels = c("EpiSC", "EpiLC", "ESC"))

pdf("Enhancer_tile_plot.pdf", height = 4, width = 8)
p1 <- ggplot(Active, aes(x=mod, y=State, fill=RPGC))+geom_tile()+coord_equal()+
  theme_classic()+theme(axis.title=element_text(size=16), axis.text= element_text(size=16, color = "black"), axis.text.x = element_text(angle=45, hjust = 1), axis.line = element_blank(), legend.title = element_text(size = 12), legend.text=element_text(size=12), legend.position="top", legend.key.size = unit(0.5, "cm"))+labs(x="",y="", fill="RPGC")+
  scale_fill_gradient2(low="blue", midpoint = 6, high="red", breaks = c(3,7.5,12))
p2 <- ggplot(Repressive, aes(x=mod, y=State, fill=RPGC))+geom_tile()+coord_equal()+
  theme_classic()+theme(axis.title=element_text(size=16), axis.text= element_text(size=16, color = "black"), axis.text.x = element_text(angle=45, hjust = 1), axis.line = element_blank(), legend.title = element_text(size = 12), legend.text=element_text(size=12), legend.position="top", legend.key.size = unit(0.5, "cm"))+labs(x="",y="", fill="RPGC")+
  scale_fill_gradient2(low="blue", midpoint = 1, high="red", breaks = c(0.5, 1, 1.5))
p1+p2
dev.off()

#Example of ChIP-seq signal comparision within 1kb window and paired wilcoxon test
mCpG <- mCpG %>% dplyr::select(mCpG_ESC:mCpG_EpiSC)
rise <- nrow(mCpG)*3
mCpG[is.na(mCpG)] <- 0
mCpG <- gather(mCpG, Epigenetic, RPGC)

mCpG_comparisons <- list( c("mCpG_ESC", "mCpG_EpiLC"), c("mCpG_EpiSC", "mCpG_EpiLC"))
title <- "mCpG - in vitro"
pdf("mCpG_invitro_PGCLC_Enhancer.pdf", height = 5, width = 3)
labeling <- c("ESC", "EpiLC", "EpiSC")
ggboxplot(mCpG, x = "Epigenetic", y = "RPGC", fill= "#0101DF", xlab="", ylab="CpG methylation [%]", title = title)+stat_compare_means(method = "wilcox.test", paired=F, comparisons = mCpG_comparisons,  label = "p.signif", size=8, label.y=c(105,115))+theme_classic()+theme(legend.position="", axis.text= element_text(color = "black", size=18), axis.line = element_line(colour = "white"), axis.title=element_text(size=21), axis.text.x = element_text(angle=45, hjust = 1, size=21))+scale_x_discrete(labels= labeling)+ylim(0,119)
dev.off()


#Comparision between EpiLC (germline competent) and EpiSC (not germline competent)
#Effective size of paired wilcoxon test
d <- data.frame()
r <- data.frame()
EpiLC <- list(ATAC=Enhancer$ATAC_EpiLC)#, H3K4me1=Enhancer$H3K4me1_EpiLC, H3K4me2=Enhancer$H3K4me2_EpiLC, H3K4me3=Enhancer$H3K4me3_EpiLC, H3K27ac=Enhancer$H3K27ac_EpiLC, H3K9me2=Enhancer$H3K9me2_EpiLC, H3K9me3=Enhancer$H3K9me3_EpiLC, H3K27me2=Enhancer$H3K27me2_EpiLC, H3K27me3=Enhancer$H3K27me3_EpiLC, mCpG=mCpG$mCpG_EpiLC)
EpiSC <- list(ATAC=Enhancer$ATAC_EpiSC)#, H3K4me1=Enhancer$H3K4me1_EpiSC, H3K4me2=Enhancer$H3K4me2_EpiSC, H3K4me3=Enhancer$H3K4me3_EpiSC, H3K27ac=Enhancer$H3K27ac_EpiSC, H3K9me2=Enhancer$H3K9me2_EpiSC, H3K9me3=Enhancer$H3K9me3_EpiSC, H3K27me2=Enhancer$H3K27me2_EpiSC, H3K27me3=Enhancer$H3K27me3_EpiSC, mCpG=mCpG$mCpG_EpiSC)
for(i in 1:length(EpiLC)){
  print(names(EpiLC[i]))
  LC <- data.frame(RPGC=EpiLC[[i]], State="EpiLC")
  SC <- data.frame(RPGC=EpiSC[[i]], State="EpiSC")
  eff <- rbind(LC, SC)
  tmp2 <- wilcox_effsize(eff, RPGC~State, paired = T, ci=T, ci.type = "basic")
  tmp2$ChIP <- names(EpiLC[i])
  print(qnorm(wilcox.test(eff$RPGC~eff$State, paired = T)$p.value/2))
  if(mean(LC$RPGC/SC$RPGC, na.rm=T)<1){
    tmp2 <- tmp2 %>% mutate(effsize=effsize*(-1), conf.low=conf.low*(-1), conf.high=conf.high*(-1))}
  if(names(EpiLC[i])=="ATAC"|names(EpiLC[i])=="H3K4me1"|names(EpiLC[i])=="H3K4me2"|names(EpiLC[i])=="H3K4me3"|names(EpiLC[i])=="H3K27ac"){
    tmp2$Group <- "Active/Primed"}
  if(names(EpiLC[i])=="H3K27me2"|names(EpiLC[i])=="H3K27me3"|names(EpiLC[i])=="H3K9me2"|names(EpiLC[i])=="H3K9me3"|names(EpiLC[i])=="mCpG"){
    tmp2$Group <- "Poised/Repressive"}
  r <- rbind(r, tmp2)
}

#Arrange by effective size and plot
r$ChIP <- factor(r$ChIP, levels=(r %>% arrange(effsize))$ChIP)
pdf("Effective size PGCLC enhancer.pdf", height = 4, width = 5)
ggplot(r, aes(y=effsize, x=ChIP, fill=Group))+geom_bar(stat = "identity")+coord_flip()+scale_y_reverse()+theme_classic()+
  theme(axis.title=element_text(size=18), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=14), legend.position="top", legend.key.size = unit(0.5, "cm"))+
  scale_fill_manual(values = c("#0101DF", "#DF0101"))+labs(x="ChIP-seq samples", y="Effective Size (EpiLC vs. EpiSC)")+geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2,position=position_dodge(.9), size=1) 
dev.off()