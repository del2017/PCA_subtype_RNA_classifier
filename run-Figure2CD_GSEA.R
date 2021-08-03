library(MASS)
library(RColorBrewer)
library(eulerr)
library(RColorBrewer)
library(gdata)
library(stringr)
library(Rtsne)
library(ggrepel)
library(ggpmisc)

#Color schemes
col.DEG.both <- brewer.pal(9,"Set1")[1]
col.DEG.1 <- brewer.pal(9,"Set1")[2]
col.DEG.2 <- brewer.pal(9,"Set1")[3]
col.DEG.no <- adjustcolor( brewer.pal(12,"Set3")[9] , alpha.f = 0.75)
pval_cutoff <- 0.05

#Figure s4: Correlation between SPOP_CHD1wt/N and ERG_PTENwt/N in TCGA
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=10, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1wt_N", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENwt_N", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP_CHD1_ERG_PTEN_N_hallmark-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1wt_N, y=NES.ERG.PTENwt_N)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1wt_N+0.000000001), colour=-log10(NOM.p.val.ERG.PTENwt_N+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()


#Figure 2C: Correlation between SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt in TCGA
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=11, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP_CHD1_ERG_PTEN_wt_hallmark-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()


#Figure 2C: Correlation between SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt in Taylor
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=12, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_SPOP_CHD1_ERG_PTEN_wt_hallmark-cor.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()


#Figure S5: Correlation between SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt in ICGC
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=13, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
GSEA_merge <- merge(GSEA_ICGC_ERG_PTENdel_PTENwt[,c(1,6,7,8)], GSEA_ICGC_SPOP_CHD1del_CHD1wt[,c(1,6,7,8)], by="NAME")
colnames(GSEA_merge) <- c(colnames(GSEA_merge)[1], 
                          paste(colnames(GSEA_ICGC_ERG_PTENdel_PTENwt[, c(6,7,8)]), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(colnames(GSEA_ICGC_SPOP_CHD1del_CHD1wt[, c(6,7,8)]), ".ERG.PTENdel_PTENwt", sep=""))
#write.csv(GSEA_merge, 
#          file="/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_SPOP_CHD1_ERG_PTEN_wt_hallmark_merge.csv", row.names = F)
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_SPOP_CHD1_ERG_PTEN_wt_hallmark-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()



#GSEA c6 signatures
#Figure S3: Correlation c6 GSEA_T_SPOPmut_CHD1wt_N and GSEA_T_ERG_PTENwt_N ggplot 
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=14, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1wt_N", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENwt_N", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP_CHD1_ERG_PTEN_N_c6-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1wt_N, y=NES.ERG.PTENwt_N)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1wt_N+0.000000001), colour=-log10(NOM.p.val.ERG.PTENwt_N+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()


#Figure 2C: Correlation c6 SPOP_CHD1del/CHD1wt, ERG_PTENdel/PTENwt ggplot in TCGA ggplot
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=15, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP_CHD1_ERG_PTEN_wt_c6-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()



#Figure 2C: Correlation c6 SPOP_CHD1del/CHD1wt, ERG_PTENdel/PTENwt ggplot in Taylor ggplot 
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=16, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_SPOP_CHD1_ERG_PTEN_wt_c6-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()


#Figure S5: Correlation c6 SPOP_CHD1del/CHD1wt, ERG_PTENdel/PTENwt ggplot in ICGC ggplot 
GSEA_merge <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=17, check.names=F)
colnames(GSEA_merge) <- c("NAME", 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".SPOP.CHD1del_CHD1wt", sep=""), 
                          paste(c("NES", "NOM.p.val", "FDR.q.val"), ".ERG.PTENdel_PTENwt", sep=""))
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_SPOP_CHD1_ERG_PTEN_wt_c6-cor_v2.pdf", width=9, height=5, useDingbats=FALSE)
ggplot(GSEA_merge, aes(x= NES.SPOP.CHD1del_CHD1wt, y=NES.ERG.PTENdel_PTENwt)) +
  geom_point(aes(size=-log10(NOM.p.val.SPOP.CHD1del_CHD1wt+0.000000001), colour=-log10(NOM.p.val.ERG.PTENdel_PTENwt+0.000000001)), shape=19) +
  scale_color_gradient(low=brewer.pal(9,"Greys")[3], high=brewer.pal(9,"Greys")[8])+
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE)
dev.off()




#GSEA signatures: hallmark in TCGA
GSEA_T_SPOPmut_CHD1del_CHD1wt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_SPOPmut_CHD1wt_CHD1del_rowSum_10_DESeq2.rnk_hallmark.GseaPreranked.1594351043940/gsea_report_for_na_pos_1594351043940.xls", sep="\t", header=T)
GSEA_T_SPOPmut_CHD1del_CHD1wt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_SPOPmut_CHD1wt_CHD1del_rowSum_10_DESeq2.rnk_hallmark.GseaPreranked.1594351043940/gsea_report_for_na_neg_1594351043940.xls", sep="\t", header=T)
GSEA_T_ERG_PTENdel_PTENwt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_ERG_PTENwt_PTENdel_rowSum_10_DESeq2.rnk_hallmark.GseaPreranked.1594352172838/gsea_report_for_na_pos_1594352172838.xls", sep="\t", header=T)
GSEA_T_ERG_PTENdel_PTENwt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_ERG_PTENwt_PTENdel_rowSum_10_DESeq2.rnk_hallmark.GseaPreranked.1594352172838/gsea_report_for_na_neg_1594352172838.xls", sep="\t", header=T)
GSEA_T_SPOPmut_CHD1del_CHD1wt <- rbind(GSEA_T_SPOPmut_CHD1del_CHD1wt_pos, GSEA_T_SPOPmut_CHD1del_CHD1wt_neg)
GSEA_T_ERG_PTENdel_PTENwt <- rbind(GSEA_T_ERG_PTENdel_PTENwt_pos, GSEA_T_ERG_PTENdel_PTENwt_neg)
GSEA_T_SPOPmut_CHD1del_CHD1wt$Input<-"T_SPOPmut_CHD1del/CHD1wt"
GSEA_T_ERG_PTENdel_PTENwt$Input<-"T_ERG_PTENdel/PTENwt"
#GSEA signatures: hallmark in ICGC
GSEA_ICGC_ERG_PTENdel_PTENwt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_ERG_PTENdel_ERG_PTENwt_DESeq2.rnk_hallmark.GseaPreranked.1595280218132/gsea_report_for_na_pos_1595280218132.xls", sep="\t", header=T)
GSEA_ICGC_ERG_PTENdel_PTENwt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_ERG_PTENdel_ERG_PTENwt_DESeq2.rnk_hallmark.GseaPreranked.1595280218132/gsea_report_for_na_neg_1595280218132.xls", sep="\t", header=T)
GSEA_ICGC_SPOP_CHD1del_CHD1wt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_SPOP_CHD1del_SPOP_CHD1wt_DESeq2.rnk_hallmark.GseaPreranked.1595280304181/gsea_report_for_na_pos_1595280304181.xls", sep="\t", header=T)
GSEA_ICGC_SPOP_CHD1del_CHD1wt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_SPOP_CHD1del_SPOP_CHD1wt_DESeq2.rnk_hallmark.GseaPreranked.1595280304181/gsea_report_for_na_neg_1595280304181.xls", sep="\t", header=T)
GSEA_ICGC_ERG_PTENdel_PTENwt <- rbind(GSEA_ICGC_ERG_PTENdel_PTENwt_pos, GSEA_ICGC_ERG_PTENdel_PTENwt_neg)
GSEA_ICGC_SPOP_CHD1del_CHD1wt <- rbind(GSEA_ICGC_SPOP_CHD1del_CHD1wt_pos, GSEA_ICGC_SPOP_CHD1del_CHD1wt_neg)
GSEA_ICGC_ERG_PTENdel_PTENwt$Input<-"ICGC_ERG.PTENdel/PTENwt"
GSEA_ICGC_SPOP_CHD1del_CHD1wt$Input<-"ICGC_SPOP.CHD1del/CHD1wt"
#GSEA signatures: hallmark in Taylor
GSEA_Taylor_ERG_PTENdel_PTENwt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_ERG_PTENdel_PTENwt__log2FC.rnk_hallmark.GseaPreranked.1597892651196/gsea_report_for_na_pos_1597892651196.tsv", sep="\t", header=T)
GSEA_Taylor_ERG_PTENdel_PTENwt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_ERG_PTENdel_PTENwt__log2FC.rnk_hallmark.GseaPreranked.1597892651196/gsea_report_for_na_neg_1597892651196.tsv", sep="\t", header=T)
GSEA_Taylor_SPOP_CHD1del_CHD1wt_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_SPOP_CHD1del_CHD1wt__log2FC.rnk_hallmark.GseaPreranked.1597892580098/gsea_report_for_na_pos_1597892580098.tsv", sep="\t", header=T)
GSEA_Taylor_SPOP_CHD1del_CHD1wt_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_SPOP_CHD1del_CHD1wt__log2FC.rnk_hallmark.GseaPreranked.1597892580098/gsea_report_for_na_neg_1597892580098.tsv", sep="\t", header=T)
GSEA_Taylor_ERG_PTENdel_PTENwt <- rbind(GSEA_Taylor_ERG_PTENdel_PTENwt_pos, GSEA_Taylor_ERG_PTENdel_PTENwt_neg)
GSEA_Taylor_SPOP_CHD1del_CHD1wt <- rbind(GSEA_Taylor_SPOP_CHD1del_CHD1wt_pos, GSEA_Taylor_SPOP_CHD1del_CHD1wt_neg)
GSEA_Taylor_ERG_PTENdel_PTENwt$Input<-"Taylor_ERG.PTENdel/PTENwt"
GSEA_Taylor_SPOP_CHD1del_CHD1wt$Input<-"Taylor_SPOP.CHD1del/CHD1wt"

#GSEA output comparison between CHD1_mouse and ERG+_PTENloss_mouse
#Hallmark signatures
GSEA_mouse_CHD1_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Michael4038_2016_09_22-RNAseq-DESeq_results-human.csv-log2FC.rnk_hallmark.GseaPreranked.1574871820637/gsea_report_for_na_pos_1574871820637.xls", sep="\t", header=T)
GSEA_mouse_CHD1_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/Michael4038_2016_09_22-RNAseq-DESeq_results-human.csv-log2FC.rnk_hallmark.GseaPreranked.1574871820637/gsea_report_for_na_neg_1574871820637.xls", sep="\t", header=T)
GSEA_mouse_ERG_PTEN_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/GSE46799_WT_ERG_PTEN_log2FC.rnk_hallmark.GseaPreranked.1594655288718/gsea_report_for_na_pos_1594655288718.xls", sep="\t", header=T)
GSEA_mouse_ERG_PTEN_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/GSE46799_WT_ERG_PTEN_log2FC.rnk_hallmark.GseaPreranked.1594655288718/gsea_report_for_na_neg_1594655288718.xls", sep="\t", header=T)
GSEA_mouse_PTEN_pos<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/GSE46799_WT_PTEN_log2FC.rnk_hallmark.GseaPreranked.1594655359196/gsea_report_for_na_pos_1594655359196.xls", sep="\t", header=T)
GSEA_mouse_PTEN_neg<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/GSE46799_WT_PTEN_log2FC.rnk_hallmark.GseaPreranked.1594655359196/gsea_report_for_na_neg_1594655359196.xls", sep="\t", header=T)
GSEA_mouse_CHD1<-rbind(GSEA_mouse_CHD1_pos, GSEA_mouse_CHD1_neg[nrow(GSEA_mouse_CHD1_neg):1,])
GSEA_mouse_ERG_PTEN<-rbind(GSEA_mouse_ERG_PTEN_pos, GSEA_mouse_ERG_PTEN_neg[nrow(GSEA_mouse_ERG_PTEN_neg):1,])
GSEA_mouse_PTEN<-rbind(GSEA_mouse_PTEN_pos, GSEA_mouse_PTEN_neg)
GSEA_mouse_CHD1$Input<-"Mouse_CHD1/wt"
GSEA_mouse_ERG_PTEN$Input<-"Mouse_ERG.PTEN/wt"
GSEA_mouse_PTEN$Input<-"Mouse_PTEN/wt"
#GSEA_merge_mouse<-rbind(GSEA_mouse_CHD1[,c(1,6,7,8,13)], GSEA_mouse_ERG_PTEN[,c(1,6,7,8,13)])
#Combind human and mouse GSEA output
GSEA_merge_mouse <- rbind(GSEA_mouse_CHD1[,c(1,6,7,8,13)], GSEA_mouse_ERG_PTEN[,c(1,6,7,8,13)])

#Heatmap of human TCGA, ICGA, Taylor and mouse GSEA hallmark signatures
myPalette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
GSEA_merge_TCGA <- rbind(#GSEA_T_SPOPmut_CHD1wt_N[,c(1,6,7,8,13)], 
                         GSEA_T_SPOPmut_CHD1del_CHD1wt[,c(1,6,7,8,13)], 
                         #GSEA_T_ERG_PTENwt_N[,c(1,6,7,8,13)], 
                         GSEA_T_ERG_PTENdel_PTENwt[,c(1,6,7,8,13)])
GSEA_merge_ICGC <- rbind(GSEA_ICGC_ERG_PTENdel_PTENwt[,c(1,6,7,8,13)], 
                         GSEA_ICGC_SPOP_CHD1del_CHD1wt[,c(1,6,7,8,13)])
GSEA_merge_Taylor <- rbind(GSEA_Taylor_ERG_PTENdel_PTENwt[,c(1,6,7,8,13)], 
                         GSEA_Taylor_SPOP_CHD1del_CHD1wt[,c(1,6,7,8,13)])
GSEA_merge_mouse <- rbind(GSEA_mouse_CHD1[,c(1,6,7,8,13)], GSEA_mouse_ERG_PTEN[,c(1,6,7,8,13)])
#GSEA_merge_mouse <- rbind(GSEA_mouse_CHD1[,c(1,6,7,8,13)], GSEA_mouse_PTEN[,c(1,6,7,8,13)])
GSEA_merge_combine <- rbind(GSEA_merge_TCGA[, c(1,2,5)], 
                            GSEA_merge_ICGC[, c(1,2,5)], 
                            GSEA_merge_Taylor[, c(1,2,5)],
                            GSEA_merge_mouse[, c(1,2,5)])
table(GSEA_merge_combine$Input)
GSEA_merge_combine$NAME <- sub("_OXIGEN_", "_OXYGEN_", GSEA_merge_combine$NAME)
GSEA_merge_combine$NAME_v2 <- sub("HALLMARK_", "", GSEA_merge_combine$NAME)
GSEA_merge_combine$NAME_v3 <- lapply(GSEA_merge_combine$NAME_v2, tolower) 
GSEA_merge_combine$Input_v2 <- factor(GSEA_merge_combine$Input, 
                                      levels = c("T_ERG_PTENdel/PTENwt", "Taylor_ERG.PTENdel/PTENwt", "ICGC_ERG.PTENdel/PTENwt", "Mouse_ERG.PTEN/wt", 
                                                 "T_SPOPmut_CHD1del/CHD1wt", "Taylor_SPOP.CHD1del/CHD1wt", "ICGC_SPOP.CHD1del/CHD1wt", "Mouse_CHD1/wt"))
#data <- GSEA_merge_combine
pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Figures/Human_TCGA_ICGC_Taylor_mouse_hallmark-hm.pdf", width=12, height=14)
ggplot(data=GSEA_merge_combine, aes(y=NAME_v2, x=Input_v2, fill=NES)) +
  geom_tile(color="white") +
  coord_equal() + 
  scale_fill_gradientn(colours = myPalette(100), limits=c(min(data$NES), max(data$NES))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



