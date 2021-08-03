library(survival)
library(survminer)

#write.csv(svm.pred.info_clas_cohort, 
#          file="/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher.csv", row.names=F)

nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",])

#MET-free survival
#SPOP_signature and CHD1_signature
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-CHD1_SPOP-SUR_MET.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos" | svm.pred.info_clas_cohort$CHD1_signature=="Pos",]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$CHD1_signature=="Pos", "CHD1_signature", "SPOP_signature")
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(met_time, met) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_CHD1o, sur_col_SPOPo),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="MET time (month)")
dev.off()

#ERG_PTENwt and ERG_PTENdel
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-ERG_PTEN-SUR_MET.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENwt=="Pos" | svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$ERG_PTENdel=="Pos", "ERG_PTENdel", "ERG_PTENwt")
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(met_time, met) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_ERG_PTENdel, sur_col_ERG_PTENwt),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="MET time (month)")
dev.off()

#SPOP_signature and CHD1_signature + ERG_PTENwt and ERG_PTENdel
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-CHD1_SPOP_ERG_PTEN-SUR_MET.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos" | svm.pred.info_clas_cohort$CHD1_signature=="Pos" | 
                                                           svm.pred.info_clas_cohort$ERG_PTENwt=="Pos" | svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$CHD1_signature=="Pos", "CHD1_signature", 
                                              ifelse(svm.pred.info_clas_cohort_sur$SPOP_signature=="Pos", "SPOP_signature",
                                                     ifelse(svm.pred.info_clas_cohort_sur$ERG_PTENdel=="Pos", "ERG_PTENdel", "ERG_PTENwt")))
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(met_time, met) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_CHD1o, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_SPOPo),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="MET time (month)")
dev.off()



#BCR-free survival
#SPOP_signature and CHD1_signature + ERG_PTENwt and ERG_PTENdel
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-CHD1_SPOP_ERG_PTEN-SUR_BCR.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos" | svm.pred.info_clas_cohort$CHD1_signature=="Pos" | 
                                                           svm.pred.info_clas_cohort$ERG_PTENwt=="Pos" | svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$CHD1_signature=="Pos", "CHD1_signature", 
                                              ifelse(svm.pred.info_clas_cohort_sur$SPOP_signature=="Pos", "SPOP_signature",
                                                     ifelse(svm.pred.info_clas_cohort_sur$ERG_PTENdel=="Pos", "ERG_PTENdel", "ERG_PTENwt")))
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(bcr_time, bcr) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_CHD1o, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_SPOPo),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="BCR time (month)")
dev.off()


#PCSM-free survival
#SPOP_signature and CHD1_signature + ERG_PTENwt and ERG_PTENdel
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-CHD1_SPOP_ERG_PTEN-SUR_PCSM.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos" | svm.pred.info_clas_cohort$CHD1_signature=="Pos" | 
                                                           svm.pred.info_clas_cohort$ERG_PTENwt=="Pos" | svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$CHD1_signature=="Pos", "CHD1_signature", 
                                              ifelse(svm.pred.info_clas_cohort_sur$SPOP_signature=="Pos", "SPOP_signature",
                                                     ifelse(svm.pred.info_clas_cohort_sur$ERG_PTENdel=="Pos", "ERG_PTENdel", "ERG_PTENwt")))
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(pcsm_time, pcsm) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_CHD1o, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_SPOPo),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="PCSM time (month)")
dev.off()



