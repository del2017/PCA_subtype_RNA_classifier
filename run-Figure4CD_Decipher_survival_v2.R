library(survival)
library(survminer)


#Call the ERG/ETS first, then SPOP/CHD1/PTEN calling
svm.pred.info_clas_cohort$ERG_status<-ifelse(svm.pred.info_clas_cohort$ERG_overexp=="Pos" & svm.pred.info_clas_cohort$ETS_overexp=="Neg", "Pos", "Neg")
svm.pred.info_clas_cohort$ETS_status<-ifelse(svm.pred.info_clas_cohort$ETS_overexp=="Pos" & svm.pred.info_clas_cohort$ERG_overexp=="Neg", "Pos", "Neg")

svm.pred.info_clas_cohort$SPOPmut<-ifelse(svm.pred.info_clas_cohort$pred_SPOP==1 
                                          & (svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$ETS_status=="Neg")
                                          & svm.pred.info_clas_cohort$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1del<-ifelse(svm.pred.info_clas_cohort$pred_CHD1==1 
                                          & (svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$ETS_status=="Neg")
                                          & svm.pred.info_clas_cohort$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1_SPOP_share<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Pos" & svm.pred.info_clas_cohort$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1_signature<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Neg" & svm.pred.info_clas_cohort$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_cohort$SPOP_signature<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Pos" & svm.pred.info_clas_cohort$CHD1del=="Neg", "Pos", "Neg")

svm.pred.info_clas_cohort$ERG_PTEN_share <- ifelse(svm.pred.info_clas_cohort$ERG_status=="Pos" & svm.pred.info_clas_cohort$pred_PTEN==1
                                                   & svm.pred.info_clas_cohort$pred_CHD1==0 & svm.pred.info_clas_cohort$pred_SPOP==0, "Pos", "Neg")
svm.pred.info_clas_cohort$ERG_PTENwt <- ifelse(svm.pred.info_clas_cohort$ERG_status=="Pos" & svm.pred.info_clas_cohort$pred_PTEN==0
                                               & svm.pred.info_clas_cohort$pred_CHD1==0 & svm.pred.info_clas_cohort$pred_SPOP==0, "Pos", "Neg")
svm.pred.info_clas_cohort$ERGneg_PTENdel <- ifelse(svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$pred_PTEN==1
                                                   & svm.pred.info_clas_cohort$ETS_status=="Neg" 
                                                   & svm.pred.info_clas_cohort$pred_CHD1==0 & svm.pred.info_clas_cohort$pred_SPOP==0, "Pos", "Neg")


#MET-free survival
svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos" 
                                                         | svm.pred.info_clas_cohort$CHD1_signature=="Pos" 
                                                         | svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos" 
                                                         | svm.pred.info_clas_cohort$ERG_PTEN_share=="Pos"
                                                         | svm.pred.info_clas_cohort$ERG_PTENwt=="Pos"
                                                         | svm.pred.info_clas_cohort$ERGneg_PTENdel=="Pos", ]
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$CHD1_signature=="Pos", "CHD1_signature", 
                                              ifelse(svm.pred.info_clas_cohort_sur$SPOP_signature=="Pos", "SPOP_signature",
                                                     ifelse(svm.pred.info_clas_cohort_sur$CHD1_SPOP_share=="Pos", "SPOP_CHD1_share",
                                                            ifelse(svm.pred.info_clas_cohort_sur$ERG_PTEN_share=="Pos", "ERG_PTEN_share", 
                                                                   ifelse(svm.pred.info_clas_cohort_sur$ERG_PTENwt=="Pos", "ERG_PTENwt", 
                                                                          ifelse(svm.pred.info_clas_cohort_sur$ERGneg_PTENdel=="Pos", "ERGneg_PTENdel", "NA"))))))
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(pcsm_time, pcsm) ~ Subtype, data = svm.pred.info_clas_cohort_sur)

pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Figures/Decipher_retro_svm-clas-CHD1_SPOP_ERG_PTEN_6groups-SUR_MET.pdf", width=6, height=6, onefile=F)
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(met_time, met) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_CHD1o, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ERGneg_PTENdel, sur_col_CHD1_SPOP, sur_col_SPOPo),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="MET time (month)")
dev.off()




