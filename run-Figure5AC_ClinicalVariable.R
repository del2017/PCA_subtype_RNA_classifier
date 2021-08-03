library(MASS)
library(forestplot)
library(alluvial)
library(tidyr)
library(dplyr)

#Subtype prediction from Figure 4
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",])


#Figure 5A: Clinical and pathological characteristics between PTEN and CHD1 deleted tumors in Decipher retrospective cohort (n=1,626) via univariable analyses, with Other samples as reference.

#Clinical outcome comparison between subgroups CHD1 and SPOP
#exclude NA in met_time 
svm.pred.info_clas_cohortm<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos" | svm.pred.info_clas_cohort$SPOP_signature=="Pos",]
svm.pred.info_clas_cohortm$CHD1_del<-ifelse(svm.pred.info_clas_cohortm$CHD1del=="Pos", 1,0)

#CHD1_signature and EGR_PTENdel_signature clinical outcomes, by using reference from the rest of samples
svm.pred.info_clas_cohortm<-svm.pred.info_clas_cohort
#Generate subtype reference
svm.pred.info_clas_cohortm$Subtype<-ifelse(svm.pred.info_clas_cohortm$CHD1_signature=="Pos", 1, 
                                           ifelse(svm.pred.info_clas_cohortm$ERG_PTENdel=="Pos", 2, 3))
svm.pred.info_clas_cohortm$Subtype<-as.factor(svm.pred.info_clas_cohortm$Subtype)
table(svm.pred.info_clas_cohortm$Subtype)
#Univariate analysis 
#GS<7 (GS=7 as ref) 
svm.pred.info_clas_cohortm_GS7less<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7less$pathgs))<7 | 
                                                     (as.numeric(as.character(svm.pred.info_clas_cohortm_GS7less$pathgs_p))==3 & as.numeric(as.character(svm.pred.info_clas_cohortm_GS7less$pathgs_p))==4), 1, 0)
fit <- glm(GS7less ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7
svm.pred.info_clas_cohortm_GS7more<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7more$pathgs))>=8 | 
                                                     (as.numeric(as.character(svm.pred.info_clas_cohortm_GS7more$pathgs_p))==4 & as.numeric(as.character(svm.pred.info_clas_cohortm_GS7more$pathgs_p))==3), 1, 0)
fit <- glm(GS7more ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 or GS=4+3 vs. GS<7 or GS=3+4 (GS<7 or GS=3+4 as ref) vs. ERG_PTEN
svm.pred.info_clas_cohortm_GS7m_l<-svm.pred.info_clas_cohortm[complete.cases(svm.pred.info_clas_cohortm$pathgs), ]
svm.pred.info_clas_cohortm_GS7m_l$GS7_4_3_more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7m_l$pathgs))>7
                                                       | (svm.pred.info_clas_cohortm_GS7m_l$pathgs_p==4 & svm.pred.info_clas_cohortm_GS7m_l$pathgs_s==3), 1, 0)
fit <- glm(GS7_4_3_more ~ C(Subtype, base=3),  data=svm.pred.info_clas_cohortm_GS7m_l, family = "binomial")
GS7_4_3_more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_cohortm$lni_status<-ifelse(svm.pred.info_clas_cohortm$lni==1, 1, 0)
fit <- glm(lni_status ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_cohortm$sm_status<-ifelse(svm.pred.info_clas_cohortm$sm==1, 1, 0)
fit <- glm(sm_status ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#ece
svm.pred.info_clas_cohortm$ece_status<-ifelse(svm.pred.info_clas_cohortm$ece==1, 1, 0)
fit <- glm(ece_status ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_cohortm$svi_status<-ifelse(svm.pred.info_clas_cohortm$svi==1, 1, 0)
fit <- glm(svi_status ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_cohortm$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_cohortm$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#race
svm.pred.info_clas_cohortm$race_BAA<-ifelse(svm.pred.info_clas_cohortm$race=="African American", 1, 0)
fit <- glm(race_BAA ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_race<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_cohortm$age70m<-ifelse(svm.pred.info_clas_cohortm$age>=70, 1, 0)
fit <- glm(age70m ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#adt
fit <- glm(adt ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_adt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#rt
fit <- glm(rt ~ C(Subtype, base=3), data=svm.pred.info_clas_cohortm, family = "binomial")
m_rt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]


table(svm.pred.info_clas_cohortm$Subtype)
results<-NULL
for (i in 1:2)
{
  results<-cbind(results, 
                 rbind(cbind(paste(formatC(m_gs7less[i,1],digits=3,format = "f"), "(", formatC(m_gs7less[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_gs7less[i,3],digits=3,format = "f"),")", sep=""),formatC(m_gs7less[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_gs7more[i,1],digits=3,format = "f"), "(", formatC(m_gs7more[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_gs7more[i,3],digits=3,format = "f"),")", sep=""),formatC(m_gs7more[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_lni[i,1],digits=3,format = "f"), "(", formatC(m_lni[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_lni[i,3],digits=3,format = "f"),")", sep=""),formatC(m_lni[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_sm[i,1],digits=3,format = "f"), "(", formatC(m_sm[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_sm[i,3],digits=3,format = "f"),")", sep=""),formatC(m_sm[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_ece[i,1],digits=3,format = "f"), "(", formatC(m_ece[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_ece[i,3],digits=3,format = "f"),")", sep=""),formatC(m_ece[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_svi[i,1],digits=3,format = "f"), "(", formatC(m_svi[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_svi[i,3],digits=3,format = "f"),")", sep=""),formatC(m_svi[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_psa10l[i,1],digits=3,format = "f"), "(", formatC(m_psa10l[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_psa10l[i,3],digits=3,format = "f"),")", sep=""),formatC(m_psa10l[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_psa20m[i,1],digits=3,format = "f"), "(", formatC(m_psa20m[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_psa20m[i,3],digits=3,format = "f"),")", sep=""),formatC(m_psa20m[i,4],digits=3,format = "f")),
                       #cbind(paste(formatC(m_race[i,1],digits=3,format = "f"), "(", formatC(m_race[i,2],digits=3,format = "f"), ",", 
                       #             formatC(m_race[i,3],digits=3,format = "f"),")", sep=""),formatC(m_race[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_age[i,1],digits=3,format = "f"), "(", formatC(m_age[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_age[i,3],digits=3,format = "f"),")", sep=""),formatC(m_age[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_adt[i,1],digits=3,format = "f"), "(", formatC(m_adt[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_adt[i,3],digits=3,format = "f"),")", sep=""),formatC(m_adt[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_rt[i,1],digits=3,format = "f"), "(", formatC(m_rt[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_rt[i,3],digits=3,format = "f"),")", sep=""),formatC(m_rt[i,4],digits=3,format = "f"))
                 )
  )
}
rownames(results)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                     "Seminal vesicle invasion", "PSA<10","PSA>20","Age>=70","ADT","RT")
results
#write.csv(results, file="/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-CHD1_ERG_PTENdel-UV-forestplot.csv", row.names=T)


m_CHD1_sig_ref_other<-data.frame(rbind(m_gs7less[1,], m_gs7more[1,], m_lni[1,], m_sm[1,], m_ece[1,], m_svi[1,], 
                                       m_psa10l[1,], m_psa20m[1,], m_age[1,], m_adt[1,], m_rt[1,]))
row.names(m_CHD1_sig_ref_other)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                   "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70","ADT","RT")
m_ERG_PTENdel_sig_ref_other<-data.frame(rbind(m_gs7less[2,], m_gs7more[2,], m_lni[2,], m_sm[2,], m_ece[2,], m_svi[2,], 
                                              m_psa10l[2,], m_psa20m[2,], m_age[2,], m_adt[2,], m_rt[2,]))
row.names(m_ERG_PTENdel_sig_ref_other)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                          "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70","ADT","RT")
#forest plot of odd ratio
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-CHD1_other-UV-forestplot.pdf", height=8, width=6)
forestplot(m_CHD1_sig_ref_other[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_CHD1_sig_ref_other[,4])/7.5,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-ERG_PTENdel_other-UV-forestplot.pdf", height=8, width=6)
forestplot(m_ERG_PTENdel_sig_ref_other[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_ERG_PTENdel_sig_ref_other[,4])/7.5,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()



#Figure 5C: Clinical and pathological characteristics between PTEN and CHD1 deleted status in Decipher prospective cohort (n=6,532) via univariable analyses, with other samples as reference. 
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ETS_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$other=="Pos",])

#CHD1_signature and EGR_PTENdel_signature clinical outcomes, by using reference from the rest of samples
svm.pred.info_clas_6532m<-svm.pred.info_clas_6532
#Generate subtype reference
svm.pred.info_clas_6532m$Subtype<-ifelse(svm.pred.info_clas_6532m$CHD1_signature=="Pos", 1, 
                                         ifelse(svm.pred.info_clas_6532m$ERG_PTENdel=="Pos", 2, 3))
svm.pred.info_clas_6532m$Subtype<-as.factor(svm.pred.info_clas_6532m$Subtype)
table(svm.pred.info_clas_6532m$Subtype)
#Univariate analysis 
#GS<7 (GS=7 as ref) 
svm.pred.info_clas_6532m_GS7less<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7less$pathgs))<7 | 
                                                   (as.numeric(as.character(svm.pred.info_clas_6532m_GS7less$pathgs_p))==3 & as.numeric(as.character(svm.pred.info_clas_6532m_GS7less$pathgs_p))==4), 1, 0)
fit <- glm(GS7less ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7
svm.pred.info_clas_6532m_GS7more<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7more$pathgs))>=8 | 
                                                   (as.numeric(as.character(svm.pred.info_clas_6532m_GS7more$pathgs_p))==4 & as.numeric(as.character(svm.pred.info_clas_6532m_GS7more$pathgs_p))==3), 1, 0)
fit <- glm(GS7more ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_6532m$lni_status<-ifelse(svm.pred.info_clas_6532m$lni==1, 1, 0)
fit <- glm(lni_status ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_6532m$sm_status<-ifelse(svm.pred.info_clas_6532m$sm==1, 1, 0)
fit <- glm(sm_status ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#ece
svm.pred.info_clas_6532m$ece_status<-ifelse(svm.pred.info_clas_6532m$epe==1, 1, 0)
fit <- glm(ece_status ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_6532m$svi_status<-ifelse(svm.pred.info_clas_6532m$svi==1, 1, 0)
fit <- glm(svi_status ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_6532m$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_6532m$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_6532m$age70m<-ifelse(svm.pred.info_clas_6532m$age>=70, 1, 0)
fit <- glm(age70m ~ C(Subtype, base=3), data=svm.pred.info_clas_6532m, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]


table(svm.pred.info_clas_6532m$Subtype)
results<-NULL
for (i in 1:2)
{
  results<-cbind(results, 
                 rbind(cbind(paste(formatC(m_gs7less[i,1],digits=3,format = "f"), "(", formatC(m_gs7less[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_gs7less[i,3],digits=3,format = "f"),")", sep=""),formatC(m_gs7less[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_gs7more[i,1],digits=3,format = "f"), "(", formatC(m_gs7more[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_gs7more[i,3],digits=3,format = "f"),")", sep=""),formatC(m_gs7more[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_lni[i,1],digits=3,format = "f"), "(", formatC(m_lni[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_lni[i,3],digits=3,format = "f"),")", sep=""),formatC(m_lni[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_sm[i,1],digits=3,format = "f"), "(", formatC(m_sm[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_sm[i,3],digits=3,format = "f"),")", sep=""),formatC(m_sm[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_ece[i,1],digits=3,format = "f"), "(", formatC(m_ece[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_ece[i,3],digits=3,format = "f"),")", sep=""),formatC(m_ece[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_svi[i,1],digits=3,format = "f"), "(", formatC(m_svi[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_svi[i,3],digits=3,format = "f"),")", sep=""),formatC(m_svi[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_psa10l[i,1],digits=3,format = "f"), "(", formatC(m_psa10l[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_psa10l[i,3],digits=3,format = "f"),")", sep=""),formatC(m_psa10l[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_psa20m[i,1],digits=3,format = "f"), "(", formatC(m_psa20m[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_psa20m[i,3],digits=3,format = "f"),")", sep=""),formatC(m_psa20m[i,4],digits=3,format = "f")),
                       cbind(paste(formatC(m_age[i,1],digits=3,format = "f"), "(", formatC(m_age[i,2],digits=3,format = "f"), ",", 
                                   formatC(m_age[i,3],digits=3,format = "f"),")", sep=""),formatC(m_age[i,4],digits=3,format = "f"))
                 )
  )
}
rownames(results)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                     "Seminal vesicle invasion", "PSA<10","PSA>20","Age>=70")
results
#write.csv(results, file="/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-CHD1_ERG_PTENdel-UV-forestplot.csv", row.names=T)

m_CHD1_sig_ref_other<-data.frame(rbind(m_gs7less[1,], m_gs7more[1,], m_lni[1,], m_sm[1,], m_ece[1,], m_svi[1,], 
                                       m_psa10l[1,], m_psa20m[1,], m_age[1,]))
row.names(m_CHD1_sig_ref_other)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                   "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70")
m_ERG_PTENdel_sig_ref_other<-data.frame(rbind(m_gs7less[2,], m_gs7more[2,], m_lni[2,], m_sm[2,], m_ece[2,], m_svi[2,], 
                                              m_psa10l[2,], m_psa20m[2,], m_age[2,]))
row.names(m_ERG_PTENdel_sig_ref_other)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                          "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70")
#forest plot of odd ratio
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-CHD1_other-UV-forestplot.pdf", height=8, width=6)
forestplot(m_CHD1_sig_ref_other[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_CHD1_sig_ref_other[,4])/15,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-ERG_PTENdel_other-UV-forestplot.pdf", height=8, width=6)
forestplot(m_ERG_PTENdel_sig_ref_other[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_ERG_PTENdel_sig_ref_other[,4])/15,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()



#Figure S18. Pathological outcome difference between progression and early event from each tumor lineage model in Decipher retrospective cohort (n=1,626) via univariable analyses.

#Univariate analysis of between CHD1del and SPOPmut (reference: SPOP)
#GS<7 (GS=7 as ref) vs. CHD1del
svm.pred.info_clas_cohortm_GS7less<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7less$pathgs))<7, 1, 0)
fit <- glm(GS7less ~ CHD1_del, data=svm.pred.info_clas_cohortm_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 (GS=7 as ref) vs. CHD1del
svm.pred.info_clas_cohortm_GS7more<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7more$pathgs))>7, 1, 0)
fit <- glm(GS7more ~ CHD1_del, data=svm.pred.info_clas_cohortm_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 or GS=4+3 vs. GS<7 or GS=3+4 (GS<7 or GS=3+4 as ref) vs. CHD1del
svm.pred.info_clas_cohortm_GS7m_l<-svm.pred.info_clas_cohortm[complete.cases(svm.pred.info_clas_cohortm$pathgs), ]
svm.pred.info_clas_cohortm_GS7m_l$GS7_4_3_more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7m_l$pathgs))>7
                                                       | (svm.pred.info_clas_cohortm_GS7m_l$pathgs_p==4 & svm.pred.info_clas_cohortm_GS7m_l$pathgs_s==3), 1, 0)
fit <- glm(GS7_4_3_more ~ CHD1_del, data=svm.pred.info_clas_cohortm_GS7m_l, family = "binomial")
GS7_4_3_more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_cohortm$lni_status<-ifelse(svm.pred.info_clas_cohortm$lni==1, 1, 0)
fit <- glm(lni_status ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_cohortm$sm_status<-ifelse(svm.pred.info_clas_cohortm$sm==1, 1, 0)
fit <- glm(sm_status ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#ece
svm.pred.info_clas_cohortm$ece_status<-ifelse(svm.pred.info_clas_cohortm$ece==1, 1, 0)
fit <- glm(ece_status ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_cohortm$svi_status<-ifelse(svm.pred.info_clas_cohortm$svi==1, 1, 0)
fit <- glm(svi_status ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_cohortm$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_cohortm$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_cohortm$age70m<-ifelse(svm.pred.info_clas_cohortm$age>=70, 1, 0)
fit <- glm(age70m ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#adt
fit <- glm(adt ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_adt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#rt
fit <- glm(rt ~ CHD1_del, data=svm.pred.info_clas_cohortm, family = "binomial")
m_rt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]

m_all_refSPOP<-data.frame(rbind(m_gs7less, m_gs7more, GS7_4_3_more, m_lni, m_sm, m_ece, m_svi, 
                                m_psa10l, m_psa20m, m_age))#, m_adt, m_rt))
row.names(m_all_refSPOP)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "GS>7(4+3) (GS<7(3+4) as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                            "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70") #,"ADT","RT")
m_all_refSPOP
#write.csv(m_all_refSPOP, file="/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-UV-CHD1_SPOP-forestplot.csv", row.names=T)
#forest plot of odd ratio
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-UV-CHD1_SPOP-forestplot.pdf", height=8, width=6)
forestplot(m_all_refSPOP[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_all_refSPOP[,4])/5,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()



#Clinical outcome comparison between subgroups ERG_PTENdel and ERG_PTENwt
#exclude NA in met_time 
svm.pred.info_clas_cohortm<-svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENdel=="Pos" | svm.pred.info_clas_cohort$ERG_PTENwt=="Pos",]
svm.pred.info_clas_cohortm$ERG_PTENdel_status<-ifelse(svm.pred.info_clas_cohortm$ERG_PTENdel=="Pos", 1,0)
#Univariate analysis of between ERG_PTENdel and ERG_PTENwt (reference: ERG_PTENwt)
#GS<7 (GS=7 as ref) vs. ERG_PTEN
svm.pred.info_clas_cohortm_GS7less<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7less$pathgs))<7, 1, 0)
fit <- glm(GS7less ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 (GS=7 as ref) vs. ERG_PTEN
svm.pred.info_clas_cohortm_GS7more<-svm.pred.info_clas_cohortm[as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_cohortm$pathgs))==7, ]
svm.pred.info_clas_cohortm_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7more$pathgs))>7, 1, 0)
fit <- glm(GS7more ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 or GS=4+3 vs. GS<7 or GS=3+4 (GS<7 or GS=3+4 as ref) vs. ERG_PTEN
svm.pred.info_clas_cohortm_GS7m_l<-svm.pred.info_clas_cohortm[complete.cases(svm.pred.info_clas_cohortm$pathgs), ]
svm.pred.info_clas_cohortm_GS7m_l$GS7_4_3_more<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm_GS7m_l$pathgs))>7
                                                       | (svm.pred.info_clas_cohortm_GS7m_l$pathgs_p==4 & svm.pred.info_clas_cohortm_GS7m_l$pathgs_s==3), 1, 0)
fit <- glm(GS7_4_3_more ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm_GS7m_l, family = "binomial")
GS7_4_3_more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_cohortm$lni_status<-ifelse(svm.pred.info_clas_cohortm$lni==1, 1, 0)
fit <- glm(lni_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_cohortm$sm_status<-ifelse(svm.pred.info_clas_cohortm$sm==1, 1, 0)
fit <- glm(sm_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#ece
svm.pred.info_clas_cohortm$ece_status<-ifelse(svm.pred.info_clas_cohortm$ece==1, 1, 0)
fit <- glm(ece_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_cohortm$svi_status<-ifelse(svm.pred.info_clas_cohortm$svi==1, 1, 0)
fit <- glm(svi_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_cohortm$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_cohortm$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_cohortm$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_cohortm$age70m<-ifelse(svm.pred.info_clas_cohortm$age>=70, 1, 0)
fit <- glm(age70m ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#adt
fit <- glm(adt ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_adt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#rt
fit <- glm(rt ~ ERG_PTENdel_status, data=svm.pred.info_clas_cohortm, family = "binomial")
m_rt<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]

m_all_ref_ERG_PTENwt<-data.frame(rbind(m_gs7less, m_gs7more, GS7_4_3_more, m_lni, m_sm, m_ece, m_svi, 
                                       m_psa10l, m_psa20m, m_age))#, m_adt, m_rt))
row.names(m_all_ref_ERG_PTENwt)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "GS>7(4+3) (GS<7(3+4) as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                   "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70")#,"ADT","RT")
m_all_ref_ERG_PTENwt
#write.csv(m_all_ref_ERG_PTENwt, file="/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-UV-ERG-PTENdel_PTENwt-forestplot.csv", row.names=T)
#forest plot of odd ratio
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.03_CHD1-sig167-svm-cost_0.1-clas-Decipher-variable-UV-ERG-PTENdel_PTENwt-forestplot.pdf", height=8, width=6)
forestplot(m_all_ref_ERG_PTENwt[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_all_ref_ERG_PTENwt[,4])/7,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()



#Figure S19: Pathological outcome difference between progression events and early events in Decipher prospective cohort (n=6,532) via univariable analyses.

#Clinical outcome comparison between subgroups CHD1 and SPOP
#exclude NA in met_time 
svm.pred.info_clas_6532m<-svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_signature=="Pos" | svm.pred.info_clas_6532$SPOP_signature=="Pos",]
svm.pred.info_clas_6532m$CHD1_del<-ifelse(svm.pred.info_clas_6532m$CHD1del=="Pos", 1,0)

#Univariate analysis of between CHD1del and SPOPmut (reference: SPOP)
#GS<7 (GS=7 as ref) vs. CHD1del
svm.pred.info_clas_6532m_GS7less<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7less$pathgs))<7, 1, 0)
fit <- glm(GS7less ~ CHD1_del, data=svm.pred.info_clas_6532m_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 (GS=7 as ref) vs. CHD1del
svm.pred.info_clas_6532m_GS7more<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7more$pathgs))>7, 1, 0)
fit <- glm(GS7more ~ CHD1_del, data=svm.pred.info_clas_6532m_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_6532m$lni_status<-ifelse(svm.pred.info_clas_6532m$lni==1, 1, 0)
fit <- glm(lni_status ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_6532m$sm_status<-ifelse(svm.pred.info_clas_6532m$sm==1, 1, 0)
fit <- glm(sm_status ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#epe
svm.pred.info_clas_6532m$epe_status<-ifelse(svm.pred.info_clas_6532m$epe==1, 1, 0)
fit <- glm(epe_status ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_6532m$svi_status<-ifelse(svm.pred.info_clas_6532m$svi==1, 1, 0)
fit <- glm(svi_status ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_6532m$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_6532m$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_6532m$age70m<-ifelse(svm.pred.info_clas_6532m$age>=70, 1, 0)
fit <- glm(age70m ~ CHD1_del, data=svm.pred.info_clas_6532m, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]

m_all_refSPOP<-data.frame(rbind(m_gs7less, m_gs7more, m_lni, m_sm, m_ece, m_svi, 
                                m_psa10l, m_psa20m, m_age))
row.names(m_all_refSPOP)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                            "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70")
#write.csv(m_all_refSPOP, file="/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-UV-CHD1_SPOP-forestplot.csv", row.names=T)
#forest plot of odd ratio
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-UV-CHD1_SPOP-forestplot.pdf", height=8, width=6)
forestplot(m_all_refSPOP[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_all_refSPOP[,4])/10,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()



#Clinical outcome comparison between subgroups ERG_PTENdel and ERG_PTENwt
#exclude NA in met_time 
svm.pred.info_clas_6532m<-svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENdel=="Pos" | svm.pred.info_clas_6532$ERG_PTENwt=="Pos",]
svm.pred.info_clas_6532m$ERG_PTENdel_status<-ifelse(svm.pred.info_clas_6532m$ERG_PTENdel=="Pos", 1,0)
#Univariate analysis of between ERG_PTENdel and ERG_PTENwt (reference: ERG_PTENwt)
#GS<7 (GS=7 as ref) vs. ERG_PTEN
svm.pred.info_clas_6532m_GS7less<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))<7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7less$GS7less<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7less$pathgs))<7, 1, 0)
fit <- glm(GS7less ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m_GS7less, family = "binomial")
m_gs7less<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#GS>7 (GS=7 as ref) vs. ERG_PTEN
svm.pred.info_clas_6532m_GS7more<-svm.pred.info_clas_6532m[as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))>7 | as.numeric(as.character(svm.pred.info_clas_6532m$pathgs))==7, ]
svm.pred.info_clas_6532m_GS7more$GS7more<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m_GS7more$pathgs))>7, 1, 0)
fit <- glm(GS7more ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m_GS7more, family = "binomial")
m_gs7more<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#lni
svm.pred.info_clas_6532m$lni_status<-ifelse(svm.pred.info_clas_6532m$lni==1, 1, 0)
fit <- glm(lni_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_lni<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#sm
svm.pred.info_clas_6532m$sm_status<-ifelse(svm.pred.info_clas_6532m$sm==1, 1, 0)
fit <- glm(sm_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_sm<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#epe
svm.pred.info_clas_6532m$epe_status<-ifelse(svm.pred.info_clas_6532m$epe==1, 1, 0)
fit <- glm(epe_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_ece<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#svi
svm.pred.info_clas_6532m$svi_status<-ifelse(svm.pred.info_clas_6532m$svi==1, 1, 0)
fit <- glm(svi_status ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_svi<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA<10
svm.pred.info_clas_6532m$PSA10<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))<10, 1, 0)
fit <- glm(PSA10 ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_psa10l<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#PSA>20
svm.pred.info_clas_6532m$PSA20<-ifelse(as.numeric(as.character(svm.pred.info_clas_6532m$preop_psa))>20, 1, 0)
fit <- glm(PSA20 ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_psa20m<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]
#age
svm.pred.info_clas_6532m$age70m<-ifelse(svm.pred.info_clas_6532m$age>=70, 1, 0)
fit <- glm(age70m ~ ERG_PTENdel_status, data=svm.pred.info_clas_6532m, family = "binomial")
m_age<-cbind(or = exp(coef(fit)), exp(confint(fit, level=0.9)), "p-value" = coef(summary(fit))[, "Pr(>|z|)"])[-1,]

m_all_ref_ERG_PTENwt<-data.frame(rbind(m_gs7less, m_gs7more, m_lni, m_sm, m_ece, m_svi, 
                                       m_psa10l, m_psa20m, m_age))
row.names(m_all_ref_ERG_PTENwt)<-c("GS<7 (GS=7 as ref)", "GS>7 (GS=7 as ref)", "Lymph node invasion", "Surgical margins", "Extracapsular extension", 
                                   "Seminal vesicle invasion", "PSA<10","PSA>20", "Age>=70")
#write.csv(m_all_ref_ERG_PTENwt, file="/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-UV-ERG-PTENdel_PTENwt-forestplot.csv", row.names=T)
#forest plot of odd ratio
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-variable-UV-ERG-PTENdel_PTENwt-forestplot.pdf", height=8, width=6)
forestplot(m_all_ref_ERG_PTENwt[,1:3], 
           is.summary=F, xlog=T, zero=1, boxsize=-log10(m_all_ref_ERG_PTENwt[,4])/10,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), new_page = FALSE)
dev.off()

