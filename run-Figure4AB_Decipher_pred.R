library(survival)
library(xlsx)
library(e1071)
library(RColorBrewer)
library(gdata)
library(forestplot)
library(ggplot2)
library(ggbeeswarm)
library("survminer")

#Subclass color
sur_col_ERG<-brewer.pal(9,"Set1")[1]
sur_col_ERG_PTENdel<-brewer.pal(9,"Set1")[1]
sur_col_ERG_PTENwt<-brewer.pal(12,"Set3")[8]
sur_col_ERGneg_PTENdel<-brewer.pal(12,"Set3")[3]
sur_col_ETS<-brewer.pal(9,"Set1")[6]
sur_col_ETS_PTENdel<-brewer.pal(9,"Set1")[6]
sur_col_ETS_PTENwt<-brewer.pal(12,"Set3")[2]
sur_col_CHD1<-brewer.pal(9,"Set1")[2]
sur_col_other<-brewer.pal(9,"Set1")[3]
sur_col_SPOPo<-brewer.pal(9,"Set1")[5]
sur_col_CHD1o<-brewer.pal(9,"Set1")[2]
sur_col_CHD1_SPOP<-brewer.pal(9,"Set1")[4]

#GENCODE gene type
dgene<-read.table("/Users/deli/Desktop/Chris/annotation/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gene", sep=" ", header=F)
dgene_type<-read.table("/Users/deli/Desktop/Chris/annotation/human/gencode.v19.annotation.gtf.genes_position_type", sep=" ", header=F, check.names = F)
dgene_type<-merge(dgene, dgene_type, by.x="V1", by.y="V5", sort=F)
dgene_type_coding <- dgene_type[dgene_type$V7=="protein_coding", ]
#colnames(dgene_type_coding) <- "gene"


#Decipher data input
dat_gdx_expm<-readRDS("/Users/deli/Dropbox/Deli_LabMeeting/SPOP_signature/GenomeDX/gdx_20180530_gdxExon_expression-gene.rds")
dat_gdx_expm[1:3, 1:3]
#Decipher data info.
datinfo<-read.xls("/Users/deli/Desktop/Chris/TCGA-PCA//GenomeDX/gndx-2016-04-26_clinical_expression_data.xlsx", 1) 
datinfo2<-read.xls("/Users/deli/Desktop/Chris/TCGA-PCA//GenomeDX/gndx-2016-04-26_clinical_data_key.xlsx", 1)
#Decipher exp subtype
dat_gdx_expm_sub<-t(merge(dat_gdx_expm, data.frame(c("ERG", "ETV1", "ETV4", "ETV5", "FLI1", "SPINK1")), by.x="gene", by.y="c..ERG....ETV1....ETV4....ETV5....FLI1....SPINK1.."))
dat_gdx_expm_sub<-data.frame(dat_gdx_expm_sub[2:nrow(dat_gdx_expm_sub),])
colnames(dat_gdx_expm_sub)<-c("ERG","ETV1","ETV4","ETV5","FLI1","SPINK1")
#Molecular subtype ERG, ETS, SPINK1 and TripleNeg
datinfom<-merge(datinfo, dat_gdx_expm_sub, by.x="celfile_name", by.y="row.names", sort=F)
#Exclude 116 duplicates cases from JHMI cohort
datinfom<-unique(datinfom[,c(1,3:(ncol(datinfom)-9), (ncol(datinfom)-7):ncol(datinfom))])

datinfom$ERG_overexp<-ifelse(as.numeric(as.character(datinfom$ERG))>0.6, "Pos", "Neg")
datinfom$ETS_overexp<-ifelse(as.numeric(as.character(datinfom$ETV1))>0.41| as.numeric(as.character(datinfom$ETV4))>0.32|as.numeric(as.character(datinfom$ETV5))>0.48|as.numeric(as.character(datinfom$FLI1))>0.53, "Pos", "Neg")
datinfom$SPINK1_overexp<-ifelse(as.numeric(as.character(datinfom$SPINK1))>1.03, "Pos", "Neg")
datinfom$TripleNeg<-ifelse(as.numeric(as.character(datinfom$ERG))<=0.6 & as.numeric(as.character(datinfom$ETV1))<=0.41& as.numeric(as.character(datinfom$ETV4))<=0.32&as.numeric(as.character(datinfom$ETV5))<=0.48&as.numeric(as.character(datinfom$FLI1))<=0.53
                           &as.numeric(as.character(datinfom$SPINK1))<=1.03,"Pos", "Neg")


#TCGA PCA 333 sample information
dinfo<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//TCGA-RNA-Seq/PRAD_TCGA_annotation_20150303.txt",sep="\t",header=T, check.names=F)
dinfo_TCGA<-read.xls("/Users/deli/Dropbox/Papers/TCGA-Prostate/mmc2.xls", check.names=F)
#Tumor ID
dTCGAid<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31.tsv", sep="\t", header=T)
dTCGAid_ID<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31-ID.tsv", sep="\t", header=T, check.names=F)
dTCGAid<-cbind(dTCGAid, dTCGAid_ID)
dTCGAid<-dTCGAid[!grepl("TCGA-HC-7740-01B", dTCGAid$Sample.ID),]
dinfo_TCGA<-merge(dinfo_TCGA, dTCGAid, by.x="SAMPLE_ID", by.y="Sample ID")

#TCGA annotation SPOP, CHD1, FOXA1_wing2
dinfo_TCGA_SPOP<-dinfo_TCGA[, c("SAMPLE_ID", "Sample.ID", "SPOP_mut")]
dinfo_TCGA_CHD1<-dinfo_TCGA[, c("SAMPLE_ID", "Sample.ID", "CHD1_CNA")]
#dinfo_TCGA_CHD1$CHD1_del<-ifelse(dinfo_TCGA_CHD1$CHD1_CNA!="diploid" & dinfo_TCGA$ERG_status=="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none", 1, 0)
dinfo_TCGA_CHD1$CHD1_del<-ifelse(dinfo_TCGA_CHD1$CHD1_CNA!="diploid", 1, 0)
dinfo_TCGA_PTEN<-dinfo_TCGA[, c("SAMPLE_ID", "Sample.ID", "PTEN_CNA")]
dinfo_TCGA_PTEN$PTEN_del<-ifelse(dinfo_TCGA_PTEN$PTEN_CNA!="diploid", 1, 0)

#SPOP sig 212 genes
sig212<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//TCGA-RNA-Seq/TCGA-freeze-333-ETSNeg-SPOPMut+Wild-FDR0.0001.txt-212genes.txt", sep="\t", header=T, check.names=F)
sig212_zscore<-t(scale(t(sig212[,2:(ncol(sig212)-3)]), center=T, scale=T)[,1:nrow(sig212)])
rownames(sig212_zscore)<-sig212$Gene

#CHD1 sig 148 genes: FDR<9e-04
sig148<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR9e-04.csv", check.names = F)
colnames(sig148)<-c("Gene", colnames(sig148[,2:ncol(sig148)]))
sig148_zscore<-t(scale(t(sig148[,2:ncol(sig148)]), center=T, scale=T)[,1:nrow(sig148)])
rownames(sig148_zscore)<-sig148$Gene
sigCHD1_zscore<-sig148_zscore

#PTEN sig 45 genes: FDR<6e-07
sigPTEN<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt-Wilcox-FDR6e-07.csv", check.names = F)
colnames(sigPTEN)<-c("Gene", colnames(sigPTEN[,2:ncol(sigPTEN)]))
sigPTEN_zscore<-t(scale(t(sigPTEN[,2:ncol(sigPTEN)]), center=T, scale=T)[,1:nrow(sigPTEN)])
rownames(sigPTEN_zscore)<-sigPTEN$Gene


#tesing data from GenomeDx after Z-score scale from each cohort 
dinfo_gdx<-datinfom
#CHD1 
dat_gdx_CHD1<-merge(dat_gdx_expm, data.frame(rownames(sigCHD1_zscore)), by.x="gene", by.y="rownames.sigCHD1_zscore.")
dat_gdx_CHD1<-data.frame(dat_gdx_CHD1$gene, data.frame(apply(dat_gdx_CHD1[,2:ncol(dat_gdx_CHD1)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_CHD1)<-colnames(dat_gdx_expm)
#SPOP
dat_gdx_SPOP<-merge(dat_gdx_expm, data.frame(rownames(sig212_zscore)), by.x="gene", by.y="rownames.sig212_zscore.")
dat_gdx_SPOP<-data.frame(dat_gdx_SPOP$gene, data.frame(apply(dat_gdx_SPOP[,2:ncol(dat_gdx_SPOP)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_SPOP)<-colnames(dat_gdx_expm)
#PTEN sig
dat_gdx_PTEN<-merge(dat_gdx_expm, data.frame(rownames(sigPTEN_zscore)), by.x="gene", by.y="rownames.sigPTEN_zscore.")
dat_gdx_PTEN<-data.frame(dat_gdx_PTEN$gene, data.frame(apply(dat_gdx_PTEN[,2:ncol(dat_gdx_PTEN)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_PTEN)<-colnames(dat_gdx_expm)


cost_value_CHD1<-0.3
cost_value_SPOP<-0.03
cost_value_FOXA1<-1
cost_value_ERG<-0.01
cost_value_PTEN<-0.04

#Identify the subtypes from each cohort
svm.pred.info_clas_cohort<-NULL
for (i in 1: length(summary(datinfom$study_name)))
{
  datinfom_i<-datinfom[datinfom$'study_name'==names(summary(datinfom$study_name)[i]),]
  
  #SVM model to predict PTEN
  dat_gdx_PTEN_i<-cbind(dat_gdx_PTEN$gene, dat_gdx_PTEN[, as.matrix(datinfom_i[,1])])
  dat_gdx_PTEN_i_zscore<-t(scale(t(dat_gdx_PTEN_i[,2:ncol(dat_gdx_PTEN_i)]), center=T, scale=T)[,1:nrow(dat_gdx_PTEN_i)])
  rownames(dat_gdx_PTEN_i_zscore)<-dat_gdx_PTEN_i[,1]
  dat_gdx_PTEN_i_zscore<-dat_gdx_PTEN_i_zscore[complete.cases(dat_gdx_PTEN_i_zscore),]
  datm<-merge(sigPTEN_zscore, dat_gdx_PTEN_i_zscore, by="row.names")
  #traning data from TCGA 
  traindat<-t(datm[,2:(ncol(sigPTEN_zscore)+1)])
  colnames(traindat)<-datm[,1]
  traindat<-merge(traindat, dinfo_TCGA_PTEN[,c(2,ncol(dinfo_TCGA_PTEN))] , by.x="row.names", by.y="Sample.ID")
  traindatm<-traindat[,2:ncol(traindat)]
  rownames(traindatm)<-rownames(traindat)
  #testing data from GenomeDx 
  testdat<-t(datm[,(ncol(sigPTEN_zscore)+2):ncol(datm)])
  colnames(testdat)<-datm[,1]
  testdatm<-testdat
  # svm predict linear cost=1 type="C-classification"
  svm.model_cost_PTEN <- svm(PTEN_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_PTEN)
  svm.pred_clas_PTEN <- predict(svm.model_cost_PTEN, testdatm)
  
  #SVM model to predict CHD1 by CHD1del signature
  dat_gdx_CHD1_i<-cbind(dat_gdx_CHD1$gene, dat_gdx_CHD1[, as.matrix(datinfom_i[,1])])
  dat_gdx_CHD1_i_zscore<-t(scale(t(dat_gdx_CHD1_i[,2:ncol(dat_gdx_CHD1_i)]), center=T, scale=T)[,1:nrow(dat_gdx_CHD1_i)])
  rownames(dat_gdx_CHD1_i_zscore)<-dat_gdx_CHD1_i[,1]
  dat_gdx_CHD1_i_zscore<-dat_gdx_CHD1_i_zscore[complete.cases(dat_gdx_CHD1_i_zscore),]
  datm<-merge(sigCHD1_zscore, dat_gdx_CHD1_i_zscore, by="row.names")
  #traning data from TCGA 
  traindat<-t(datm[,2:(ncol(sigCHD1_zscore)+1)])
  colnames(traindat)<-datm[,1]
  traindat<-merge(traindat, dinfo_TCGA_CHD1[,c(1,ncol(dinfo_TCGA_CHD1))] , by.x="row.names", by.y="SAMPLE_ID")
  traindatm<-traindat[,2:ncol(traindat)]
  rownames(traindatm)<-rownames(traindat)
  #testing data from GenomeDx 
  testdat<-t(datm[,(ncol(sigCHD1_zscore)+2):ncol(datm)])
  colnames(testdat)<-datm[,1]
  testdatm<-testdat
  # svm predict linear cost=1 type="C-classification"
  svm.model_cost_CHD1 <- svm(CHD1_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_CHD1)
  svm.pred_clas_CHD1 <- predict(svm.model_cost_CHD1, testdatm)
  
  #SVM model to predict SPOP
  dat_gdx_SPOP_i<-cbind(dat_gdx_SPOP$gene, dat_gdx_SPOP[, as.matrix(datinfom[datinfom$'study_name'==names(summary(datinfom$study_name)[i]),1])])
  dat_gdx_SPOP_i_zscore<-t(scale(t(dat_gdx_SPOP_i[,2:ncol(dat_gdx_SPOP_i)]), center=T, scale=T)[,1:nrow(dat_gdx_SPOP_i)])
  rownames(dat_gdx_SPOP_i_zscore)<-dat_gdx_SPOP_i[,1]
  dat_gdx_SPOP_i_zscore<-dat_gdx_SPOP_i_zscore[complete.cases(dat_gdx_SPOP_i_zscore),]
  datm<-merge(sig212_zscore, dat_gdx_SPOP_i_zscore, by="row.names")
  #traning data from TCGA 
  traindat<-t(datm[,2:(ncol(sig212_zscore)+1)])
  colnames(traindat)<-datm[,1]
  traindat<-merge(traindat, dinfo_TCGA_SPOP[,c(1,ncol(dinfo_TCGA_SPOP))] , by.x="row.names", by.y="SAMPLE_ID")
  traindatm<-traindat[,2:ncol(traindat)]
  rownames(traindatm)<-rownames(traindat)
  #testing data from GenomeDx 
  testdat<-t(datm[,(ncol(sig212_zscore)+2):ncol(datm)])
  colnames(testdat)<-datm[,1]
  testdatm<-testdat
  # svm predict linear cost=1 type="C-classification"
  svm.model_cost_SPOP <- svm(SPOP_mut ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_SPOP)
  svm.pred_clas_SPOP <- predict(svm.model_cost_SPOP, testdatm)
  
  #Combine prediction
  svm.pred.info_clas_i<-merge(datinfom_i, testdatm[,1], by.x="celfile_name", by.y="row.names", sort=F)
  #svm.pred.info_clas_i$pred_ERG<-as.matrix(svm.pred_clas_ERG)
  #svm.pred.info_clas_i$pred_ETS<-as.matrix(svm.pred_clas_ETS)
  svm.pred.info_clas_i$pred_CHD1<-as.matrix(svm.pred_clas_CHD1)
  svm.pred.info_clas_i$pred_SPOP<-as.matrix(svm.pred_clas_SPOP)
  svm.pred.info_clas_i$pred_PTEN<-as.matrix(svm.pred_clas_PTEN)
  
  svm.pred.info_clas_cohort<-rbind(svm.pred.info_clas_cohort, svm.pred.info_clas_i)
}


svm.pred.info_clas_cohort$ERG_overexp_status<-ifelse(svm.pred.info_clas_cohort$ERG_overexp=="Pos", 1, 0)
svm.pred.info_clas_cohort$ETS_overexp_status<-ifelse(svm.pred.info_clas_cohort$ETS_overexp=="Pos", 1, 0)

#Call the ERG/ETS first, then SPOP/CHD1/PTEN calling
svm.pred.info_clas_cohort$ERG_status<-ifelse(svm.pred.info_clas_cohort$ERG_overexp=="Pos" & svm.pred.info_clas_cohort$ETS_overexp=="Neg", "Pos", "Neg")
svm.pred.info_clas_cohort$ETS_status<-ifelse(svm.pred.info_clas_cohort$ETS_overexp=="Pos" & svm.pred.info_clas_cohort$ERG_overexp=="Neg", "Pos", "Neg")

svm.pred.info_clas_cohort$SPOPmut<-ifelse(svm.pred.info_clas_cohort$pred_SPOP==1 & (svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$ETS_status=="Neg"), "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1del<-ifelse(svm.pred.info_clas_cohort$pred_CHD1==1 & (svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$ETS_status=="Neg"), "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1_SPOP_share<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Pos" & svm.pred.info_clas_cohort$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_cohort$CHD1_signature<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Neg" & svm.pred.info_clas_cohort$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_cohort$SPOP_signature<-ifelse(svm.pred.info_clas_cohort$SPOPmut=="Pos" & svm.pred.info_clas_cohort$CHD1del=="Neg", "Pos", "Neg")

svm.pred.info_clas_cohort$ERG_PTENdel<-ifelse(svm.pred.info_clas_cohort$ERG_status=="Pos" & svm.pred.info_clas_cohort$pred_PTEN==1, "Pos", "Neg")
svm.pred.info_clas_cohort$ERG_PTENwt<-ifelse(svm.pred.info_clas_cohort$ERG_status=="Pos" & svm.pred.info_clas_cohort$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_cohort$ETS_PTENdel<-ifelse((svm.pred.info_clas_cohort$ERG_status=="Pos"|svm.pred.info_clas_cohort$ETS_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==1, "Pos", "Neg")
svm.pred.info_clas_cohort$ETS_PTENwt<-ifelse((svm.pred.info_clas_cohort$ERG_status=="Pos"|svm.pred.info_clas_cohort$ETS_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_cohort$other<-ifelse(svm.pred.info_clas_cohort$ETS_status=="Neg" & svm.pred.info_clas_cohort$ERG_status=="Neg" 
                                        &svm.pred.info_clas_cohort$SPOPmut=="Neg" & svm.pred.info_clas_cohort$CHD1del=="Neg", "Pos", "Neg")

nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",])

#Barplot of molecular subgroup in each cohort
svm.pred.info_clas_CCF<-svm.pred.info_clas_cohort[grepl("CCF-", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_EMC<-svm.pred.info_clas_cohort[grepl("EMC-", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_JHMI<-svm.pred.info_clas_cohort[grepl("JHMI-", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_Mayo1<-svm.pred.info_clas_cohort[grepl("Mayo-Jenkins-MetsDisc", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_Mayo2<-svm.pred.info_clas_cohort[grepl("Mayo-Jenkins-MetsVal", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_MSKCC<-svm.pred.info_clas_cohort[grepl("MSKCC-", svm.pred.info_clas_cohort$study_name),]
svm.pred.info_clas_TJU<-svm.pred.info_clas_cohort[grepl("TJU-", svm.pred.info_clas_cohort$study_name),]

stat_CCF<-c(nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$SPOP_signature=="Pos",]),
            nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$CHD1_SPOP_share=="Pos",]),
            #nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$ERG_status=="Pos",]), 
            nrow(svm.pred.info_clas_CCF[(svm.pred.info_clas_CCF$ERG_status=="Pos") & svm.pred.info_clas_CCF$pred_PTEN==1,]),
            nrow(svm.pred.info_clas_CCF[(svm.pred.info_clas_CCF$ERG_status=="Pos") & svm.pred.info_clas_CCF$pred_PTEN==0,]),
            nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$ETS_status=="Pos",]),
            nrow(svm.pred.info_clas_CCF[svm.pred.info_clas_CCF$other=="Pos",]))/nrow(svm.pred.info_clas_CCF)
stat_EMC<-c(nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$SPOP_signature=="Pos",]),
            nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$CHD1_SPOP_share=="Pos",]),
            #nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$ERG_status=="Pos",]), 
            nrow(svm.pred.info_clas_EMC[(svm.pred.info_clas_EMC$ERG_status=="Pos") & svm.pred.info_clas_EMC$pred_PTEN==1,]),
            nrow(svm.pred.info_clas_EMC[(svm.pred.info_clas_EMC$ERG_status=="Pos") & svm.pred.info_clas_EMC$pred_PTEN==0,]),
            nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$ETS_status=="Pos",]),
            nrow(svm.pred.info_clas_EMC[svm.pred.info_clas_EMC$other=="Pos",]))/nrow(svm.pred.info_clas_EMC)
stat_JHMI<-c(nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$SPOP_signature=="Pos",]),
             nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$CHD1_SPOP_share=="Pos",]),
             #nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$ERG_status=="Pos",]), 
             nrow(svm.pred.info_clas_JHMI[(svm.pred.info_clas_JHMI$ERG_status=="Pos") & svm.pred.info_clas_JHMI$pred_PTEN==1,]),
             nrow(svm.pred.info_clas_JHMI[(svm.pred.info_clas_JHMI$ERG_status=="Pos") & svm.pred.info_clas_JHMI$pred_PTEN==0,]),
             nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$ETS_status=="Pos",]),
             nrow(svm.pred.info_clas_JHMI[svm.pred.info_clas_JHMI$other=="Pos",]))/nrow(svm.pred.info_clas_JHMI)
stat_Mayo1<-c(nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$SPOP_signature=="Pos",]),
              nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$CHD1_SPOP_share=="Pos",]),
              #nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$ERG_status=="Pos",]), 
              nrow(svm.pred.info_clas_Mayo1[(svm.pred.info_clas_Mayo1$ERG_status=="Pos") & svm.pred.info_clas_Mayo1$pred_PTEN==1,]),
              nrow(svm.pred.info_clas_Mayo1[(svm.pred.info_clas_Mayo1$ERG_status=="Pos") & svm.pred.info_clas_Mayo1$pred_PTEN==0,]),
              nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$ETS_status=="Pos",]),
              nrow(svm.pred.info_clas_Mayo1[svm.pred.info_clas_Mayo1$other=="Pos",]))/nrow(svm.pred.info_clas_Mayo1)
stat_Mayo2<-c(nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$SPOP_signature=="Pos",]),
              nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$CHD1_SPOP_share=="Pos",]),
              #nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$ERG_status=="Pos",]), 
              nrow(svm.pred.info_clas_Mayo2[(svm.pred.info_clas_Mayo2$ERG_status=="Pos") & svm.pred.info_clas_Mayo2$pred_PTEN==1,]),
              nrow(svm.pred.info_clas_Mayo2[(svm.pred.info_clas_Mayo2$ERG_status=="Pos") & svm.pred.info_clas_Mayo2$pred_PTEN==0,]),
              nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$ETS_status=="Pos",]),
              nrow(svm.pred.info_clas_Mayo2[svm.pred.info_clas_Mayo2$other=="Pos",]))/nrow(svm.pred.info_clas_Mayo2)
stat_MSKCC<-c(nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$SPOP_signature=="Pos",]),
              nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$CHD1_SPOP_share=="Pos",]),
              #nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$ERG_status=="Pos",]), 
              nrow(svm.pred.info_clas_MSKCC[(svm.pred.info_clas_MSKCC$ERG_status=="Pos") & svm.pred.info_clas_MSKCC$pred_PTEN==1,]),
              nrow(svm.pred.info_clas_MSKCC[(svm.pred.info_clas_MSKCC$ERG_status=="Pos") & svm.pred.info_clas_MSKCC$pred_PTEN==0,]),
              nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$ETS_status=="Pos",]),
              nrow(svm.pred.info_clas_MSKCC[svm.pred.info_clas_MSKCC$other=="Pos",]))/nrow(svm.pred.info_clas_MSKCC)
stat_TJU<-c(nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$SPOP_signature=="Pos",]),
            nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$CHD1_SPOP_share=="Pos",]),
            #nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$ERG_status=="Pos",]), 
            nrow(svm.pred.info_clas_TJU[(svm.pred.info_clas_TJU$ERG_status=="Pos") & svm.pred.info_clas_TJU$pred_PTEN==1,]),
            nrow(svm.pred.info_clas_TJU[(svm.pred.info_clas_TJU$ERG_status=="Pos") & svm.pred.info_clas_TJU$pred_PTEN==0,]),
            nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$ETS_status=="Pos",]),
            nrow(svm.pred.info_clas_TJU[svm.pred.info_clas_TJU$other=="Pos",]))/nrow(svm.pred.info_clas_TJU)
stat_combine<-c(nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",]),
            nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",]),  nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",]),
            #nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",]), 
            nrow(svm.pred.info_clas_cohort[(svm.pred.info_clas_cohort$ERG_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==1,]),
            nrow(svm.pred.info_clas_cohort[(svm.pred.info_clas_cohort$ERG_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==0,]),
            nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",]),
            nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",]))/nrow(svm.pred.info_clas_cohort)
#plot.col=c(brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[2], brewer.pal(9,"Set1")[4], 
#           brewer.pal(9,"Set1")[1],brewer.pal(9,"Set1")[8],brewer.pal(9,"Set1")[3])
plot.col=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS, sur_col_other)
#summary(svm.pred.info_clas_cohort$'study_name')
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-barplot.pdf", 
#    width=6, height=6, useDingbats=FALSE)
barplot(c(stat_CCF, 0, stat_EMC, 0, stat_JHMI, 0, 
          stat_Mayo1, 0, stat_Mayo2, 0, stat_MSKCC, 0, stat_TJU, 0, 
          stat_combine), col=c(plot.col,1), xaxt="n", yaxt="n",
        ylim=c(0,0.5), ylab="Subclass %")
legend("topright",c("SPOP_signature", "CHD1_signature", "CHD1del+SPOPmut","ERG+_PTENdel", "ERG+_PTENwt", "ETS+", "Other"),
       pch=16, col=plot.col, bty="n", cex=0.75)
axis_start<-2
axis_start_2<-4
axis_ins<-10
#seq(axis_start, (axis_start+axis_ins)*7, by=axis_ins)
axis(1, at=seq(axis_start, (axis_start+axis_ins)*6, by=axis_ins), label=c("CCF","EMC","JHMI","Mayo-Disc","Mayo-Val","MSKCC","TJU", "Total"), las=2, col="white", cex=0.75)
axis(1, at=seq(axis_start_2, ((axis_start+axis_ins)*6+axis_start), by=axis_ins), label=c(summary(svm.pred.info_clas_cohort$'study_name'), nrow(svm.pred.info_clas_cohort)), las=2, col="white", cex=0.75)
axis(2, at=c(0,0.5), label=c(0, 50))
#abline(h=0.05)
dev.off()

stat_cohort<-data.frame(rbind(stat_CCF, stat_EMC, stat_JHMI, stat_Mayo1, stat_Mayo2, stat_MSKCC, stat_TJU, stat_TJU))
colnames(stat_cohort)<-c("SPOP_signature", "CHD1_signature", "CHD1del+SPOPmut","ERG+_PTENdel", "ERG+_PTENwt", "ETS+", "Other")


#Pie of molecular subclasses
svm.pred.info_clas_cohort_slices <- c(nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",]),
                                      nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",]), 
                                      nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",]),
                                      nrow(svm.pred.info_clas_cohort[(svm.pred.info_clas_cohort$ERG_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==1,]),
                                      nrow(svm.pred.info_clas_cohort[(svm.pred.info_clas_cohort$ERG_status=="Pos") & svm.pred.info_clas_cohort$pred_PTEN==0,]),
                                      nrow(svm.pred.info_clas_cohort[(svm.pred.info_clas_cohort$ETS_status=="Pos"),]),
                                      nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",]))
svm.pred.info_clas_cohort_labels <- c("SPOPmut_signature", "CHD1del_signature", "CHD1del+SPOPmut", "ERG+_PTENdel", "ERG+_PTENwt", "ETS+", "Other")
svm.pred.info_clas_cohort_pct <- round(svm.pred.info_clas_cohort_slices/sum(svm.pred.info_clas_cohort_slices)*100)
svm.pred.info_clas_cohort_labels <- paste(svm.pred.info_clas_cohort_labels, svm.pred.info_clas_cohort_pct) # add percents to labels 
svm.pred.info_clas_cohort_labels <- paste(svm.pred.info_clas_cohort_labels,"%", sep="") # ad % to labels 
#plot.col=c(brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[2], brewer.pal(9,"Set1")[4], 
#           brewer.pal(12,"Set3")[8], brewer.pal(9,"Set1")[1], brewer.pal(12,"Set3")[10], brewer.pal(9,"Set1")[8], brewer.pal(9,"Set1")[3])
plot.col=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS_PTENdel, sur_col_other)
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher_PTEN-pie.pdf", width=6, height=6)
pie(svm.pred.info_clas_cohort_slices, labels =svm.pred.info_clas_cohort_labels, cex=0.7,
    col=plot.col,
    main="Decipher-1626-molecular-subclass-pie_chart")
dev.off()



#GenomeDx 6532 all gene expression dataset
data.prospective_all<-readRDS("/Users/deli/Desktop/Chris/TCGA-PCA/GenomeDX/Prospective_cohort_20190424/gdx_20190422_Prospective_gdxExon_expression.rds")
#gene<-rownames(data.prospective_all)
#data.prospective_all_data<-cbind(gene, data.prospective_all)
data.prospective_all[1:3, 1:3]
#GenomeDx 6532 sample clinical info.
data.prospective_info<-read.csv("/Users/deli/Desktop/Chris/TCGA-PCA/GenomeDX/gndx-6532-data.prospective.csv")
data.prospective_clnic<-data.prospective_info[,1:11]
#Add sample info
dat_gdx_expm_sub<-t(merge(data.prospective_all, data.frame(c("ERG", "ETV1", "ETV4", "ETV5", "FLI1", "SPINK1")), by.x="row.names", by.y="c..ERG....ETV1....ETV4....ETV5....FLI1....SPINK1.."))
dat_gdx_expm_sub<-data.frame(dat_gdx_expm_sub[2:nrow(dat_gdx_expm_sub),])
colnames(dat_gdx_expm_sub)<-c("ERG","ETV1","ETV4","ETV5","FLI1","SPINK1")
datinfo_6532<-merge(data.prospective_clnic, dat_gdx_expm_sub, by.x="celfile_name", by.y="row.names", sort=F)

datinfo_6532$ERG_overexp<-ifelse(as.numeric(as.character(datinfo_6532$ERG))>0.6, "Pos", "Neg")
datinfo_6532$ETS_overexp<-ifelse(as.numeric(as.character(datinfo_6532$ETV1))>0.41| as.numeric(as.character(datinfo_6532$ETV4))>0.32|as.numeric(as.character(datinfo_6532$ETV5))>0.48|as.numeric(as.character(datinfo_6532$FLI1))>0.53, "Pos", "Neg")
datinfo_6532$SPINK1_overexp<-ifelse(as.numeric(as.character(datinfo_6532$SPINK1))>1.03, "Pos", "Neg")
datinfo_6532$TripleNeg<-ifelse(as.numeric(as.character(datinfo_6532$ERG))<=0.6 & as.numeric(as.character(datinfo_6532$ETV1))<=0.41& as.numeric(as.character(datinfo_6532$ETV4))<=0.32&as.numeric(as.character(datinfo_6532$ETV5))<=0.48&as.numeric(as.character(datinfo_6532$FLI1))<=0.53
                               &as.numeric(as.character(datinfo_6532$SPINK1))<=1.03,"Pos", "Neg")

#tesing data from GenomeDx after Z-score scale from each cohort 
#dinfo_gdx<-datinfo_6532
#CHD1 
dat_gdx_CHD1<-merge(data.prospective_all, data.frame(rownames(sigCHD1_zscore)), by.x="row.names", by.y="rownames.sigCHD1_zscore.")
dat_gdx_CHD1<-data.frame(dat_gdx_CHD1$Row.names, data.frame(apply(dat_gdx_CHD1[,2:ncol(dat_gdx_CHD1)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_CHD1)<-c("gene", colnames(data.prospective_all))
#SPOP
dat_gdx_SPOP<-merge(data.prospective_all, data.frame(rownames(sig212_zscore)), by.x="row.names", by.y="rownames.sig212_zscore.")
dat_gdx_SPOP<-data.frame(dat_gdx_SPOP$Row.names, data.frame(apply(dat_gdx_SPOP[,2:ncol(dat_gdx_SPOP)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_SPOP)<-c("gene", colnames(data.prospective_all))
#PTEN sig
dat_gdx_PTEN<-merge(data.prospective_all, data.frame(rownames(sigPTEN_zscore)), by.x="row.names", by.y="rownames.sigPTEN_zscore.")
dat_gdx_PTEN<-data.frame(dat_gdx_PTEN$Row.names, data.frame(apply(dat_gdx_PTEN[,2:ncol(dat_gdx_PTEN)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_PTEN)<-c("gene", colnames(data.prospective_all))

#Identify the subtypes from each sample
cost_value_CHD1<-0.3
cost_value_SPOP<-0.03
cost_value_FOXA1<-1
cost_value_ERG<-0.01
cost_value_PTEN<-0.04

#SVM model to predict PTEN
dat_gdx_PTEN_i<-cbind(dat_gdx_PTEN$gene, dat_gdx_PTEN[, as.matrix(datinfo_6532[,1])])
dat_gdx_PTEN_i_zscore<-t(scale(t(dat_gdx_PTEN_i[,2:ncol(dat_gdx_PTEN_i)]), center=T, scale=T)[,1:nrow(dat_gdx_PTEN_i)])
rownames(dat_gdx_PTEN_i_zscore)<-dat_gdx_PTEN_i[,1]
dat_gdx_PTEN_i_zscore<-dat_gdx_PTEN_i_zscore[complete.cases(dat_gdx_PTEN_i_zscore),]
datm<-merge(sigPTEN_zscore, dat_gdx_PTEN_i_zscore, by="row.names")
#traning data from TCGA 
traindat<-t(datm[,2:(ncol(sigPTEN_zscore)+1)])
colnames(traindat)<-datm[,1]
traindat<-merge(traindat, dinfo_TCGA_PTEN[,c(2,ncol(dinfo_TCGA_PTEN))] , by.x="row.names", by.y="Sample.ID")
traindatm<-traindat[,2:ncol(traindat)]
rownames(traindatm)<-rownames(traindat)
#testing data from GenomeDx 
testdat<-t(datm[,(ncol(sigPTEN_zscore)+2):ncol(datm)])
colnames(testdat)<-datm[,1]
testdatm<-testdat
# svm predict linear cost=1 type="C-classification"
svm.model_cost_PTEN <- svm(PTEN_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_PTEN)
svm.pred_clas_PTEN <- predict(svm.model_cost_PTEN, testdatm)
# svm predict linear cost=1 type="C-classification" probability
svm.model_cost_PTEN_prob <- svm(PTEN_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_PTEN, probability=TRUE)
svm.pred_clas_PTEN_prob <- predict(svm.model_cost_PTEN_prob, testdatm, probability=TRUE)
svm.pred_clas_PTEN_prob<-attr(svm.pred_clas_PTEN_prob, "probabilities")

#SVM model to predict CHD1 
dat_gdx_CHD1_i<-cbind(dat_gdx_CHD1$gene, dat_gdx_CHD1[, as.matrix(datinfo_6532[,1])])
dat_gdx_CHD1_i_zscore<-t(scale(t(dat_gdx_CHD1_i[,2:ncol(dat_gdx_CHD1_i)]), center=T, scale=T)[,1:nrow(dat_gdx_CHD1_i)])
rownames(dat_gdx_CHD1_i_zscore)<-dat_gdx_CHD1_i[,1]
dat_gdx_CHD1_i_zscore<-dat_gdx_CHD1_i_zscore[complete.cases(dat_gdx_CHD1_i_zscore),]
datm<-merge(sigCHD1_zscore, dat_gdx_CHD1_i_zscore, by="row.names")
#traning data from TCGA 
traindat<-t(datm[,2:(ncol(sigCHD1_zscore)+1)])
colnames(traindat)<-datm[,1]
#traindat<-merge(traindat, dinfo_TCGA_CHD1[,c(2,ncol(dinfo_TCGA_CHD1))] , by.x="row.names", by.y="Sample.ID")
traindat<-merge(traindat, dinfo_TCGA_CHD1[,c(1,ncol(dinfo_TCGA_CHD1))] , by.x="row.names", by.y="SAMPLE_ID")
traindatm<-traindat[,2:ncol(traindat)]
rownames(traindatm)<-rownames(traindat)
#testing data from GenomeDx 
testdat<-t(datm[,(ncol(sigCHD1_zscore)+2):ncol(datm)])
colnames(testdat)<-datm[,1]
testdatm<-testdat
# svm predict linear cost=1 type="C-classification"
svm.model_cost_CHD1 <- svm(CHD1_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_CHD1)
svm.pred_clas_CHD1 <- predict(svm.model_cost_CHD1, testdatm)
# svm predict linear cost=1 type="C-classification" probability
svm.model_cost_CHD1_prob <- svm(CHD1_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_CHD1, probability=TRUE)
svm.pred_clas_CHD1_prob <- predict(svm.model_cost_CHD1_prob, testdatm, probability=TRUE)
svm.pred_clas_CHD1_prob<-attr(svm.pred_clas_CHD1_prob, "probabilities")

#SVM model to predict SPOP
dat_gdx_SPOP_i<-cbind(dat_gdx_SPOP$gene, dat_gdx_SPOP[, as.matrix(datinfo_6532[,1])])
dat_gdx_SPOP_i_zscore<-t(scale(t(dat_gdx_SPOP_i[,2:ncol(dat_gdx_SPOP_i)]), center=T, scale=T)[,1:nrow(dat_gdx_SPOP_i)])
rownames(dat_gdx_SPOP_i_zscore)<-dat_gdx_SPOP_i[,1]
dat_gdx_SPOP_i_zscore<-dat_gdx_SPOP_i_zscore[complete.cases(dat_gdx_SPOP_i_zscore),]
datm<-merge(sig212_zscore, dat_gdx_SPOP_i_zscore, by="row.names")
#traning data from TCGA 
traindat<-t(datm[,2:(ncol(sig212_zscore)+1)])
colnames(traindat)<-datm[,1]
traindat<-merge(traindat, dinfo_TCGA_SPOP[,c(1,ncol(dinfo_TCGA_SPOP))] , by.x="row.names", by.y="SAMPLE_ID")
traindatm<-traindat[,2:ncol(traindat)]
rownames(traindatm)<-rownames(traindat)
#testing data from GenomeDx 
testdat<-t(datm[,(ncol(sig212_zscore)+2):ncol(datm)])
colnames(testdat)<-datm[,1]
testdatm<-testdat
# svm predict linear cost=1 type="C-classification"
svm.model_cost_SPOP <- svm(SPOP_mut ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_SPOP)
svm.pred_clas_SPOP <- predict(svm.model_cost_SPOP, testdatm)
# svm predict linear cost=1 type="C-classification" probability
svm.model_cost_SPOP_prob <- svm(SPOP_mut ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_SPOP, probability=TRUE)
svm.pred_clas_SPOP_prob <- predict(svm.model_cost_SPOP_prob, testdatm, probability=TRUE)
svm.pred_clas_SPOP_prob<-attr(svm.pred_clas_SPOP_prob, "probabilities")

#Combine prediction
#svm.pred.info_clas_6532<-merge(datinfo_6532, data.frame(testdatm[,1]), by.x="celfile_name", by.y="row.names", sort=F)
svm.pred.info_clas_6532<-datinfo_6532
svm.pred.info_clas_6532$pred_CHD1<-as.matrix(svm.pred_clas_CHD1)
svm.pred.info_clas_6532$pred_SPOP<-as.matrix(svm.pred_clas_SPOP)
svm.pred.info_clas_6532$pred_PTEN<-as.matrix(svm.pred_clas_PTEN)
svm.pred.info_clas_6532$pred_CHD1_prob<-as.matrix(svm.pred_clas_CHD1_prob[,2])
svm.pred.info_clas_6532$pred_SPOP_prob<-as.matrix(svm.pred_clas_SPOP_prob[,2])
svm.pred.info_clas_6532$pred_PTEN_prob<-as.matrix(svm.pred_clas_PTEN_prob[,2])  
svm.pred.info_clas_6532$ERG_overexp_status<-ifelse(svm.pred.info_clas_6532$ERG_overexp=="Pos", 1, 0)
svm.pred.info_clas_6532$ETS_overexp_status<-ifelse(svm.pred.info_clas_6532$ETS_overexp=="Pos", 1, 0)

#Call the ERG/ETS first, then SPOP/CHD1/PTEN calling
svm.pred.info_clas_6532$ERG_status<-ifelse(svm.pred.info_clas_6532$ERG_overexp=="Pos" & svm.pred.info_clas_6532$ETS_overexp=="Neg", "Pos", "Neg")
svm.pred.info_clas_6532$ETS_status<-ifelse(svm.pred.info_clas_6532$ETS_overexp=="Pos" & svm.pred.info_clas_6532$ERG_overexp=="Neg", "Pos", "Neg")

svm.pred.info_clas_6532$SPOPmut<-ifelse(svm.pred.info_clas_6532$pred_SPOP==1 & (svm.pred.info_clas_6532$ERG_status=="Neg" & svm.pred.info_clas_6532$ETS_status=="Neg"), "Pos", "Neg")
svm.pred.info_clas_6532$CHD1del<-ifelse(svm.pred.info_clas_6532$pred_CHD1==1 & (svm.pred.info_clas_6532$ERG_status=="Neg" & svm.pred.info_clas_6532$ETS_status=="Neg"), "Pos", "Neg")
svm.pred.info_clas_6532$CHD1_SPOP_share<-ifelse(svm.pred.info_clas_6532$SPOPmut=="Pos" & svm.pred.info_clas_6532$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_6532$CHD1_signature<-ifelse(svm.pred.info_clas_6532$SPOPmut=="Neg" & svm.pred.info_clas_6532$CHD1del=="Pos", "Pos", "Neg")
svm.pred.info_clas_6532$SPOP_signature<-ifelse(svm.pred.info_clas_6532$SPOPmut=="Pos" & svm.pred.info_clas_6532$CHD1del=="Neg", "Pos", "Neg")

svm.pred.info_clas_6532$ERG_PTENdel<-ifelse(svm.pred.info_clas_6532$ERG_status=="Pos" & svm.pred.info_clas_6532$pred_PTEN==1, "Pos", "Neg")
svm.pred.info_clas_6532$ERG_PTENwt<-ifelse(svm.pred.info_clas_6532$ERG_status=="Pos" & svm.pred.info_clas_6532$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_6532$ETS_PTENdel<-ifelse((svm.pred.info_clas_6532$ERG_status=="Pos"|svm.pred.info_clas_6532$ETS_status=="Pos") & svm.pred.info_clas_6532$pred_PTEN==1, "Pos", "Neg")
svm.pred.info_clas_6532$ETS_PTENwt<-ifelse((svm.pred.info_clas_6532$ERG_status=="Pos"|svm.pred.info_clas_6532$ETS_status=="Pos") & svm.pred.info_clas_6532$pred_PTEN==0, "Pos", "Neg")
svm.pred.info_clas_6532$other<-ifelse(svm.pred.info_clas_6532$ETS_status=="Neg" & svm.pred.info_clas_6532$ERG_status=="Neg" 
                                      &svm.pred.info_clas_6532$SPOPmut=="Neg" & svm.pred.info_clas_6532$CHD1del=="Neg", "Pos", "Neg")

nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ETS_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$other=="Pos",])


#Pie of molecular subclasses with PTENdel subtypes
svm.pred.info_clas_6532_slices <- c(nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$SPOP_signature=="Pos",]),
                                    nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_signature=="Pos",]), 
                                    nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_SPOP_share=="Pos",]),
                                    nrow(svm.pred.info_clas_6532[(svm.pred.info_clas_6532$ERG_status=="Pos") & svm.pred.info_clas_6532$pred_PTEN==1,]),
                                    nrow(svm.pred.info_clas_6532[(svm.pred.info_clas_6532$ERG_status=="Pos") & svm.pred.info_clas_6532$pred_PTEN==0,]),
                                    nrow(svm.pred.info_clas_6532[(svm.pred.info_clas_6532$ETS_status=="Pos"),]),
                                    nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$other=="Pos",]))
svm.pred.info_clas_6532_labels <- c("SPOPmut_signature", "CHD1del_signature", "CHD1del+SPOPmut", "ERG+_PTENdel", "ERG+_PTENwt", "ETS+", "Other")
svm.pred.info_clas_6532_pct <- round(svm.pred.info_clas_6532_slices/sum(svm.pred.info_clas_6532_slices)*100)
svm.pred.info_clas_6532_labels <- paste(svm.pred.info_clas_6532_labels, svm.pred.info_clas_6532_pct) # add percents to labels 
svm.pred.info_clas_6532_labels <- paste(svm.pred.info_clas_6532_labels,"%", sep="") # ad % to labels 
#plot.col=c(brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[2], brewer.pal(9,"Set1")[4], 
#           brewer.pal(12,"Set3")[8], brewer.pal(9,"Set1")[1], brewer.pal(12,"Set3")[10], brewer.pal(9,"Set1")[8], brewer.pal(9,"Set1")[3])
plot.col=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS_PTENdel, sur_col_other)
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-pie.pdf", width=6, height=6)
pie(svm.pred.info_clas_6532_slices, labels =svm.pred.info_clas_6532_labels, cex=0.7,
    col=plot.col,
    main="Decipher-6532-molecular-subclass-pie_chart")
dev.off()





