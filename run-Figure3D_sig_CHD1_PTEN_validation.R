library(MASS)
library(RColorBrewer)
library(gdata)
library(e1071)
library(beeswarm)
library(forestplot)
library(gridExtra)
library(ggplot2)

#Subclass color
sur_col_ERG<-brewer.pal(9,"Set1")[1]
sur_col_CHD1_homo<-brewer.pal(9,"Set1")[2]
sur_col_CHD1_hete<-brewer.pal(12,"Set3")[5]
sur_col_other<-brewer.pal(9,"Set1")[3]

#SPOP sig 212 genes
sig212<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//TCGA-RNA-Seq/TCGA-freeze-333-ETSNeg-SPOPMut+Wild-FDR0.0001.txt-212genes.txt", sep="\t", header=T, check.names=F)
sig212_zscore<-t(scale(t(sig212[,2:(ncol(sig212)-3)]), center=T, scale=T)[,1:nrow(sig212)])
rownames(sig212_zscore)<-sig212$Gene
sigSPOP_zscore <- sig212_zscore

#CHD1 sig 148 genes: FDR<9e-04
sig148<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR9e-04.csv", check.names = F)
colnames(sig148)<-c("Gene", colnames(sig148[,2:ncol(sig148)]))
sig148_zscore<-t(scale(t(sig148[,2:ncol(sig148)]), center=T, scale=T)[,1:nrow(sig148)])
rownames(sig148_zscore)<-sig148$Gene
sigCHD1_zscore<-sig148_zscore

#CHD1 sig 167 genes
#sig167<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR0.001-333.csv", check.names = F)
#colnames(sig167)<-c("Gene", colnames(sig167[,2:ncol(sig167)]))
#sig167_zscore<-t(scale(t(sig167[,2:ncol(sig167)]), center=T, scale=T)[,1:nrow(sig167)])
#rownames(sig167_zscore)<-sig167$Gene

#CHD1_homdel sig 179 genes
#sig179<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1_homdel+CHD1wt-Wilcox-FDR0.005.csv", check.names = F)
#colnames(sig179)<-c("Gene", colnames(sig179[,2:ncol(sig179)]))
#sig179_zscore<-t(scale(t(sig179[,2:ncol(sig179)]), center=T, scale=T)[,1:nrow(sig179)])
#rownames(sig179_zscore)<-sig179$Gene

#FOXA1 sig 67 genes
sig67<-read.csv("/Users/deli/Desktop/Chris/TCGA-PCA/FOXA1/TCGA-freeze-333-FOXA1mut_wing2+FOXA1wt_ETSneg-wilcox-FDR_0.05.csv", check.names = F)
colnames(sig67)<-c("Gene", colnames(sig67[,2:ncol(sig67)]))
sig67_zscore<-t(scale(t(sig67[,2:ncol(sig67)]), center=T, scale=T)[,1:nrow(sig67)])
rownames(sig67_zscore)<-sig67$Gene

#PTEN sig 45 genes: FDR<6e-07
sigPTEN<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt-Wilcox-FDR6e-07.csv", check.names = F)
colnames(sigPTEN)<-c("Gene", colnames(sigPTEN[,2:ncol(sigPTEN)]))
sigPTEN_zscore<-t(scale(t(sigPTEN[,2:ncol(sigPTEN)]), center=T, scale=T)[,1:nrow(sigPTEN)])
rownames(sigPTEN_zscore)<-sigPTEN$Gene

#TCGA sample info.
dinfo_TCGA<-read.xls("/Users/deli/Dropbox/Papers/TCGA-Prostate/mmc2.xls", check.names=F)
dTCGAid<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31.tsv", sep="\t", header=T)
dTCGAid_ID<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31-ID.tsv", sep="\t", header=T, check.names=F)
dTCGAid<-cbind(dTCGAid, dTCGAid_ID)
dTCGAid<-dTCGAid[!grepl("TCGA-HC-7740-01B", dTCGAid$Sample.ID),]
dinfo_TCGA<-merge(dinfo_TCGA, dTCGAid, by.x="SAMPLE_ID", by.y="Sample ID")
dinfo_TCGA_SPOP<-dinfo_TCGA
dinfo_TCGA_CHD1<-dinfo_TCGA
dinfo_TCGA_CHD1$CHD1_del<-ifelse(dinfo_TCGA$CHD1_CNA!="diploid", 1, 0)
dinfo_TCGA_PTEN<-dinfo_TCGA
dinfo_TCGA_PTEN$PTEN_del<-ifelse(dinfo_TCGA$PTEN_CNA!="diploid", 1, 0)

#GPL10264 annotation
dgene_GPL10164<- read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/GSE21034_RAW/GPL10264_ref.txt", sep="\t", header=T, check.names = F)
dref_hg19<-read.table("/Users/deli/Desktop/Chris/annotation/human/refGene.txt-hg19.gene", sep="\t", header=F, check.names = F)
dgene_GPL10164<-merge(dgene_GPL10164[grepl("NM", dgene_GPL10164$GB_ACC),], dref_hg19, by.x="GB_ACC", by.y="V1")

#GSEA21034 data
#Gene expression results from cluster
dat.expgse21034<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/GSE21034_RAW/GSE21034-GPL10264_series_matrix.txt", sep="\t", header=T, check.names=F)
dat.info_GSE21034<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/GSE21034_RAW/GSE21034-GPL10264_info.txt", sep="\t", header=F, check.names = F)
dat.info_cbio_taylor_CHD1<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/GSE21034_RAW/cBio_CHD1_cna.txt", sep="\t", header=T, check.names = F)
dat.info_cbio_taylor_PTEN<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/GSE21034_RAW/cBio_PTEN_cna.txt", sep="\t", header=T, check.names = F)
#dat.exp<-dat.expgse46691
dat.infor_gse21034<-read.xls("/Users/deli/Desktop/Chris/Papers/Integrative Genomic Profiling of Human Prostate Cancer/Table-S1.xls")
dat.infor_gse21034<-dat.infor_gse21034[grepl("PCA", dat.infor_gse21034$Sample.ID), ]
dat.infor_gse21034<-dat.infor_gse21034[dat.infor_gse21034$Type=="PRIMARY", 1:37]
dat.info_cbio_taylor<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/HumanExonArray/prad_mskcc_clinical_data-Taylor-20170427.tsv", check.names=F, header=T, sep="\t")
dat.info_gse21034<-merge(merge(merge(merge(dat.infor_gse21034, dat.info_cbio_taylor, by.x="Sample.ID", by.y="Sample ID"), dat.info_GSE21034, by.x="Sample.ID", by.y="V2"), 
                         dat.info_cbio_taylor_CHD1, by.x="Sample.ID", by.y="SAMPLE_ID"), dat.info_cbio_taylor_PTEN, by.x="Sample.ID", by.y="SAMPLE_ID")

dat.expgse21034_primary<-cbind(dat.expgse21034$ID_REF, dat.expgse21034[,as.matrix(dat.info_gse21034$V1)])
colnames(dat.expgse21034_primary)<-c("ID_REF", colnames(dat.expgse21034_primary[,2:ncol(dat.expgse21034_primary)]))
dat.expgse21034_primary<-merge(dat.expgse21034_primary, dgene_GPL10164[,c(2,3)], by.x="ID_REF", by.y="ID")
dat.info_gse21034$CHD1_exp<-t(dat.expgse21034_primary[match("CHD1", dat.expgse21034_primary$V2), 2:(ncol(dat.expgse21034_primary)-1)])


#Predict SPOP, CHD1, PTEN in Taylor GSE21034 cohort
dat.expgse21034_primary_zscore<-t(scale(t(dat.expgse21034_primary[,2:(ncol(dat.expgse21034_primary)-1)]), center=T, scale=T)[,1:nrow(dat.expgse21034_primary)])
rownames(dat.expgse21034_primary_zscore)<-dat.expgse21034_primary$V2
dat.expgse21034_primary_zscore<-dat.expgse21034_primary_zscore[complete.cases(dat.expgse21034_primary_zscore),]


#CHD1 model
sig_zscore <- sigCHD1_zscore
#SVM model 
datm<-merge(sig_zscore, dat.expgse21034_primary_zscore, by="row.names")
#traning data from TCGA 
traindat<-t(datm[,2:(ncol(sig_zscore)+1)])
colnames(traindat)<-datm[,1]
traindat<-merge(traindat, dinfo_TCGA_CHD1[,c("SAMPLE_ID", "CHD1_del")], by.x="row.names", by.y="SAMPLE_ID")
traindatm<-traindat[,2:ncol(traindat)]
rownames(traindatm)<-rownames(traindat)
#testing data from GEO data
testdat<-t(datm[,(ncol(sig_zscore)+2):ncol(datm)])
colnames(testdat)<-colnames(traindatm[,1:(ncol(traindatm)-1)])
testdatm<-testdat
## svm predict linear type="C-classification"
cost_value_CHD1 <- 0.3
svm.model_clas <- svm(CHD1_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_CHD1)
svm.pred_clas<- predict(svm.model_clas, testdatm)
svm.pred.info_clas<-merge(dat.info_gse21034, testdatm[,1], by.x="V1", by.y="row.names")
svm.pred.info_clas<-merge(svm.pred.info_clas, data.frame(svm.pred_clas), by.x="V1", by.y="row.names")
svm.pred.info_clas_CHD1 <- svm.pred.info_clas[, c(1,2,ncol(svm.pred.info_clas))]
colnames(svm.pred.info_clas_CHD1) <- c(colnames(svm.pred.info_clas_CHD1[,1:2]), "svm_pred_CHD1")


#PTEN model
sig_zscore <- sigPTEN_zscore
#SVM model 
datm<-merge(sig_zscore, dat.expgse21034_primary_zscore, by="row.names")
#traning data from TCGA 
traindat<-t(datm[,2:(ncol(sig_zscore)+1)])
colnames(traindat)<-datm[,1]
traindat<-merge(traindat, dinfo_TCGA_PTEN[,c("Sample.ID", "PTEN_del")], by.x="row.names", by.y="Sample.ID")
traindatm<-traindat[,2:ncol(traindat)]
rownames(traindatm)<-rownames(traindat)
#testing data from GEO data
testdat<-t(datm[,(ncol(sig_zscore)+2):ncol(datm)])
colnames(testdat)<-colnames(traindatm[,1:(ncol(traindatm)-1)])
testdatm<-testdat
## svm predict linear type="C-classification"
#cost_value_PTEN<-0.04
cost_value_PTEN<-0.05
svm.model_clas <- svm(PTEN_del ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_PTEN)
svm.pred_clas<- predict(svm.model_clas, testdatm)
svm.pred.info_clas<-merge(dat.info_gse21034, testdatm[,1], by.x="V1", by.y="row.names")
svm.pred.info_clas<-merge(svm.pred.info_clas, data.frame(svm.pred_clas), by.x="V1", by.y="row.names")
svm.pred.info_clas<-svm.pred.info_clas[complete.cases(svm.pred.info_clas$PTEN),]
svm.pred.info_clas_PTEN <- svm.pred.info_clas[, c(1,2,ncol(svm.pred.info_clas))]
colnames(svm.pred.info_clas_PTEN) <- c(colnames(svm.pred.info_clas_PTEN[,1:2]), "svm_pred_PTEN")

#Combine prediction results
dinfo_GSE21034_pred <- merge(dat.info_gse21034, svm.pred.info_clas_CHD1[,c(2:3)], by="Sample.ID")
dinfo_GSE21034_pred <- merge(dinfo_GSE21034_pred, svm.pred.info_clas_PTEN[,c(2:3)], by="Sample.ID")

#CHD1 prediction results
dinfo_GSE21034_pred$CHD1_del_ETSneg<-ifelse(dinfo_GSE21034_pred$CHD1<0 & dinfo_GSE21034_pred$`ERG Fusion ACGH`=="negative", 1, 0)
dinfo_GSE21034_pred$CHD1homo_ETSneg<-ifelse(dinfo_GSE21034_pred$CHD1==-1 & dinfo_GSE21034_pred$`ERG Fusion ACGH`=="negative", 1, 0)
dinfo_GSE21034_pred$CHD1hete_ETSneg<-ifelse(dinfo_GSE21034_pred$CHD1==-2 & dinfo_GSE21034_pred$`ERG Fusion ACGH`=="negative", 1, 0)
TP<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_CHD1==1 & dinfo_GSE21034_pred$CHD1_del_ETSneg==1,])
FN<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$CHD1_del_ETSneg==1,]) - TP
FP<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_CHD1==1,]) - TP
TN<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$CHD1_del_ETSneg==0,]) - FP
TPR<-TP/(TP+FN)
TNR<-TN/(FP+TN)
TPR 
TNR

#PTEN prediction results
dinfo_GSE21034_pred$PTEN_del<-ifelse(dinfo_GSE21034_pred$PTEN<0 , 1, 0)
#dinfo_GSE21034_pred$PTEN_del<-ifelse(dinfo_GSE21034_pred$PTEN<0 & dinfo_GSE21034_pred$`ERG Fusion ACGH`=="positive", 1, 0)
TP<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==1 & dinfo_GSE21034_pred$PTEN_del==1,])
FN<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$PTEN_del==1,]) - TP
FP<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==1,]) - TP
TN<-nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$PTEN_del==0,]) - FP
TPR<-TP/(TP+FN)
TNR<-TN/(FP+TN)
TPR 
TNR

#Case count about PTENdel prediction 
count<-data.frame(c(nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==1&dinfo_GSE21034_pred$PTEN_del==1,]), 
                    nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==1&dinfo_GSE21034_pred$PTEN_del==0,]),
                    nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==0&dinfo_GSE21034_pred$PTEN_del==1,]),
                    nrow(dinfo_GSE21034_pred[dinfo_GSE21034_pred$svm_pred_PTEN==0&dinfo_GSE21034_pred$PTEN_del==0,])))
colnames(count)<-"Counts"
count$PTEN_predict<-c(rep("pred=1",2), rep("pred=0",2))
count$PTEN_del<-c(rep(c("PTEN_del", "PTEN_wt"), 2))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt_ETSpos-Wilcox-sig45-svm-cost-GSE21034_subtype_PTENdel-barplot.pdf", height=8, width=8)
ggplot(data=count, aes(x=PTEN_predict, y=Counts, fill=PTEN_del)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("PTEN_del"=sur_col_CHD1,  "PTEN_wt"=sur_col_other))
#dev.off()

