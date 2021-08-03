library(MASS)
library(RUVSeq)
library(RColorBrewer)
library(rgl)
library(DESeq2)
library("pheatmap")
library(eulerr)
library(RColorBrewer)
library(ImpulseDE2)
#library(ComplexHeatmap)
library(VennDiagram)
library(gdata)


#gene<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/TCGA-RNA-Seq/TCGA-PC-gene",sep=" ",header=T)
dTCGAid<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//TCGA-RNA-Seq/Download/FILE_SAMPLE_MAP_ID.txt", header=F)
dTCGAnid1<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//Priyanka/final_normal_27_meth_id.txt", header=F)
dinfo_HTSeq<-read.table(file = '/Users/deli/Desktop/Chris/TCGA-PCA/HTSeq/gdc_sample_sheet.2018-03-09.tsv', sep = '\t', header = TRUE)
dinfo_HTSeq_n<-dinfo_HTSeq[grepl("Normal", dinfo_HTSeq$Sample.Type),]
#dat<-read.table(gzfile("/Users/deli/Desktop/Chris/TCGA-PCA//HTSeq/120710_UNC12-SN629_0215_BC0WRNACXX_TGACCA_L005.sorted.genecode.htscount.gz"), header=F)

#TCGA PCA 333 sample information
dinfo_TCGA<-read.xls("/Users/del2017/Dropbox/Papers/TCGA-Prostate/mmc2.xls", check.names=F)
#Tumor ID
dTCGAid<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31.tsv", sep="\t", header=T)
dTCGAid_ID<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31-ID.tsv", sep="\t", header=T, check.names=F)
dTCGAid<-cbind(dTCGAid, dTCGAid_ID)
dTCGAid<-dTCGAid[!grepl("TCGA-HC-7740-01B", dTCGAid$Sample.ID),]
dinfo_TCGA<-merge(dinfo_TCGA, dTCGAid, by.x="SAMPLE_ID", by.y="Sample ID")

#CHD1 del and ERG/ETS negative
dinfo_TCGA_CHD1del<-dinfo_TCGA[dinfo_TCGA$CHD1_CNA!="diploid" & dinfo_TCGA$ERG_status=="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none",]
#SPOP mut and ERG/ETS negative
dinfo_TCGA_SPOPmut<-dinfo_TCGA[dinfo_TCGA$SPOP_mut==1 & dinfo_TCGA$ERG_status=="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none",]
#SPOP mut and CHD1_del ERG/ETS negative
dinfo_TCGA_SPOPmut_CHD1del<-dinfo_TCGA[(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_SPOPmut$Sample.ID)  & dinfo_TCGA$CHD1_CNA!="diploid" ,]
#SPOP mut and CHD1_wt ERG/ETS negative
dinfo_TCGA_SPOPmut_CHD1wt<-dinfo_TCGA[(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_SPOPmut$Sample.ID)  & dinfo_TCGA$CHD1_CNA=="diploid" ,]
#CHD1 del and SPOP_wt ERG/ETS negative
dinfo_TCGA_SPOPwt_CHD1del<-dinfo_TCGA[(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_CHD1del$Sample.ID)  & dinfo_TCGA$SPOP_mut==0 ,]
#ERGpos 
dinfo_TCGA_ERGpos<-dinfo_TCGA[dinfo_TCGA$ERG_status!="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none" ,]
#ETSpos 
dinfo_TCGA_ETSpos<-dinfo_TCGA[dinfo_TCGA$ERG_status=="none" & (dinfo_TCGA$ETV1_status!="none" | dinfo_TCGA$ETV4_status!="none" | dinfo_TCGA$FLI1_status!="none") ,]
#ERGpos + PTEN_del
dinfo_TCGA_ERG_PTENdel<-dinfo_TCGA[(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_ERGpos$Sample.ID)  & dinfo_TCGA$PTEN_CNA!="diploid" ,]
#ERGpos + PTEN_wt
dinfo_TCGA_ERG_PTENwt<-dinfo_TCGA[(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_ERGpos$Sample.ID)  & dinfo_TCGA$PTEN_CNA=="diploid" ,]
#other
dinfo_TCGA_other<-dinfo_TCGA[!(dinfo_TCGA$Sample.ID %in% dinfo_TCGA_CHD1del$Sample.ID | dinfo_TCGA$Sample.ID %in% dinfo_TCGA_SPOPmut$Sample.ID 
                               | dinfo_TCGA$Sample.ID %in% dinfo_TCGA_ERGpos$Sample.ID | dinfo_TCGA$Sample.ID %in% dinfo_TCGA_ETSpos$Sample.ID),]

#Read count input downloaded from GDC TCGA PCA RNA-seq data
dinfo_GDC<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//HTSeq/gdc_sample_sheet.2018-03-09.tsv", sep="\t", header=T, check.names=F)
dinfo_GDC_CHD1del<-merge(dinfo_GDC, dinfo_TCGA_CHD1del, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_SPOPmut<-merge(dinfo_GDC, dinfo_TCGA_SPOPmut, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_SPOPmut_CHD1del<-merge(dinfo_GDC, dinfo_TCGA_SPOPmut_CHD1del, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_SPOPmut_CHD1wt<-merge(dinfo_GDC, dinfo_TCGA_SPOPmut_CHD1wt, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_SPOPwt_CHD1del<-merge(dinfo_GDC, dinfo_TCGA_SPOPwt_CHD1del, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_ERGpos<-merge(dinfo_GDC, dinfo_TCGA_ERGpos, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_ETSpos<-merge(dinfo_GDC, dinfo_TCGA_ETSpos, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_other<-merge(dinfo_GDC, data.frame(dinfo_TCGA_other), by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_N<-merge(dinfo_GDC, data.frame(dinfo_HTSeq_n), by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_ERGpos_PTENdel<-merge(dinfo_GDC, dinfo_TCGA_ERG_PTENdel, by.x="Sample ID", by.y="Sample.ID")
dinfo_GDC_ERGpos_PTENwt<-merge(dinfo_GDC, dinfo_TCGA_ERG_PTENwt, by.x="Sample ID", by.y="Sample.ID")

setwd("/Users/deli/Desktop/Chris/TCGA-PCA//HTSeq/")
file_list <- list.files("/Users/deli/Desktop/Chris/TCGA-PCA//HTSeq/")

file_list_CHD1del <- file_list[grepl(paste(dinfo_GDC_CHD1del$`File Name`,collapse="|"), file_list)]
file_list_SPOPmut <- file_list[grepl(paste(dinfo_GDC_SPOPmut$`File Name`,collapse="|"), file_list)]
file_list_ERGpos <- file_list[grepl(paste(dinfo_GDC_ERGpos$`File Name`,collapse="|"), file_list)]
file_list_ETSpos <- file_list[grepl(paste(dinfo_GDC_ETSpos$`File Name`,collapse="|"), file_list)]
file_list_other <- file_list[grepl(paste(dinfo_GDC_other$`File Name`,collapse="|"), file_list)]
file_list_n <- file_list[grepl(paste(dinfo_HTSeq_n[,2],collapse="|"), file_list)]
file_list_ERGpos_PTENdel <- file_list[grepl(paste(dinfo_GDC_ERGpos_PTENdel$`File Name`,collapse="|"), file_list)]
file_list_ERGpos_PTENwt <- file_list[grepl(paste(dinfo_GDC_ERGpos_PTENwt$`File Name`,collapse="|"), file_list)]
file_list_SPOPmut_CHD1wt <- file_list[grepl(paste(dinfo_GDC_SPOPmut_CHD1wt$`File Name`,collapse="|"), file_list)]
file_list_SPOPwt_CHD1del <- file_list[grepl(paste(dinfo_GDC_SPOPwt_CHD1del$`File Name`,collapse="|"), file_list)]
file_list_SPOPmut_CHD1del <- file_list[grepl(paste(dinfo_GDC_SPOPmut_CHD1del$`File Name`,collapse="|"), file_list)]

#Read read count files
dataset_CHD1del<-NULL
for (i in 1:length(file_list_CHD1del)){ 
  if (i==1){
    dataset_CHD1del <- read.table(gzfile(file_list_CHD1del[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_CHD1del[i]), header=F)
    dataset_CHD1del<-cbind(dataset_CHD1del, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_SPOPmut<-NULL
for (i in 1:length(file_list_SPOPmut)){ 
  if (i==1){
    dataset_SPOPmut <- read.table(gzfile(file_list_SPOPmut[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_SPOPmut[i]), header=F)
    dataset_SPOPmut<-cbind(dataset_SPOPmut, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_ERGpos<-NULL
for (i in 1:length(file_list_ERGpos)){ 
  if (i==1){
    dataset_ERGpos <- read.table(gzfile(file_list_ERGpos[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_ERGpos[i]), header=F)
    dataset_ERGpos<-cbind(dataset_ERGpos, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_ERGpos_PTENdel<-NULL
for (i in 1:length(file_list_ERGpos_PTENdel)){ 
  if (i==1){
    dataset_ERGpos_PTENdel <- read.table(gzfile(file_list_ERGpos_PTENdel[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_ERGpos_PTENdel[i]), header=F)
    dataset_ERGpos_PTENdel<-cbind(dataset_ERGpos_PTENdel, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_ERGpos_PTENwt<-NULL
for (i in 1:length(file_list_ERGpos)){ 
  if (i==1){
    dataset_ERGpos_PTENwt <- read.table(gzfile(file_list_ERGpos_PTENwt[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_ERGpos_PTENwt[i]), header=F)
    dataset_ERGpos_PTENwt<-cbind(dataset_ERGpos_PTENwt, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_ETSpos<-NULL
for (i in 1:length(file_list_ETSpos)){ 
  if (i==1){
    dataset_ETSpos <- read.table(gzfile(file_list_ETSpos[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_ETSpos[i]), header=F)
    dataset_ETSpos<-cbind(dataset_ETSpos, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_other<-NULL
for (i in 1:length(file_list_other)){ 
  if (i==1){
    dataset_other <- read.table(gzfile(file_list_other[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_other[i]), header=F)
    dataset_other<-cbind(dataset_other, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_Normal<-NULL
for (i in 1:length(file_list_n)){ 
  if (i==1){
    dataset_Normal <- read.table(gzfile(file_list_n[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_n[i]), header=F)
    dataset_Normal<-cbind(dataset_Normal, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_SPOPmut_CHD1wt<-NULL
for (i in 1:length(file_list_SPOPmut_CHD1wt)){ 
  if (i==1){
    dataset_SPOPmut_CHD1wt <- read.table(gzfile(file_list_SPOPmut_CHD1wt[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_SPOPmut_CHD1wt[i]), header=F)
    dataset_SPOPmut_CHD1wt<-cbind(dataset_SPOPmut_CHD1wt, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_SPOPmut_CHD1del<-NULL
for (i in 1:length(file_list_SPOPmut_CHD1del)){ 
  if (i==1){
    dataset_SPOPmut_CHD1del <- read.table(gzfile(file_list_SPOPmut_CHD1del[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_SPOPmut_CHD1del[i]), header=F)
    dataset_SPOPmut_CHD1del<-cbind(dataset_SPOPmut_CHD1del, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_SPOPwt_CHD1del<-NULL
for (i in 1:length(file_list_SPOPmut_CHD1del)){ 
  if (i==1){
    dataset_SPOPwt_CHD1del <- read.table(gzfile(file_list_SPOPwt_CHD1del[i]), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(file_list_SPOPwt_CHD1del[i]), header=F)
    dataset_SPOPwt_CHD1del<-cbind(dataset_SPOPwt_CHD1del, temp_dataset)
    temp_dataset<-NULL
  }
}

dataset_CHD1del<-dataset_CHD1del[,c(1,seq(2,ncol(dataset_CHD1del),by=2))]
dataset_SPOPmut<-dataset_SPOPmut[,c(1,seq(2,ncol(dataset_SPOPmut),by=2))]
dataset_ERGpos<-dataset_ERGpos[,c(1,seq(2,ncol(dataset_ERGpos),by=2))]
dataset_ETSpos<-dataset_ETSpos[,c(1,seq(2,ncol(dataset_ETSpos),by=2))]
dataset_other<-dataset_other[,c(1,seq(2,ncol(dataset_other),by=2))]
dataset_Normal<-dataset_Normal[,c(1,seq(2,ncol(dataset_Normal),by=2))]
dataset_ERGpos_PTENdel<-dataset_ERGpos_PTENdel[,c(1,seq(2,ncol(dataset_ERGpos_PTENdel),by=2))]
dataset_ERGpos_PTENwt<-dataset_ERGpos_PTENwt[,c(1,seq(2,ncol(dataset_ERGpos_PTENwt),by=2))]
dataset_SPOPmut_CHD1wt<-dataset_SPOPmut_CHD1wt[,c(1,seq(2,ncol(dataset_SPOPmut_CHD1wt),by=2))]
dataset_SPOPwt_CHD1del<-dataset_SPOPwt_CHD1del[,c(1,seq(2,ncol(dataset_SPOPwt_CHD1del),by=2))]
dataset_SPOPmut_CHD1del<-dataset_SPOPmut_CHD1del[,c(1,seq(2,ncol(dataset_SPOPmut_CHD1del),by=2))]

colnames(dataset_CHD1del)<-c("ENSG", as.matrix(merge(dinfo_GDC_CHD1del, data.frame(file_list_CHD1del), by.x="File Name", by.y="file_list_CHD1del")[,2]))
colnames(dataset_SPOPmut)<-c("ENSG", as.matrix(merge(dinfo_GDC_SPOPmut, data.frame(file_list_SPOPmut), by.x="File Name", by.y="file_list_SPOPmut")[,2]))
colnames(dataset_ERGpos)<-c("ENSG", as.matrix(merge(dinfo_GDC_ERGpos, data.frame(file_list_ERGpos), by.x="File Name", by.y="file_list_ERGpos")[,2]))
colnames(dataset_ETSpos)<-c("ENSG", as.matrix(merge(dinfo_GDC_ETSpos, data.frame(file_list_ETSpos), by.x="File Name", by.y="file_list_ETSpos")[,2]))
colnames(dataset_other)<-c("ENSG", as.matrix(merge(dinfo_GDC_other, data.frame(file_list_other), by.x="File Name", by.y="file_list_other")[,2]))
colnames(dataset_Normal)<-c("ENSG", as.matrix(merge(dinfo_HTSeq_n, data.frame(file_list_n), by.x="File.Name", by.y="file_list_n")$Sample.ID))
colnames(dataset_ERGpos_PTENdel)<-c("ENSG", as.matrix(merge(dinfo_GDC_ERGpos_PTENdel, data.frame(file_list_ERGpos_PTENdel), by.x="File Name", by.y="file_list_ERGpos_PTENdel")[,2]))
colnames(dataset_ERGpos_PTENwt)<-c("ENSG", as.matrix(merge(dinfo_GDC_ERGpos_PTENwt, data.frame(file_list_ERGpos_PTENwt), by.x="File Name", by.y="file_list_ERGpos_PTENwt")[,2]))
colnames(dataset_SPOPmut_CHD1wt)<-c("ENSG", as.matrix(merge(dinfo_GDC_SPOPmut_CHD1wt, data.frame(file_list_SPOPmut_CHD1wt), by.x="File Name", by.y="file_list_SPOPmut_CHD1wt")[,2]))
colnames(dataset_SPOPmut_CHD1del)<-c("ENSG", as.matrix(merge(dinfo_GDC_SPOPmut_CHD1del, data.frame(file_list_SPOPmut_CHD1del), by.x="File Name", by.y="file_list_SPOPmut_CHD1del")[,2]))
colnames(dataset_SPOPwt_CHD1del)<-c("ENSG", as.matrix(merge(dinfo_GDC_SPOPwt_CHD1del, data.frame(file_list_SPOPwt_CHD1del), by.x="File Name", by.y="file_list_SPOPwt_CHD1del")[,2]))


#Add gene symbol
dgene<-read.table("/Users/deli/Desktop/Chris/annotation/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gene", sep=" ", header=F)
d_ENSG<-NULL
for (i in 1:nrow(dgene))
{
  d_ENSG<-c(d_ENSG, strsplit(as.character(dgene[i,1]), ".", 2)[[1]][1])
}
dgene$ENSG<-d_ENSG

dataset_CHD1del_ENSG<-NULL
for (i in 1:nrow(dataset_CHD1del))
{
  dataset_CHD1del_ENSG<-c(dataset_CHD1del_ENSG, strsplit(as.character(dataset_CHD1del[i,1]), ".", 2)[[1]][1])
}
dataset_CHD1del$ENSG<-dataset_CHD1del_ENSG

dataset_SPOPmut_ENSG<-NULL
for (i in 1:nrow(dataset_SPOPmut))
{
  dataset_SPOPmut_ENSG<-c(dataset_SPOPmut_ENSG, strsplit(as.character(dataset_SPOPmut[i,1]), ".", 2)[[1]][1])
}
dataset_SPOPmut$ENSG<-dataset_SPOPmut_ENSG

dataset_SPOPmut_CHD1wt_ENSG<-NULL
for (i in 1:nrow(dataset_SPOPmut_CHD1wt))
{
  dataset_SPOPmut_CHD1wt_ENSG<-c(dataset_SPOPmut_CHD1wt_ENSG, strsplit(as.character(dataset_SPOPmut_CHD1wt[i,1]), ".", 2)[[1]][1])
}
dataset_SPOPmut_CHD1wt$ENSG<-dataset_SPOPmut_CHD1wt_ENSG

dataset_SPOPmut_CHD1del_ENSG<-NULL
for (i in 1:nrow(dataset_SPOPmut_CHD1del))
{
  dataset_SPOPmut_CHD1del_ENSG<-c(dataset_SPOPmut_CHD1del_ENSG, strsplit(as.character(dataset_SPOPmut_CHD1del[i,1]), ".", 2)[[1]][1])
}
dataset_SPOPmut_CHD1del$ENSG<-dataset_SPOPmut_CHD1del_ENSG

dataset_SPOPwt_CHD1del_ENSG<-NULL
for (i in 1:nrow(dataset_SPOPwt_CHD1del))
{
  dataset_SPOPwt_CHD1del_ENSG<-c(dataset_SPOPwt_CHD1del_ENSG, strsplit(as.character(dataset_SPOPwt_CHD1del[i,1]), ".", 2)[[1]][1])
}
dataset_SPOPwt_CHD1del$ENSG<-dataset_SPOPwt_CHD1del_ENSG

dataset_ERGpos_ENSG<-NULL
for (i in 1:nrow(dataset_ERGpos))
{
  dataset_ERGpos_ENSG<-c(dataset_ERGpos_ENSG, strsplit(as.character(dataset_ERGpos[i,1]), ".", 2)[[1]][1])
}
dataset_ERGpos$ENSG<-dataset_ERGpos_ENSG

dataset_ERGpos_PTENdel_ENSG<-NULL
for (i in 1:nrow(dataset_ERGpos_PTENdel))
{
  dataset_ERGpos_PTENdel_ENSG<-c(dataset_ERGpos_PTENdel_ENSG, strsplit(as.character(dataset_ERGpos_PTENdel[i,1]), ".", 2)[[1]][1])
}
dataset_ERGpos_PTENdel$ENSG<-dataset_ERGpos_PTENdel_ENSG

dataset_ERGpos_PTENwt_ENSG<-NULL
for (i in 1:nrow(dataset_ERGpos))
{
  dataset_ERGpos_PTENwt_ENSG<-c(dataset_ERGpos_PTENwt_ENSG, strsplit(as.character(dataset_ERGpos_PTENwt[i,1]), ".", 2)[[1]][1])
}
dataset_ERGpos_PTENwt$ENSG<-dataset_ERGpos_PTENwt_ENSG

dataset_ETSpos_ENSG<-NULL
for (i in 1:nrow(dataset_ETSpos))
{
  dataset_ETSpos_ENSG<-c(dataset_ETSpos_ENSG, strsplit(as.character(dataset_ETSpos[i,1]), ".", 2)[[1]][1])
}
dataset_ETSpos$ENSG<-dataset_ETSpos_ENSG

dataset_other_ENSG<-NULL
for (i in 1:nrow(dataset_other))
{
  dataset_other_ENSG<-c(dataset_other_ENSG, strsplit(as.character(dataset_other[i,1]), ".", 2)[[1]][1])
}
dataset_other$ENSG<-dataset_other_ENSG

dataset_Normal_ENSG<-NULL
for (i in 1:nrow(dataset_Normal))
{
  dataset_Normal_ENSG<-c(dataset_Normal_ENSG, strsplit(as.character(dataset_Normal[i,1]), ".", 2)[[1]][1])
}
dataset_Normal$ENSG<-dataset_Normal_ENSG


dataset_CHD1del_v2<-merge(dgene[,2:3], dataset_CHD1del, by="ENSG")
dataset_SPOPmut_v2<-merge(dgene[,2:3], dataset_SPOPmut, by="ENSG")
dataset_SPOPmut_CHD1wt_v2<-merge(dgene[,2:3], dataset_SPOPmut_CHD1wt, by="ENSG")
dataset_SPOPwt_CHD1del_v2<-merge(dgene[,2:3], dataset_SPOPwt_CHD1del, by="ENSG")
dataset_SPOPmut_CHD1del_v2<-merge(dgene[,2:3], dataset_SPOPmut_CHD1del, by="ENSG")
dataset_ERGpos_v2<-merge(dgene[,2:3], dataset_ERGpos, by="ENSG")
dataset_ERGpos_PTENdel_v2<-merge(dgene[,2:3], dataset_ERGpos_PTENdel, by="ENSG")
dataset_ERGpos_PTENwt_v2<-merge(dgene[,2:3], dataset_ERGpos_PTENwt, by="ENSG")
dataset_ETSpos_v2<-merge(dgene[,2:3], dataset_ETSpos, by="ENSG")
dataset_other_v2<-merge(dgene[,2:3], dataset_other, by="ENSG")
dataset_Normal_v2<-merge(dgene[,2:3], dataset_Normal, by="ENSG")

#GENCODE gene type
dgene_type<-read.table("/Users/deli/Desktop/Chris/annotation/human/gencode.v19.annotation.gtf.genes_position_type", sep=" ", header=F, check.names = F)
dgene_type<-merge(dgene, dgene_type, by.x="V1", by.y="V5", sort=F)


#Figure 2A: heatmap of two distinct tumor lineage models of PCa progression: ERG/PTEN (N->ERG->PTEN) and SPOP/CHD1 (N->POP->CHD1) via ImpulseDE2
#Combine Normal, SPOPmut_CHD1wt, SPOPmut_CHD1del samples
dataset<-merge(merge(dataset_Normal_v2, dataset_SPOPmut_CHD1wt_v2[,c(1,3:ncol(dataset_SPOPmut_CHD1wt_v2))], by="ENSG"), 
               dataset_SPOPmut_CHD1del_v2[,c(1,3:ncol(dataset_SPOPmut_CHD1del_v2))], by="ENSG")
dataset_coding<-merge(dataset, dgene_type[dgene_type$V7=="protein_coding",], by.x="ENSG", by.y="ENSG")
dataset_coding_count<-dataset_coding[,3:(ncol(dataset_coding)-9)]
rownames(dataset_coding_count)<-dataset_coding[,1]
filtered <- dataset_coding_count[apply(dataset_coding_count[,1:ncol(dataset_coding_count)],1,mean)>=1,]
#dataset_count<-dataset[,3:ncol(dataset)]
#rownames(dataset_count)<-dataset[,1]
#filtered <- dataset_count[apply(dataset_count[,1:ncol(dataset_count)],1,mean)>=1,]
x <- as.factor(c(rep("Normal",nrow(dinfo_GDC_N)), rep("SPOPmut_CHD1wt",nrow(dinfo_GDC_SPOPmut_CHD1wt)), rep("SPOPmut_CHD1del",nrow(dinfo_GDC_SPOPmut_CHD1del))))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set1 <- betweenLaneNormalization(set, which="upper")
GeneData<-normCounts(set1)
#ImpulseDE2 test
GeneData_N_SPOP_CHD1<-GeneData[,c(1:nrow(dinfo_GDC_N), (nrow(dinfo_GDC_N)+1):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1wt)-1), 
                                  ((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1wt)):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1wt)+nrow(dinfo_GDC_SPOPmut_CHD1del)-1))]
CondVector<- c(rep("Normal",nrow(dinfo_GDC_N)), rep("SPOPmut_CHD1wt",nrow(dinfo_GDC_SPOPmut_CHD1wt)), rep("SPOPmut_CHD1del",nrow(dinfo_GDC_SPOPmut_CHD1del)))
Time <- factor(CondVector, levels=c("Normal","SPOPmut_CHD1wt","SPOPmut_CHD1del"))
annotation<-data.frame(Time)
rownames(annotation)<-colnames(GeneData_N_SPOP_CHD1)
annotation$Sample<-rownames(annotation)
annotation$Condition<-"PRCA"
annotation$Time <- as.numeric(annotation$Time)
annotation$Condition <- "case"
GeneData_N_SPOP_CHD1_impulse2_results <- runImpulseDE2(
                       matCountData    = GeneData_N_SPOP_CHD1, 
                       dfAnnotation    = annotation,
                       boolCaseCtrl    = FALSE,
                       vecConfounders  = NULL,
                       boolIdentifyTransients = TRUE, 
                       scaNProc        = 1 )
FDR_cutoff<-1e-10
pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_", FDR_cutoff, "-hm.pdf", sep=""), width=12, height=12)
GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps <- plotHeatmap(GeneData_N_SPOP_CHD1_impulse2_results,
                                                                strCondition           = "case",
                                                                boolIdentifyTransients = TRUE,
                                                                scaQThres              = FDR_cutoff)
draw(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$complexHeatmapRaw)
dev.off()
#GeneData_N_SPOP_CHD1_NormConst<-computeNormConst(matCountData    = GeneData_N_SPOP_CHD1)
#Genes in each pattern
length(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up)
length(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down)
length(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up)
length(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_down)
#Output protein_coding genes
write.csv(merge(data.frame(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up), dgene_type[dgene_type$V7=="protein_coding",], 
                by.y="ENSG", by.x="GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_up", sort=F), 
            file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_", FDR_cutoff, "-transition_up_coding.csv", sep=""), row.names=F)
write.csv(merge(data.frame(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down), dgene_type[dgene_type$V7=="protein_coding",], 
                  by.y="ENSG", by.x="GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_down", sort=F), 
            file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_", FDR_cutoff, "-transition_down_coding.csv", sep=""), row.names=F)
write.csv(merge(data.frame(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up), dgene_type[dgene_type$V7=="protein_coding",], 
                  by.y="ENSG", by.x="GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps.lsvecGeneGroups.transient_up", sort=F), 
            file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_", FDR_cutoff, "-transient_up_coding.csv", sep=""), row.names=F)


#Combine Normal, ERG_PTENwt, ERG_PTENdel samples
dataset<-merge(merge(dataset_Normal_v2, dataset_ERGpos_PTENwt_v2[,c(1,3:ncol(dataset_ERGpos_PTENwt_v2))], by="ENSG"), 
               dataset_ERGpos_PTENdel_v2[,c(1,3:ncol(dataset_ERGpos_PTENdel_v2))], by="ENSG")
dataset_coding<-merge(dataset, dgene_type[dgene_type$V7=="protein_coding",], by.x="ENSG", by.y="ENSG")
dataset_coding_count<-dataset_coding[,3:(ncol(dataset_coding)-9)]
rownames(dataset_coding_count)<-dataset_coding[,1]
filtered <- dataset_coding_count[apply(dataset_coding_count[,1:ncol(dataset_coding_count)],1,mean)>=1,]
#dataset_count<-dataset[,3:ncol(dataset)]
#rownames(dataset_count)<-dataset[,1]
#filtered <- dataset_count[apply(dataset_count[,1:ncol(dataset_count)],1,mean)>=1,]
x <- as.factor(c(rep("Normal",nrow(dinfo_GDC_N)), rep("ERG_PTENwt",nrow(dinfo_GDC_ERGpos_PTENwt)), rep("ERG_PTENdel",nrow(dinfo_GDC_ERGpos_PTENdel))))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
#set
set1 <- betweenLaneNormalization(set, which="upper")
GeneData<-normCounts(set1)
#set1_normCounts<-merge(normCounts(set1), dgene[,2:3], by.y="ENSG", by.x="row.names", sort=F)
#ImpulseDE2 test
#GeneData_N_ERG_PTEN<-GeneData[,c(1:52, 53:(53+88-1), (53+88):(53+88+60-1))]
GeneData_N_ERG_PTEN<-GeneData[,c(1:nrow(dinfo_GDC_N), (nrow(dinfo_GDC_N)+1):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENwt)-1), 
                                 ((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENwt)):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENwt)+nrow(dinfo_GDC_ERGpos_PTENdel)-1))]
CondVector<- c(rep("0",nrow(dinfo_GDC_N)), rep("1",nrow(dinfo_GDC_ERGpos_PTENwt)), rep("2",nrow(dinfo_GDC_ERGpos_PTENdel)))
Time <- factor(CondVector, levels=c("0","1","2"))
annotation<-data.frame(Time)
rownames(annotation)<-colnames(GeneData_N_ERG_PTEN)
annotation$Sample<-rownames(annotation)
annotation$Condition<-"PRCA"
annotation$Time <- as.numeric(annotation$Time)
annotation$Condition <- "case"
GeneData_N_ERG_PTEN_impulse2_results <- runImpulseDE2(
  matCountData    = GeneData_N_ERG_PTEN, 
  dfAnnotation    = annotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  boolIdentifyTransients = TRUE, 
  scaNProc        = 1 )
FDR_cutoff<-1e-10
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_", FDR_cutoff, "-hm.pdf", sep=""), width=12, height=12)
GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps <- plotHeatmap(GeneData_N_ERG_PTEN_impulse2_results,
                                                               strCondition           = "case",
                                                               boolIdentifyTransients = TRUE,
                                                               scaQThres              = FDR_cutoff)
draw(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$complexHeatmapRaw)
dev.off()
#Genes in each pattern
length(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up)
length(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down)
length(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up)
length(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_down)
#Output protein_coding genes
write.csv(merge(data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up), dgene_type[dgene_type$V7=="protein_coding",], 
                by.y="ENSG", by.x="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_up", sort=F), 
          file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_", FDR_cutoff, "-transition_up_coding.csv", sep=""), row.names=F)
write.csv(merge(data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down), dgene_type[dgene_type$V7=="protein_coding",], 
                by.y="ENSG", by.x="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_down", sort=F), 
          file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_", FDR_cutoff, "-transition_down_coding.csv", sep=""), row.names=F)
write.csv(merge(data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up), dgene_type[dgene_type$V7=="protein_coding",], 
                by.y="ENSG", by.x="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transient_up", sort=F), 
          file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_", FDR_cutoff, "-transient_up_coding.csv", sep=""), row.names=F)
write.csv(merge(data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_down), dgene_type[dgene_type$V7=="protein_coding",], 
                by.y="ENSG", by.x="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transient_down", sort=F), 
          file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_", FDR_cutoff, "-transient_down_coding.csv", sep=""), row.names=F)



#Figure 2A: Gene number: two distinct tumor lineage models of PCa progression: ERG/PTEN (N->ERG->PTEN) and SPOP/CHD1 (N->SPOP->CHD1) via ImpulseDE2
#Impulse2 signature 
GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_up_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transition_up_coding.csv")
GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_down_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transition_down_coding.csv")
GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transient_up_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transient_up_coding.csv")
#GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transient_down_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transient_down_coding.cs")

GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_up_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_1e-10-transition_up_coding.csv")
GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_down_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_1e-10-transition_down_coding.csv")
GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transient_up_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_1e-10-transient_up_coding.csv")
GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transient_down_coding<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_ERG_PTENdel_unique-ImpulseDE2_qval_1e-10-transient_down_coding.csv")

GeneData_N_SPOP_CHD1_impulse2_results_count<-data.frame(c(nrow(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_up_coding),
                                                          nrow(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_down_coding),
                                                          nrow(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transient_up_coding),
                                                          0))
#rownames(GeneData_N_SPOP_CHD1_impulse2_results_count) <- c("transition_up", "transition_down", "transient_up", "transient_down")
colnames(GeneData_N_SPOP_CHD1_impulse2_results_count) <- "Impulse2_results_count"
GeneData_N_ERG_PTEN_impulse2_results_count<-data.frame(c(nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_up_coding),
                                                         nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_down_coding),
                                                         nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transient_up_coding),
                                                         nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transient_down_coding)))
#rownames(GeneData_N_ERG_PTEN_impulse2_results_count) <- c("transition_up", "transition_down", "transient_up", "transient_down")
colnames(GeneData_N_ERG_PTEN_impulse2_results_count) <- "Impulse2_results_count"

Impulse2_results_count<-rbind(GeneData_N_ERG_PTEN_impulse2_results_count, GeneData_N_SPOP_CHD1_impulse2_results_count)
Impulse2_results_count$Input<-c(rep("N_ERG_PTEN", 4), rep("N_SPOP_CHD1", 4))
Impulse2_results_count$Type<-rep(c("transition_up", "transition_down", "transient_up", "transient_down"),2)
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-ImpulseDE2_gene_count_v2.pdf", width=6,height=6)
ggplot(aes(x = Type, y = Impulse2_results_count, fill=Input),  data = Impulse2_results_count) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values = c(brewer.pal(9,"Set1")[1], brewer.pal(9,"Set1")[5])) +
  #coord_flip() +
  theme(text = element_text(size=12))
dev.off()


#Figure 2B: Venn diagrams of shared and uniquely upregulated and downregulated genes between two tumor lineage models.  Number of shared and unique altered genes were annotated.
#Transition up signature comparison between ERG_PTEN and SPOP_CHD1
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transition_up_coding_merge_ERG_PTEN_gene_venn.pdf", width=6, height=6)
venn.plot <- draw.pairwise.venn(area1      = nrow(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_up_coding),
                                area2      = nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_up_coding),
                                cross.area = nrow(merge(data.frame(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_up_coding),
                                                        data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_up_coding),
                                                        by.x="GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_up",
                                                        by.y="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_up", sort=F)),
                                category   = c("N_SPOP_CHD1_transition_up", "N_ERG_PTEN_transition_up"),
                                fill            = c(col_ERG, col_SPOP),
                                #fill = c("light blue", "pink"),
                                scaled     = TRUE)
dev.off()

#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_SPOP_CHD1-ImpulseDE2_qval_1e-10-transition_down_coding_merge_ERG_PTEN_gene_venn.pdf", width=6, height=6)
venn.plot <- draw.pairwise.venn(area1      = nrow(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_down_coding),
                                area2      = nrow(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_down_coding),
                                cross.area = nrow(merge(data.frame(GeneData_N_SPOP_CHD1_impulse2_results_lsHeatmaps_transition_down_coding),
                                                        data.frame(GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps_transition_down_coding),
                                                        by.x="Ensembl_ID",
                                                        by.y="GeneData_N_ERG_PTEN_impulse2_results_lsHeatmaps.lsvecGeneGroups.transition_down", sort=F)),
                                category   = c("N_SPOP_CHD1_transition_down", "N_ERG_PTEN_transition_down"),
                                fill            = c(col_ERG, col_SPOP),
                                scaled     = TRUE)
dev.off()




#Figure S2. Number of transiently and progressively upregulated and downregulated genes from normal to PTEN deletion and ERG fusion (red), or from normal to CHD1 deletion and SPOP mutation (orange).
#ImpulseDE2 test: Normal, SPOPmut_CHD1del, SPOPmut_CHD1wt samples
dataset<-merge(merge(dataset_Normal_v2, dataset_SPOPmut_CHD1del_v2[,c(1,3:ncol(dataset_SPOPmut_CHD1del_v2))], by="ENSG"), 
               dataset_SPOPmut_CHD1wt_v2[,c(1,3:ncol(dataset_SPOPmut_CHD1wt_v2))], by="ENSG")
dataset_coding<-merge(dataset, dgene_type[dgene_type$V7=="protein_coding",], by.x="ENSG", by.y="ENSG")
dataset_coding_count<-dataset_coding[,3:(ncol(dataset_coding)-9)]
rownames(dataset_coding_count)<-dataset_coding[,1]
filtered <- dataset_coding_count[apply(dataset_coding_count[,1:ncol(dataset_coding_count)],1,mean)>=1,]
#dataset_count<-dataset[,3:ncol(dataset)]
#rownames(dataset_count)<-dataset[,1]
#filtered <- dataset_count[apply(dataset_count[,1:ncol(dataset_count)],1,mean)>=1,]
x <- as.factor(c(rep("Normal",nrow(dinfo_GDC_N)), rep("SPOPmut_CHD1wt",nrow(dinfo_GDC_SPOPmut_CHD1del)), rep("SPOPmut_CHD1del",nrow(dinfo_GDC_SPOPmut_CHD1wt))))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set1 <- betweenLaneNormalization(set, which="upper")
GeneData<-normCounts(set1)
GeneData_N_CHD1_SPOP<-GeneData[,c(1:nrow(dinfo_GDC_N), (nrow(dinfo_GDC_N)+1):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1del)-1), 
                                  ((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1del)):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_SPOPmut_CHD1del)+nrow(dinfo_GDC_SPOPmut_CHD1wt)-1))]
CondVector<- c(rep("Normal",nrow(dinfo_GDC_N)), rep("SPOPmut_CHD1del",nrow(dinfo_GDC_SPOPmut_CHD1wt)), rep("SPOPmut_CHD1wt",nrow(dinfo_GDC_SPOPmut_CHD1del)))
Time <- factor(CondVector, levels=c("Normal","SPOPmut_CHD1del","SPOPmut_CHD1wt"))
annotation<-data.frame(Time)
rownames(annotation)<-colnames(GeneData_N_CHD1_SPOP)
annotation$Sample<-rownames(annotation)
annotation$Condition<-"PRCA"
annotation$Time <- as.numeric(annotation$Time)
annotation$Condition <- "case"
GeneData_N_CHD1_SPOP_impulse2_results <- runImpulseDE2(
  matCountData    = GeneData_N_CHD1_SPOP, 
  dfAnnotation    = annotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  boolIdentifyTransients = TRUE, 
  scaNProc        = 1 )
FDR_cutoff<-1e-10
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-coding-N_CHD1_SPOP-ImpulseDE2_qval_", FDR_cutoff, "-hm.pdf", sep=""), width=12, height=12)
GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps <- plotHeatmap(GeneData_N_CHD1_SPOP_impulse2_results,
                                                                strCondition           = "case",
                                                                boolIdentifyTransients = TRUE,
                                                                scaQThres              = FDR_cutoff)
draw(GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps$complexHeatmapRaw)
dev.off()
#GeneData_N_SPOP_CHD1_NormConst<-computeNormConst(matCountData    = GeneData_N_CHD1_SPOP)
#Genes in each pattern
length(GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up)
length(GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down)
length(GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up)
length(GeneData_N_CHD1_SPOP_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_down)


#ImpulseDE2 test: Normal,ERG_PTENdel, ERG_PTENwt samples
dataset<-merge(merge(dataset_Normal_v2, dataset_ERGpos_PTENdel_v2[,c(1,3:ncol(dataset_ERGpos_PTENdel_v2))], by="ENSG"), 
               dataset_ERGpos_PTENwt_v2[,c(1,3:ncol(dataset_ERGpos_PTENwt_v2))], by="ENSG")
dataset_coding<-merge(dataset, dgene_type[dgene_type$V7=="protein_coding",], by.x="ENSG", by.y="ENSG")
dataset_coding_count<-dataset_coding[,3:(ncol(dataset_coding)-9)]
rownames(dataset_coding_count)<-dataset_coding[,1]
filtered <- dataset_coding_count[apply(dataset_coding_count[,1:ncol(dataset_coding_count)],1,mean)>=1,]
#dataset_count<-dataset[,3:ncol(dataset)]
#rownames(dataset_count)<-dataset[,1]
#filtered <- dataset_count[apply(dataset_count[,1:ncol(dataset_count)],1,mean)>=1,]
x <- as.factor(c(rep("Normal",nrow(dinfo_GDC_N)), rep("ERG_PTENdel",nrow(dinfo_GDC_ERGpos_PTENdel)), rep("ERG_PTENwt",nrow(dinfo_GDC_ERGpos_PTENwt))))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
#set
set1 <- betweenLaneNormalization(set, which="upper")
GeneData<-normCounts(set1)
GeneData_N_PTEN_ERG<-GeneData[,c(1:nrow(dinfo_GDC_N), (nrow(dinfo_GDC_N)+1):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENdel)-1), 
                                 ((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENdel)):((nrow(dinfo_GDC_N)+1)+nrow(dinfo_GDC_ERGpos_PTENdel)+nrow(dinfo_GDC_ERGpos_PTENwt)-1))]
CondVector<- c(rep("0",nrow(dinfo_GDC_N)), rep("1",nrow(dinfo_GDC_ERGpos_PTENdel)), rep("2",nrow(dinfo_GDC_ERGpos_PTENwt)))
Time <- factor(CondVector, levels=c("0","1","2"))
annotation<-data.frame(Time)
rownames(annotation)<-colnames(GeneData_N_PTEN_ERG)
annotation$Sample<-rownames(annotation)
annotation$Condition<-"PRCA"
annotation$Time <- as.numeric(annotation$Time)
annotation$Condition <- "case"
GeneData_N_PTEN_ERG_impulse2_results <- runImpulseDE2(
  matCountData    = GeneData_N_PTEN_ERG, 
  dfAnnotation    = annotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  boolIdentifyTransients = TRUE, 
  scaNProc        = 1 )
FDR_cutoff<-1e-10
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-all-ImpulseDE2/TCGA-freeze-333-coding-N_ERGpPTENdel_ERGpPTENwt_unique-ImpulseDE2_qval_", FDR_cutoff, "-hm.pdf", sep=""), width=12, height=12)
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-N_ERGpPTENdel_ERGpPTENwt_unique-ImpulseDE2_qval_1e-7-hm.pdf", width=12, height=12)
GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps <- plotHeatmap(GeneData_N_PTEN_ERG_impulse2_results,
                                                               strCondition           = "case",
                                                               boolIdentifyTransients = TRUE,
                                                               scaQThres              = FDR_cutoff)
draw(GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps$complexHeatmapRaw)
dev.off()
#Genes in each pattern
length(GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_up)
length(GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps$lsvecGeneGroups$transition_down)
length(GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_up)
length(GeneData_N_PTEN_ERG_impulse2_results_lsHeatmaps$lsvecGeneGroups$transient_down)




