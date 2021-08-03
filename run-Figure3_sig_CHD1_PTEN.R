library(MASS)
library(gdata)
library(RColorBrewer)
library(ggplot2)
library(ggbeeswarm)

#Heatmap color
color2<-colorRampPalette(c("blue", "white", "red"))(256)
#Different cluster method
myclustw <- function(x) {return(hclust(x, method="ward.D") )}
myclustw2 <- function(x) {return(hclust(x, method="ward.D2") )}
myclusts <- function(x) {return(hclust(x, method="single") )}
myclustc <- function(x) {return(hclust(x, method="complete") )}
myclusta <- function(x) {return(hclust(x, method="average") )}
myclustm <- function(x) {return(hclust(x, method="median") )}
myclustmc <- function(x) {return(hclust(x, method="mcquitty") )}
myclustcen <- function(x) {return(hclust(x, method="centroid") )}

#Subclass color
sur_col_ERG<-brewer.pal(9,"Set1")[1]
sur_col_ERG_PTENdel<-brewer.pal(9,"Set1")[1]
sur_col_ERG_PTENwt<-brewer.pal(12,"Set3")[8]
sur_col_ETS<-brewer.pal(9,"Set1")[6]
sur_col_ETS_PTENdel<-brewer.pal(9,"Set1")[6]
sur_col_ETS_PTENwt<-brewer.pal(12,"Set3")[2]
sur_col_CHD1<-brewer.pal(9,"Set1")[2]
sur_col_other<-brewer.pal(9,"Set1")[3]
sur_col_SPOPo<-brewer.pal(9,"Set1")[5]
sur_col_CHD1o<-brewer.pal(9,"Set1")[2]
sur_col_CHD1_SPOP<-brewer.pal(9,"Set1")[4]

#Subtype color
col_ERG<-brewer.pal(9,"Set1")[8]
col_ETS<-sur_col_ETS
#col_ETS<-brewer.pal(9,"Greens")[5]
col_FOXA1_mut<-brewer.pal(9,"Set1")[8]
col_other<-brewer.pal(12,"Set3")[9]
col_SPOP<-brewer.pal(9,"Set1")[5]
col_CHD1<-brewer.pal(9,"Set1")[2]
col_CHD1homo<-brewer.pal(9,"Set1")[2]
col_CHD1hete<-brewer.pal(12,"Set3")[5]
col_PTEN<-brewer.pal(9,"Set1")[4]

#TCGA PCA 333 sample information from TCGA Cell paper Table S1
dinfo_TCGA<-read.xls("/Users/deli/Dropbox/Papers/TCGA-Prostate/mmc2.xls", check.names=F)

#Tumor ID
dTCGAid<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31.tsv", sep="\t", header=T)
dTCGAid_ID<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/gdc_sample_sheet.2018-05-31-ID.tsv", sep="\t", header=T, check.names=F)
dTCGAid<-cbind(dTCGAid, dTCGAid_ID)
dTCGAid<-dTCGAid[!grepl("TCGA-HC-7740-01B", dTCGAid$Sample.ID),]
dinfo_TCGA<-merge(dinfo_TCGA, dTCGAid, by.x="SAMPLE_ID", by.y="Sample ID")
#dcna_333_CNA<-read.csv("/Users/deli/Desktop/Chris/BPH/SNP6/TCGA-CNA-log2_0.3-fraction.csv")
#dinfo_TCGA_v2<-merge(dinfo_TCGA, dcna_333_CNA[,2:ncol(dcna_333_CNA)], by.x="SAMPLE_ID", by.y="V1", sort=F)


#CHD1 del and ERG/ETS negative
dinfo_TCGA_CHD1del<-dinfo_TCGA[dinfo_TCGA$CHD1_CNA!="diploid" & dinfo_TCGA$ERG_status=="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none",]
#CHD1 wt and ERG/ETS negative
dinfo_TCGA_CHD1wt_ETSneg<-dinfo_TCGA[dinfo_TCGA$CHD1_CNA=="diploid" & dinfo_TCGA$ERG_status=="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none",]
#ERGpos 
dinfo_TCGA_ERGpos<-dinfo_TCGA[dinfo_TCGA$ERG_status!="none" & dinfo_TCGA$ETV1_status=="none" & dinfo_TCGA$ETV4_status=="none" & dinfo_TCGA$FLI1_status=="none" ,]
#ETSpos 
dinfo_TCGA_ETSpos<-dinfo_TCGA[dinfo_TCGA$ERG_status=="none" & (dinfo_TCGA$ETV1_status!="none" | dinfo_TCGA$ETV4_status!="none" | dinfo_TCGA$FLI1_status!="none") ,]
#ETSpos and PTEN deletion
dinfo_TCGA_ETS_PTENdel<-dinfo_TCGA[dinfo_TCGA$PTEN_CNA!="diploid" 
                                   & dinfo_TCGA$ERG_status=="none" & (dinfo_TCGA$ETV1_status!="none" | dinfo_TCGA$ETV4_status!="none" | dinfo_TCGA$FLI1_status!="none") ,]
#ETSpos and PTEN deletion
dinfo_TCGA_ETS_PTENwt<-dinfo_TCGA[dinfo_TCGA$PTEN_CNA=="diploid" 
                                   & dinfo_TCGA$ERG_status=="none" & (dinfo_TCGA$ETV1_status!="none" | dinfo_TCGA$ETV4_status!="none" | dinfo_TCGA$FLI1_status!="none") ,]


#FPKM input from GDC TCGA PCA study
setwd("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/EXP/Download/")
file_list <- list.files("/Users/deli/Desktop/Chris/TCGA-PCA/GDC/EXP/Download/")
dinfo_GDC<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA//HTSeq/gdc_sample_sheet.2018-03-09.tsv", sep="\t", header=T, check.names=F)

dat_FPKM_333<-NULL
for (i in 1:nrow(dinfo_TCGA)){ 
  if (i==1){
    dat_FPKM_333 <- read.table(gzfile(as.matrix(dinfo_TCGA[i,]$File.Name)), header=F)
  }
  if (i>1){
    temp_dataset <-read.table(gzfile(as.matrix(dinfo_TCGA[i,]$File.Name)), header=F)
    dat_FPKM_333<-cbind(dat_FPKM_333, temp_dataset)
    temp_dataset<-NULL
  }
}
dat_FPKM_333<-dat_FPKM_333[,c(1,seq(2,ncol(dat_FPKM_333),by=2))]
colnames(dat_FPKM_333)<-c("ENSG", as.matrix(dinfo_TCGA$Sample.ID))

#Add gene symbol
dgene<-read.table("/Users/deli/Desktop/Chris/annotation/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gene", sep=" ", header=F)
d_ENSG<-NULL
for (i in 1:nrow(dgene))
{
  d_ENSG<-c(d_ENSG, strsplit(as.character(dgene[i,1]), ".", 2)[[1]][1])
}
dgene$ENSG<-d_ENSG

dat_FPKM_333_ENSG<-NULL
for (i in 1:nrow(dat_FPKM_333))
{
  dat_FPKM_333_ENSG<-c(dat_FPKM_333_ENSG, strsplit(as.character(dat_FPKM_333[i,1]), ".", 2)[[1]][1])
}
dat_FPKM_333$ENSG<-dat_FPKM_333_ENSG


dat_FPKM_333_v2<-merge(dgene[,2:3], dat_FPKM_333, by="ENSG")

dat_FPKM_CHD1del_v2 <- dat_FPKM_333_v2[, c(colnames(dat_FPKM_333_v2[, 1:2]), as.matrix(dinfo_TCGA_CHD1del$Sample.ID))]
dat_FPKM_CHD1wt_ETSneg_v2 <- dat_FPKM_333_v2[, c(colnames(dat_FPKM_333_v2[, 1:2]), as.matrix(dinfo_TCGA_CHD1wt_ETSneg$Sample.ID))]

dat_FPKM_ETS_PTENdel_v2 <- dat_FPKM_333_v2[, c(colnames(dat_FPKM_333_v2[, 1:2]), as.matrix(dinfo_TCGA_ETS_PTENdel$Sample.ID))]
dat_FPKM_ETS_PTENwt_v2 <- dat_FPKM_333_v2[, c(colnames(dat_FPKM_333_v2[, 1:2]), as.matrix(dinfo_TCGA_ETS_PTENwt$Sample.ID))]

setwd("/Users/deli/Desktop/Chris/Michael/BRCA/")


##Differentially expressed genes between CHD1del(hetloss+homdel) and CHD1wt_ETSneg
dat<-merge(dat_FPKM_CHD1del_v2, dat_FPKM_CHD1wt_ETSneg_v2[,c(1,3:ncol(dat_FPKM_CHD1wt_ETSneg_v2))], by="ENSG")
dat1<-dat[,c(2,3:ncol(dat))]
dat1$mean_CHD1del<-apply(dat1[,2:(ncol(dat_FPKM_CHD1del_v2)-1)],1,mean)
dat1$mean_CHD1wt_ETSneg<-apply(dat1[,ncol(dat_FPKM_CHD1del_v2):ncol(dat1)],1,mean)
dat1<-dat1[dat1$mean_CHD1del>=1 | dat1$mean_CHD1wt_ETSneg>=1,]
dat2<-cbind(dat1$V2, log2(dat1[,2:ncol(dat1)]+1)) #Log2 transformation
s1<-c(2:ncol(dat_FPKM_CHD1del)) #CHD1del
s2<-c((ncol(dat_FPKM_CHD1del)+1):(ncol(dat)-2)) #CHD1wt_ETSneg
wval<-NULL
pval<-NULL
qval<-NULL
for(x in 1:nrow(dat2))
{       wval<-c(wval,wilcox.test(as.matrix(dat2[x,s1]),as.matrix(dat2[x,s2]))$statistic) #use wilcoxon test to detect the significantly differentially expressed gene between two subtypes
pval<-c(pval,wilcox.test(as.matrix(dat2[x,s1]),as.matrix(dat2[x,s2]))$p.value)
}
qval<-p.adjust(pval, method ="BH", n = length(pval))
datall<-cbind(dat1,wval,pval,qval)
hist(qval, breaks=seq(0,1,by=0.01))

#Output all CHD1del sig in TCGA 333 samples
data_TCGA_333<-dat_FPKM_333_v2
data_TCGA_333_CHD1del<-merge(data_TCGA_333, data.frame(datall[,c(1,ncol(datall))]), by.x="V2", by.y="V2")
#write.csv(data_TCGA_333_CHD1del, file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR-333",".csv", sep=""), row.names=F)

#data_TCGA_333_CHD1del<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR-333.csv")
#FDR_cutoff<-9e-04
#datall_sig<-datall[datall$qval<FDR_cutoff,]
#nrow(datall_sig)
datall_sig <- read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR9e-04.csv")
nrow(datall_sig)
#Sig test in TCGA 333 samples
data_TCGA_333<-dat_FPKM_333_v2
data_TCGA_333_sig<-merge(data_TCGA_333, data.frame(datall_sig$V2), by.x="V2", by.y="datall_sig.V2")
#write.csv(data_TCGA_333_sig[,c(1,3:ncol(data_TCGA_333_sig))], file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR",FDR_cutoff,".csv", sep=""), row.names=F)
datinfo<-merge(data.frame(t(data_TCGA_333_sig[1,])), dinfo_TCGA, by.x="row.names", by.y="Sample.ID",sort=F)
datinfo$CHD1_exp<-as.numeric(data_TCGA_333[data_TCGA_333$V2=="CHD1",as.matrix(datinfo$Row.names)])
#CHD1_del<-ifelse(datinfo$CHD1_CNA=="homdel" & datinfo$ERG_status=="none" & datinfo$ETV1_status=="none"&datinfo$ETV4_status=="none"&datinfo$FLI1_status=="none", col_CHD1homo, 
#             ifelse(datinfo$CHD1_CNA=="hetloss" & datinfo$ERG_status=="none" & datinfo$ETV1_status=="none"&datinfo$ETV4_status=="none"&datinfo$FLI1_status=="none", col_CHD1hete, col_other))
CHD1_del<-ifelse(datinfo$CHD1_CNA=="homdel", col_CHD1homo, ifelse(datinfo$CHD1_CNA=="hetloss", col_CHD1hete, col_other))
#CHD1_exp<-colorRampPalette(brewer.pal(9,"Greys")[c(1,9)])(9)[as.numeric(cut(datinfo$CHD1_exp,breaks = 9))]
CNA_fraction<-colorRampPalette(brewer.pal(9,"OrRd")[c(1,9)])(9)[as.numeric(cut(datinfo$Fraction_genome_altered, breaks = 9))]
ERG_ETS<-ifelse(datinfo$ERG_status!="none", col_ERG, 
            ifelse(datinfo$ETV1_status!="none"|datinfo$ETV4_status!="none"|datinfo$FLI1_status!="none", col_ETS, col_other))
SPOP<-ifelse(datinfo$SPOP_mut==1, col_SPOP, col_other)
PTEN_del<-ifelse(datinfo$PTEN_CNA=="homdel", col_CHD1homo, ifelse(datinfo$PTEN_CNA=="hetloss", col_CHD1hete, col_other))
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR",FDR_cutoff, "-clusterW-333-hm.pdf", sep=""),width=10, height=10)
heatmap.3(as.matrix(log2(data_TCGA_333_sig[,3:ncol(data_TCGA_333_sig)]+1)), col=color2, dist=mydist, hclustfun=myclustw, 
          scale="row", key=TRUE,keysize=0.8,  symkey=FALSE, density.info="none", trace="none", 
          side.height.fraction=0.5, 
          #ColSideColors=cbind(CHD1_exp, CHD1_del, SPOP, ETS, CNA_fraction), 
          ColSideColors=cbind(CHD1_del, SPOP, PTEN_del, ERG_ETS), 
          labRow=data_TCGA_333_sig[,1],Colv=T, cexCol=0.2,cexRow=0.3,
          main=paste("TCGA-333-FPKM-wilcox_FDR", FDR_cutoff, nrow(data_TCGA_333_sig)))
legend("topright", legend=c("ERG_fusion", "ETS_fusion", "homdel", "hetloss", "SPOPmut"), col=c(col_ERG, col_ETS, col_CHD1homo, col_CHD1hete, col_SPOP), 
       pch=15, cex=0.5, bty="n")
dev.off()

#Heatmap annotation with only homdel
CHD1_del<-ifelse(datinfo$CHD1_CNA=="homdel", col_CHD1homo, col_other)
#CHD1_exp<-colorRampPalette(brewer.pal(9,"Greys")[c(1,9)])(9)[as.numeric(cut(datinfo$CHD1_exp,breaks = 9))]
CNA_fraction<-colorRampPalette(brewer.pal(9,"OrRd")[c(1,9)])(9)[as.numeric(cut(datinfo$Fraction_genome_altered, breaks = 9))]
ERG_ETS<-ifelse(datinfo$ERG_status!="none", col_ERG, 
                ifelse(datinfo$ETV1_status!="none"|datinfo$ETV4_status!="none"|datinfo$FLI1_status!="none", col_ETS, col_other))
SPOP<-ifelse(datinfo$SPOP_mut==1, col_SPOP, col_other)
PTEN_del<-ifelse(datinfo$PTEN_CNA=="homdel", col_CHD1homo, col_other)
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-CHD1del+CHD1wt_ETSneg-Wilcox-FDR",FDR_cutoff, "-clusterW-333-hm_homdel.pdf", sep=""),width=10, height=10)
heatmap.3(as.matrix(log2(data_TCGA_333_sig[,3:ncol(data_TCGA_333_sig)]+1)), col=color2, dist=mydist, hclustfun=myclustw, 
          scale="row", key=TRUE,keysize=0.8,  symkey=FALSE, density.info="none", trace="none", 
          side.height.fraction=0.5, 
          #ColSideColors=cbind(CHD1_exp, CHD1_del, SPOP, ETS, CNA_fraction), 
          ColSideColors=cbind(CHD1_del, SPOP, PTEN_del, ERG_ETS), 
          labRow=data_TCGA_333_sig[,1],Colv=T, cexCol=0.2,cexRow=0.3,
          main=paste("TCGA-333-FPKM-wilcox_FDR", FDR_cutoff, nrow(data_TCGA_333_sig)))
legend("topright", legend=c("ERG_fusion", "ETS_fusion", "homdel", "SPOPmut"), col=c(col_ERG, col_ETS, col_CHD1homo, col_SPOP), 
       pch=15, cex=0.5, bty="n")
dev.off()



##Differentially expressed genes between ETS_PTENdel and ETS_PTENwt
dat<-merge(dat_FPKM_ETS_PTENdel_v2, dat_FPKM_ETS_PTENwt_v2[,c(1,3:ncol(dat_FPKM_ETS_PTENwt_v2))], by="ENSG")
dat1<-dat[,c(2,3:ncol(dat))]
dat1$mean_ETS_PTENdel<-apply(dat1[,2:(ncol(dat_FPKM_ETS_PTENdel_v2)-1)],1,mean)
dat1$mean_ETS_PTENwt<-apply(dat1[,ncol(dat_FPKM_ETS_PTENdel_v2):ncol(dat1)],1,mean)
dat1<-dat1[dat1$mean_ETS_PTENdel>=1 | dat1$mean_ETS_PTENwt>=1,]
dat2<-cbind(dat1$V2, log2(dat1[,2:ncol(dat1)]+1)) #Log2 transformation
s1<-c(2:ncol(dat_FPKM_ETS_PTENdel)) 
s2<-c((ncol(dat_FPKM_ETS_PTENdel)+1):(ncol(dat)-2)) 
wval<-NULL
pval<-NULL
qval<-NULL
for(x in 1:nrow(dat2))
{       wval<-c(wval,wilcox.test(as.matrix(dat2[x,s1]),as.matrix(dat2[x,s2]))$statistic) #use wilcoxon test to detect the significantly differentially expressed gene between two subtypes
pval<-c(pval,wilcox.test(as.matrix(dat2[x,s1]),as.matrix(dat2[x,s2]))$p.value)
}
qval<-p.adjust(pval, method ="BH", n = length(pval))
datall<-cbind(dat1,wval,pval,qval)
hist(qval, breaks=seq(0,1,by=0.01))
#FDR_cutoff<-1e-05
#FDR_cutoff<-1e-06


#datall<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt_ETSpos-Wilcox.csv")
#FDR_cutoff<-6e-07
#datall_sig<-datall[datall$qval<FDR_cutoff,]
datall_sig <- read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt-Wilcox-FDR6e-07.csv", check.names = F)
nrow(datall_sig)

#Sig test in TCGA 333 samples
#data_TCGA_333<-merge(merge(merge(dat_FPKM_ETS_PTENdel_v2,  dat_FPKM_ETS_PTENwt_v2[,c(1,3:ncol(dat_FPKM_ETS_PTENwt_v2))], by="ENSG"), 
#                           dat_FPKM_CHD1del_v2[,c(1,3:ncol(dat_FPKM_CHD1del_v2))], by="ENSG"), dat_FPKM_other_v2[,c(1,3:ncol(dat_FPKM_other_v2))],  by="ENSG")
data_TCGA_333<-dat_FPKM_333_v2
data_TCGA_333_sig<-merge(data_TCGA_333, data.frame(datall_sig$V2), by.x="V2", by.y="datall_sig.V2")
#write.csv(data_TCGA_333_sig[,c(1,3:ncol(data_TCGA_333_sig))], file=paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt-Wilcox-FDR",FDR_cutoff,".csv", sep=""), row.names=F)

datinfo<-merge(t(data_TCGA_333_sig[1,]), dinfo_TCGA, by.x="row.names", by.y="Sample.ID",sort=F)
CHD1_del<-ifelse(datinfo$CHD1_CNA=="homdel", col_CHD1homo, ifelse(datinfo$CHD1_CNA=="hetloss", col_CHD1hete, col_other))
ERG_ETS<-ifelse(datinfo$ERG_status!="none", col_ERG, 
                ifelse(datinfo$ETV1_status!="none"|datinfo$ETV4_status!="none"|datinfo$FLI1_status!="none", col_ETS, col_other))
SPOP<-ifelse(datinfo$SPOP_mut==1, col_SPOP, col_other)
PTEN_del<-ifelse(datinfo$PTEN_CNA=="homdel", col_CHD1homo, ifelse(datinfo$PTEN_CNA=="hetloss", col_CHD1hete, col_other))
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt_ETSpos-Wilcox-FDR",FDR_cutoff, "-clusterW-333-hm.pdf", sep=""),width=10, height=10)
heatmap.3(as.matrix(log2(data_TCGA_333_sig[,3:ncol(data_TCGA_333_sig)]+1)), col=color2, dist=mydist, hclustfun=myclustw, 
          scale="row", key=TRUE,keysize=0.8,  symkey=FALSE, density.info="none", trace="none", 
          #side.height.fraction=0.4, ColSideColors=cbind(CHD1,SPOP,ETS), 
          #side.height.fraction=0.4, ColSideColors=cbind(CHD1,ETS,PTEN), 
          side.height.fraction=0.5, 
          ColSideColors=cbind(CHD1_del, SPOP, PTEN_del, ERG_ETS), 
          labRow=data_TCGA_333_sig[,1],Colv=T, cexCol=0.2,cexRow=0.75,
          main=paste("TCGA-333-FPKM-wilcox_FDR", FDR_cutoff, nrow(data_TCGA_333_sig)))
legend("topright", legend=c("ERG_fusion", "ETS_fusion", "homdel", "hetloss", "SPOPmut"), col=c(col_ERG, col_ETS, col_CHD1homo, col_CHD1hete, col_SPOP), 
       pch=15, cex=0.5, bty="n")
dev.off()

#Heatmap with homdel only
datinfo<-merge(t(data_TCGA_333_sig[1,]), dinfo_TCGA, by.x="row.names", by.y="Sample.ID",sort=F)
CHD1_del<-ifelse(datinfo$CHD1_CNA=="homdel", col_CHD1homo, col_other)
ERG_ETS<-ifelse(datinfo$ERG_status!="none", col_ERG, 
                ifelse(datinfo$ETV1_status!="none"|datinfo$ETV4_status!="none"|datinfo$FLI1_status!="none", col_ETS, col_other))
SPOP<-ifelse(datinfo$SPOP_mut==1, col_SPOP, col_other)
PTEN_del<-ifelse(datinfo$PTEN_CNA=="homdel", col_CHD1homo, col_other)
#pdf(paste("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-FPKM-PTENdel+PTENwt_ETSpos-Wilcox-FDR",FDR_cutoff, "-clusterW-333-hm_homdel.pdf", sep=""),width=10, height=10)
heatmap.3(as.matrix(log2(data_TCGA_333_sig[,3:ncol(data_TCGA_333_sig)]+1)), col=color2, dist=mydist, hclustfun=myclustw, 
          scale="row", key=TRUE,keysize=0.8,  symkey=FALSE, density.info="none", trace="none", 
          #side.height.fraction=0.4, ColSideColors=cbind(CHD1,SPOP,ETS), 
          #side.height.fraction=0.4, ColSideColors=cbind(CHD1,ETS,PTEN), 
          side.height.fraction=0.5, 
          ColSideColors=cbind(CHD1_del, SPOP, PTEN_del, ERG_ETS), 
          labRow=data_TCGA_333_sig[,1],Colv=T, cexCol=0.2,cexRow=0.75,
          main=paste("TCGA-333-FPKM-wilcox_FDR", FDR_cutoff, nrow(data_TCGA_333_sig)))
legend("topright", legend=c("ERG_fusion", "ETS_fusion", "homdel", "SPOPmut"), col=c(col_ERG, col_ETS, col_CHD1homo, col_SPOP), 
       pch=15, cex=0.5, bty="n")
dev.off()


