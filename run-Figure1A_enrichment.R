library(MASS)
library(gplots)
library(gdata)
#library(heatmap.plus)
library(beeswarm)
library("VennDiagram")
require("heatmap.plus")
library(RColorBrewer)

#Figure S1. Different enrichment of recurrent genomic alterations from ERG fusion and SPOP mutant subclasses in TCGA localized prostate cancer cohort.

#subtype color
plot.col=c(brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[2], brewer.pal(9,"Set1")[4], 
           brewer.pal(9,"Set1")[1],brewer.pal(12,"Set3")[6],brewer.pal(9,"Set1")[3])
#col_ERG<-brewer.pal(9,"Set1")[1]
col_ERG<-brewer.pal(9,"Set1")[8]
#col_ETS<-brewer.pal(9,"Set1")[8]
col_ETS<-sur_col_ETS
col_PTEN<-brewer.pal(12,"Set3")[8]
col_CHD1<-brewer.pal(9,"Set1")[2]
col_other<-brewer.pal(9,"Set1")[3]
col_SPOPo<-brewer.pal(9,"Set1")[5]
col_CHD1o<-brewer.pal(9,"Set1")[2]
col_CHD1_SPOP<-brewer.pal(9,"Set1")[4]
col_CHD1homo<-brewer.pal(9,"Set1")[2]
col_CHD1hete<-brewer.pal(12,"Set3")[5]

#Genomic alterations
col_mut<-brewer.pal(12,"Set3")[1]
col_cna_amp<-brewer.pal(12,"Set3")[4]
col_cna_homdel<-brewer.pal(12,"Set3")[5]
col_cna_hetloss<-brewer.pal(12,"Set3")[3]

#Input TCGA PCA paper TableS2 with sample cluster information
dinfo_cluster<-read.xls("/Users/deli/Desktop/Chris/Papers/TCGA-Prostate/mmc2.xls", sheet=1)
dinfo_clusterm<-dinfo_cluster
dinfo_clusterm<-dinfo_clusterm[order(dinfo_clusterm$ERG_status, dinfo_clusterm$ETV1_status, dinfo_clusterm$ETV4_status,
                                     dinfo_clusterm$SPOP_mut, dinfo_clusterm$CHD1_CNA, dinfo_clusterm$PTEN_CNA,decreasing=F),]
#dinfo_clusterm<-dinfo_clusterm[order(dinfo_clusterm$mRNA_cluster, dinfo_clusterm$methylation_cluster, dinfo_clusterm$SCNA ), ]
dinfo_clusterm_ERG<-dinfo_clusterm[dinfo_clusterm$'Subtype'=="1-ERG",]
dinfo_clusterm_ETS<-dinfo_clusterm[dinfo_clusterm$'Subtype'=="2-ETV1" | dinfo_clusterm$'Subtype'=="3-ETV4" | dinfo_clusterm$'Subtype'=="4-FLI1",]
dinfo_clusterm_SPOP<-dinfo_clusterm[dinfo_clusterm$'Subtype'=="5-SPOP",]
dinfo_clusterm_other<-dinfo_clusterm[dinfo_clusterm$'Subtype'!="1-ERG" 
                                     & dinfo_clusterm$'Subtype'!="2-ETV1" & dinfo_clusterm$'Subtype'!="3-ETV4" & dinfo_clusterm$'Subtype'!="4-FLI1" 
                                     & dinfo_clusterm$'Subtype'!="5-SPOP",]

#Output gene
pvals_TCGA_gene<-unique(gsub("_[A-z][A-z][A-z]", "", pvals_TCGA_sort$NAME))

#Gene: mut + cna 
AKT1_mut<-dinfo_clusterm[, grepl("AKT1_mut|AKT1_MUT", colnames(dinfo_clusterm))]
AKT1_cna<-rep("diploid", nrow(dinfo_clusterm))
ATM_mut<-dinfo_clusterm[, grepl("ATM_mut|ATM_MUT", colnames(dinfo_clusterm))]
ATM_cna<-rep("diploid", nrow(dinfo_clusterm))
BRAF_mut<-dinfo_clusterm[, grepl("BRAF_mut|BRAF_MUT", colnames(dinfo_clusterm))]
BRAF_cna<-rep("diploid", nrow(dinfo_clusterm))
BRCA1_mut<-dinfo_clusterm[, grepl("BRCA1_mut|BRCA1_MUT", colnames(dinfo_clusterm))]
BRCA1_cna<-dinfo_clusterm[, grepl("BRCA1_CNA", colnames(dinfo_clusterm))]
BRCA2_mut<-dinfo_clusterm[, grepl("BRCA2_mut|BRCA2_MUT", colnames(dinfo_clusterm))]
BRCA2_cna<-dinfo_clusterm[, grepl("BRCA2_CNA", colnames(dinfo_clusterm))]
CDK12_mut<-dinfo_clusterm[, grepl("CDK12_mut|CDK12_MUT", colnames(dinfo_clusterm))]
CDK12_cna<-dinfo_clusterm[, grepl("CDK12_CNA", colnames(dinfo_clusterm))]
CDKN1B_mut<-dinfo_clusterm[, grepl("CDKN1B_mut|CDKN1B_MUT", colnames(dinfo_clusterm))]
CDKN1B_cna<-dinfo_clusterm[, grepl("CDKN1B_CNA", colnames(dinfo_clusterm))]
CHD1_mut<-dinfo_clusterm[, grepl("CHD1_mut|CHD1_MUT", colnames(dinfo_clusterm))]
CHD1_cna<-dinfo_clusterm[, grepl("CHD1_CNA", colnames(dinfo_clusterm))]
CTNNB1_mut<-dinfo_clusterm[, grepl("CTNNB1_mut|CTNNB1_MUT", colnames(dinfo_clusterm))]
CTNNB1_cna<-rep("diploid", nrow(dinfo_clusterm))
FAM175A_mut<-rep(0, nrow(dinfo_clusterm))
FAM175A_cna<-dinfo_clusterm[, grepl("FAM175A_CNA", colnames(dinfo_clusterm))]
FANCC_mut<-dinfo_clusterm[, grepl("FANCC_mut|FANCC_MUT", colnames(dinfo_clusterm))]
FANCC_cna<-dinfo_clusterm[, grepl("FANCC_CNA", colnames(dinfo_clusterm))]
FANCD2_mut<-dinfo_clusterm[, grepl("FANCD2_mut|FANCD2_MUT", colnames(dinfo_clusterm))]
FANCD2_cna<-dinfo_clusterm[, grepl("FANCD2_CNA", colnames(dinfo_clusterm))]
FOXA1_mut<-dinfo_clusterm[, grepl("FOXA1_mut|FOXA1_MUT", colnames(dinfo_clusterm))]
FOXA1_cna<-rep("diploid", nrow(dinfo_clusterm))
HRAS_mut<-dinfo_clusterm[, grepl("HRAS_mut|HRAS_MUT", colnames(dinfo_clusterm))]
HRAS_cna<-rep("diploid", nrow(dinfo_clusterm))
IDH1_mut<-dinfo_clusterm[, grepl("IDH1_mut|IDH1_MUT", colnames(dinfo_clusterm))]
IDH1_cna<-rep("diploid", nrow(dinfo_clusterm))
KDM6A_mut<-dinfo_clusterm[, grepl("KDM6A_mut|KDM6A_MUT", colnames(dinfo_clusterm))]
KDM6A_cna<-rep("diploid", nrow(dinfo_clusterm))
KMT2A_mut<-dinfo_clusterm[, grepl("KMT2A_mut|KMT2A_MUT", colnames(dinfo_clusterm))]
KMT2A_cna<-rep("diploid", nrow(dinfo_clusterm))
KMT2C_mut<-dinfo_clusterm[, grepl("KMT2C_mut|KMT2C_MUT", colnames(dinfo_clusterm))]
KMT2C_cna<-rep("diploid", nrow(dinfo_clusterm))
MED12_mut<-dinfo_clusterm[, grepl("MED12_mut|MED12_MUT", colnames(dinfo_clusterm))]
MED12_cna<-rep("diploid", nrow(dinfo_clusterm))
PIK3CA_mut<-dinfo_clusterm[, grepl("PIK3CA_mut|PIK3CA_MUT", colnames(dinfo_clusterm))]
PIK3CA_cna<-rep("diploid", nrow(dinfo_clusterm))
PTEN_mut<-dinfo_clusterm[, grepl("PTEN_mut|PTEN_MUT", colnames(dinfo_clusterm))]
PTEN_cna<-dinfo_clusterm[, grepl("PTEN_CNA", colnames(dinfo_clusterm))]
RAD51C_mut<-rep(0, nrow(dinfo_clusterm))
RAD51C_cna<-dinfo_clusterm[, grepl("RAD51C_CNA", colnames(dinfo_clusterm))]
RB1_mut<-dinfo_clusterm[, grepl("RB1_mut|RB1_MUT", colnames(dinfo_clusterm))]
RB1_cna<-dinfo_clusterm[, grepl("RB1_CNA", colnames(dinfo_clusterm))]
SETD2_mut<-dinfo_clusterm[, grepl("SETD2_mut|SETD2_MUT", colnames(dinfo_clusterm))]
SETD2_cna<-rep("diploid", nrow(dinfo_clusterm))
SPOPL_mut<-rep(0, nrow(dinfo_clusterm))
SPOPL_cna<-dinfo_clusterm[, grepl("SPOPL_CNA", colnames(dinfo_clusterm))]
TP53_mut<-dinfo_clusterm[, grepl("TP53_mut|TP53_MUT", colnames(dinfo_clusterm))]
TP53_cna<-dinfo_clusterm[, grepl("TP53_CNA", colnames(dinfo_clusterm))]
ZMYM3_mut<-dinfo_clusterm[, grepl("ZMYM3_mut|ZMYM3_MUT", colnames(dinfo_clusterm))]
ZMYM3_cna<-rep("diploid", nrow(dinfo_clusterm))

dinfo_clusterm_alter<-cbind(dinfo_clusterm[,1:35], AKT1_mut, AKT1_cna, ATM_mut, ATM_cna, BRAF_mut, BRAF_cna, BRCA1_mut, BRCA1_cna, BRCA2_mut, BRCA2_cna, CDK12_mut, CDK12_cna, CDKN1B_mut, CDKN1B_cna,
                            CHD1_mut, CHD1_cna, CTNNB1_mut, CTNNB1_cna, FAM175A_mut, FAM175A_cna, FANCC_mut, FANCC_cna, FANCD2_mut, FANCD2_cna, FOXA1_mut, FOXA1_cna, HRAS_mut, HRAS_cna, IDH1_mut, IDH1_cna,
                            KDM6A_mut, KDM6A_cna, KMT2A_mut, KMT2A_cna, KMT2C_mut, KMT2C_cna, MED12_mut, MED12_cna, PIK3CA_mut, PIK3CA_cna, PTEN_mut, PTEN_cna, RAD51C_mut, RAD51C_cna, RB1_mut, RB1_cna,
                            SETD2_mut, SETD2_cna, SPOPL_mut, SPOPL_cna, TP53_mut, TP53_cna, ZMYM3_mut, ZMYM3_cna)
#write.csv(dinfo_clusterm_alter, file="/Users/deli/Desktop/Chris/TCGA-PCA/PRAD_TCGA_annotation_333_alter_gene.csv", row.names = F)

#ERG and SPOP subclasses
dinfo_clusterm_alter_ERG<-dinfo_clusterm_alter[dinfo_clusterm_alter$'Subtype'=="1-ERG",]
dinfo_clusterm_alter_ETS<-dinfo_clusterm_alter[dinfo_clusterm_alter$'Subtype'=="2-ETV1" | dinfo_clusterm_alter$'Subtype'=="3-ETV4" | dinfo_clusterm_alter$'Subtype'=="4-FLI1",]
dinfo_clusterm_alter_SPOP<-dinfo_clusterm_alter[dinfo_clusterm_alter$'Subtype'=="5-SPOP",]
dinfo_clusterm_alter_other<-dinfo_clusterm_alter[dinfo_clusterm_alter$'Subtype'!="1-ERG" 
                                     & dinfo_clusterm_alter$'Subtype'!="2-ETV1" & dinfo_clusterm_alter$'Subtype'!="3-ETV4" & dinfo_clusterm_alter$'Subtype'!="4-FLI1" 
                                     & dinfo_clusterm_alter$'Subtype'!="5-SPOP",]

#p-value for each alterations from TCGA in new barplot
n1<-nrow(dinfo_clusterm_alter_ERG)
n2<-nrow(dinfo_clusterm_alter_SPOP)
pval_alter<-NULL
data_alter<-NULL
for (j3 in c(seq(36,ncol(dinfo_clusterm_alter),by=2)))
{
  data_alter<-rbind(data_alter,
                    cbind(nrow(dinfo_clusterm_alter_ERG[dinfo_clusterm_alter_ERG[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_ERG[,j3+1]),])/n1, 
                          nrow(dinfo_clusterm_alter_SPOP[dinfo_clusterm_alter_SPOP[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_SPOP[,j3+1]),])/n2))
  pval_alter<-c(pval_alter, 
                fisher.test(rbind(c(nrow(dinfo_clusterm_alter_ERG[dinfo_clusterm_alter_ERG[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_ERG[,j3+1]),]), 
                                  n1-nrow(dinfo_clusterm_alter_ERG[dinfo_clusterm_alter_ERG[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_ERG[,j3+1]),])), 
                                c(nrow(dinfo_clusterm_alter_SPOP[dinfo_clusterm_alter_SPOP[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_SPOP[,j3+1]),]), 
                                  n2-nrow(dinfo_clusterm_alter_SPOP[dinfo_clusterm_alter_SPOP[,j3]==1 | grepl("homdel", dinfo_clusterm_alter_SPOP[,j3+1]),]))))$p.value)
}
pvals<-data.frame(cbind(data_alter, pval_alter))
colnames(pvals)<-c("ERG_group_ratio", "SPOP_group_ratio", "p.value")
pvals$info<-as.factor(paste(gsub("_mut", "", colnames(dinfo_clusterm_alter[, seq(36,ncol(dinfo_clusterm_alter),by=2)])), "_alter", sep=""))
pvals$log.pval<-0-log10(as.numeric(as.character(pvals$p.value)))
pvals<-pvals[seq(dim(pvals)[1],1),]
#pvals$NAME<-ifelse(as.numeric(as.character(pvals$p.value))<0.05, paste(pvals$info, "*",sep=""), as.matrix(pvals$info))
pvals$diff<-pvals$ERG_group_ratio-pvals$SPOP_group_ratio
#Sort by the pvals
pvals_sort<-pvals[order(pvals$p.value, decreasing = T),]
pvals_TCGA_sort<-pvals_sort

#Enrichment p-values in Figure S1
pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Figures/TCGA-freeze-333-sample-cluster-ERG_SPOP_alter_mut_amp_homdel_pval_sort.pdf", width=6, height=9, useDingbats=FALSE)
p<-ggplot(data=pvals_sort, aes(x=info, y=log.pval, fill=diff<0)) +
  geom_bar(stat="identity") + 
  scale_fill_manual("Subclass", values = c(sur_col_ERG_PTENwt, sur_col_SPOPo)) +
  coord_flip() + 
  scale_x_discrete(limits=as.factor(pvals_sort$info)) + 
  geom_hline(yintercept = 1.3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"))
p
dev.off()

#Enrichment p-values in Figure 1A
#pval<0.05 genes
pvals_sort_select<-pvals_sort[pvals_sort$p.value<0.05,]
pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Figures/TCGA-freeze-333-sample-cluster-ERG_SPOP_alter_mut_amp_homdel_pval_0.05_sort.pdf", width=6, height=4, useDingbats=FALSE)
p<-ggplot(data=pvals_sort_select, aes(x=info, y=log.pval, fill=diff<0)) +
  geom_bar(stat="identity") + 
  scale_fill_manual("Subclass", values = c(sur_col_ERG_PTENwt, sur_col_SPOPo)) +
  coord_flip() + 
  scale_x_discrete(limits=as.factor(pvals_sort_select$info)) + 
  geom_hline(yintercept = 1.3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "white"))
p
dev.off()



#Figure S1 and 1A (p-values<0.05 genes only)
pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Figures/TCGA-freeze-333-sample-cluster_mut_amp_homdel.pdf",width=12, height=8, useDingbats=FALSE)

par(mar=c(5.1,5,5,3), mgp=c(10,0.5,0))
#Add background
plot(1:(nrow(dinfo_clusterm_alter)+6), rep(1,(nrow(dinfo_clusterm_alter)+6)), 
     xlim=c(0,nrow(dinfo_clusterm_alter)+4), ylim=c(0,9),col="white", xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
     main="")

#sample info
h_gene_info=9
h_gene=8.5

pvals_sort_gene<-rev(gsub("_alter", "", pvals_sort$info))

#plot sample info ERG
dinfom<-dinfo_clusterm_alter_ERG[order(dinfo_clusterm_alter_ERG$PTEN_cna, decreasing=T),]
for (i in 1: nrow(dinfom))
{
  xstart<-0
  input<-dinfom[i,]
  ERG=input$'ERG_status'
  
  col.Age <- colorRampPalette(brewer.pal(9,"Greys")[c(1,9)])(9)[as.numeric(cut(dinfom$Age, breaks = 10))]
  
  #color 
  plot_ERG<-ifelse(grepl("none", input$'ERG_status'), "grey",  col_ERG)
  
  #sample info.
  rect(i+xstart-1, h_gene_info-0.25, i+xstart, h_gene_info, border = NA, col=plot_ERG)
  
  #CNA plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_cna<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,2]
    height<-h_gene-j1*0.25
    plot.col_CNA<-ifelse(grepl("gain", input_cna), col_cna_amp,
                         ifelse(grepl("homdel", input_cna), col_cna_homdel, brewer.pal(9,"Set3")[9]))
    #ifelse(grepl("hetloss", input_cna), col_cna_hetloss, "white")))
    rect(i+xstart-1,height-0.25, i+xstart,height, border = NA, col=plot.col_CNA)
  }
  
  #MUT plot
  for (j2 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_mut<-input[,grepl(pvals_sort_gene[j2], colnames(input))][,1]
    height<-h_gene-j2*0.25
    if(input_mut>0){rect(i+xstart-1,height-0.25, i+xstart,height-0.125, border = NA, col=col_mut)}
  }
  
}

#plot sample info ETS
dinfom<-dinfo_clusterm_alter_ETS[order(dinfo_clusterm_alter_ETS$PTEN_cna, decreasing=T),]
for (i in 1: nrow(dinfom))
{
  xstart<-nrow(dinfo_clusterm_ERG)+2
  input<-dinfom[i,]
  ETS=input$'ERG_status'
  
  #color 
  plot_ETS<-ifelse(grepl("none", input$'ETV1_status')&grepl("none", input$'ETV4_status')&grepl("none", input$'FLI1_status'), 
                   "grey", col_ETS)
  
  #sample info.
  rect(i+xstart-1, h_gene_info-0.25, i+xstart, h_gene_info, border = NA, col=plot_ETS)
  
  #CNA plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_cna<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,2]
    height<-h_gene-j1*0.25
    plot.col_CNA<-ifelse(grepl("gain", input_cna), col_cna_amp,
                         ifelse(grepl("homdel", input_cna), col_cna_homdel, brewer.pal(9,"Set3")[9]))
    #ifelse(grepl("hetloss", input_cna), col_cna_hetloss, "white")))
    rect(i+xstart-1,height-0.25, i+xstart,height, border = NA, col=plot.col_CNA)
  }
  
  #MUT plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_mut<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,1]
    height<-h_gene-j1*0.25
    if(input_mut>0){rect(i+xstart-1,height-0.25, i+xstart,height-0.125, border = NA, col=col_mut)}
  }
  
}


#plot sample info SPOP
dinfom<-dinfo_clusterm_alter_SPOP[order(dinfo_clusterm_alter_SPOP$CHD1_cna, decreasing=T),]
for (i in 1: nrow(dinfom))
{
  xstart<-nrow(dinfo_clusterm_ERG)+nrow(dinfo_clusterm_ETS)+4
  input<-dinfom[i,]
  
  #color 
  plot_SPOP<-ifelse(dinfom[i,]$SPOP_mut==1, col_SPOP, "grey")
  
  #sample info.
  rect(i+xstart-1, h_gene_info-0.25, i+xstart, h_gene_info, border = NA, col=plot_SPOP)
  
  #CNA plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_cna<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,2]
    height<-h_gene-j1*0.25
    plot.col_CNA<-ifelse(grepl("gain", input_cna), col_cna_amp,
                         ifelse(grepl("homdel", input_cna), col_cna_homdel, brewer.pal(9,"Set3")[9]))
    rect(i+xstart-1,height-0.25, i+xstart,height, border = NA, col=plot.col_CNA)
  }
  
  #MUT plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_mut<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,1]
    height<-h_gene-j1*0.25
    if(input_mut>0){rect(i+xstart-1,height-0.25, i+xstart,height-0.125, border = NA, col=col_mut)}
  }
  
}

#plot sample info other
dinfom<-dinfo_clusterm_alter_other[order(dinfo_clusterm_alter_other$CHD1_cna, decreasing=T),]
for (i in 1: nrow(dinfom))
{
  xstart<-nrow(dinfo_clusterm_ERG)+nrow(dinfo_clusterm_ETS)+nrow(dinfo_clusterm_SPOP)+6
  input<-dinfom[i,]

  #color 
  plot_other<-col_other
  
  #sample info.
  rect(i+xstart-1, h_gene_info-0.25, i+xstart, h_gene_info, border = NA, col=plot_other)
  
  #CNA plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_cna<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,2]
    height<-h_gene-j1*0.25
    plot.col_CNA<-ifelse(grepl("gain", input_cna), col_cna_amp,
                         ifelse(grepl("homdel", input_cna), col_cna_homdel, brewer.pal(9,"Set3")[9]))
    rect(i+xstart-1,height-0.25, i+xstart,height, border = NA, col=plot.col_CNA)
  }
  
  #MUT plot
  for (j1 in seq(1,length(pvals_sort_gene),by=1))
  {
    input_mut<-input[,grepl(pvals_sort_gene[j1], colnames(input))][,1]
    height<-h_gene-j1*0.25
    if(input_mut>0){rect(i+xstart-1,height-0.25, i+xstart,height-0.125, border = NA, col=col_mut)}
  }
  
}

#y axis label
axis(2, at=h_gene_info-0.125, label="Subclass", las=2, col="white", cex.axis=0.75)

#x axis label
axis(3, at=nrow(dinfo_clusterm_ERG)/2, label="ERG", las=0, col="white", cex.axis=0.75)
axis(3, at=nrow(dinfo_clusterm_ERG)+2+nrow(dinfo_clusterm_ETS)/2, label="ETS", las=0, col="white", cex.axis=0.75)
axis(3, at=nrow(dinfo_clusterm_ERG)+2+nrow(dinfo_clusterm_ETS)+2+nrow(dinfo_clusterm_SPOP)/2, label="SPOP", las=0, col="white", cex.axis=0.75)
axis(3, at=nrow(dinfo_clusterm_ERG)+2+nrow(dinfo_clusterm_ETS)+2+nrow(dinfo_clusterm_SPOP)+2+nrow(dinfo_clusterm_other)/2, label="Other", las=0, col="white", cex.axis=0.75)

#Alter gene legend
for (j1 in seq(1,length(pvals_sort_gene),by=1))
{
  height<-h_gene-j1*0.25+0.375
  axis(2, at=height-0.5, 
       label=pvals_sort_gene[j1],
       tcl=0, las=2, col="white", cex.axis=0.5)
}

#Add legend
h_last=2
h=1.25
legend(0, h_last-h, title="Genomic alteration", legend=c("Amp", "Homdel", "Mutation"), 
       col=c(col_cna_amp, col_cna_homdel, col_mut), 
       pch=15, cex=0.5, bty="n", horiz=TRUE,border=F,
       x.intersp=0.5, xjust=0, yjust=0, text.width=15)
dev.off()


