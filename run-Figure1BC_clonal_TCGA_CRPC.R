library(MASS)
library(gdata)

#Figure 1B: Clonal results
#dat_sample_summary <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/Clonal/clonality_data_TCGA_SU2C.xlsx", sheet=3)
dat_sample_summary <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", sheet=2, skip=1)
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/Clonal/clonality_data_TCGA_ggbar.pdf", width=6, height=6)
dat_sample_summary_TCGA<-dat_sample_summary[grepl("TCGA", dat_sample_summary$Study.ID), ]
dat_sample_summary_TCGA_alter<-data.frame(c(dat_sample_summary_TCGA$Fraction.clonal.aberrations, dat_sample_summary_TCGA$Fraction.subclonal.aberrations))
colnames(dat_sample_summary_TCGA_alter)<-"Alterations"
dat_sample_summary_TCGA_alter$gene<-factor(rep(dat_sample_summary_TCGA$Gene.ID, 2), levels = c("ERG", "PTEN", "SPOP", "CHD1"))
dat_sample_summary_TCGA_alter$Clonality_status<-factor(c(rep("clonal", 4), rep("subclonal", 4)), levels=c("subclonal", "clonal"))
ggplot(data=dat_sample_summary_TCGA_alter, aes(x=gene, y=Alterations, fill=Clonality_status)) +
  geom_bar(stat="identity") +
  scale_fill_manual("legend", values = c("clonal" = brewer.pal(9,"Set1")[9], "subclonal" = brewer.pal(9,"Set1")[2])) +
  coord_cartesian(ylim=c(0,0.2))
dev.off()


#Figure 1C: Enrichment of genomic alterations from localized prostate cancer to metastatic CRPC nominate progression events.
dinfo_TCGA_v2<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-sample_alter_input.csv", check.names = F)
dinfo_CRPC_v7<-read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/CRPC-sample_alter_input.csv", check.names = F)

#TCGA
tp<-dinfo_TCGA_v2[,c(1:53,56:62,66:78)]
tp<-tp[order(tp$SPOP_mut, tp$CHD1_CNA, tp$ERG_status, tp$ETV1_status, tp$ETV4_status, tp$PTEN_CNA, decreasing=T),]
#CRPC
tm<-dinfo_CRPC_v7[order(dinfo_CRPC_v7$SPOP_mut, dinfo_CRPC_v7$CHD1_CNA, dinfo_CRPC_v7$ERG_status, dinfo_CRPC_v7$ETS_status, dinfo_CRPC_v7$PTEN_CNA, decreasing=T),]


col_ERG<-brewer.pal(9,"Set1")[8]
#plot_ERG<-brewer.pal(9,"Greens")[8]
plot_ERG<-col_ERG
plot_ETS<-brewer.pal(9,"Greens")[5]
plot_back<-brewer.pal(12,"Set3")[9]

#Genomic alterations
col_mut<-brewer.pal(12,"Set3")[1]
col_cna_amp<-brewer.pal(12,"Set3")[4]
col_cna_homdel<-brewer.pal(12,"Set3")[5]
col_cna_hetloss<-brewer.pal(12,"Set3")[3]
#col_CHD1homo<-brewer.pal(9,"Set1")[2]
#col_CHD1hete<-brewer.pal(12,"Set3")[5]
plot_DEL_homo<-col_cna_homdel
plot_DEL_hete<-col_cna_hetloss
plot_AMP<-col_cna_amp
plot_MUT<-col_mut



#Plot the % and p-value between TCGA and CRPC altered genes: MUT+AMP+HomDEL
#AMP TCGA result from cBioPortal
tp_amp_3<-read.table("/Users/deli/Desktop/Chris/Michael/TCGA/cBioPortal_TCGA_oncoprint_AR_MYC_CCND1.tsv", sep="\t", header=T, check.names = F)
#AMP CRPC result from cBioPortal
tm_amp_3<-read.table("/Users/deli/Desktop/Chris/Michael/TCGA/cBioPortal_CRPC_oncoprint_AR_MYC_CCND1.tsv", sep="\t", header=T, check.names = F)
#MUT+ DEL + AMP (AR, MYC CCND1)
tp_genes<-merge(tp, tp_amp_3[,2:ncol(tp_amp_3)], by="SAMPLE_ID")
tp_genes$ETS_status<-ifelse(tp_genes$ETV1_status!="none" | tp_genes$ETV4_status!="none" | tp_genes$FLI1_status!="none", "fusion", "none")
tp_genes<-cbind(tp_genes$ERG_status, tp_genes$ETS_status, tp_genes[,seq((ncol(tp_genes)-42),(ncol(tp_genes)-1),by=1)] )
colnames(tp_genes)<-c("ERG_fusion", "ETS_fusion", "SPOP_mut", colnames(tp_genes[,4:ncol(tp_genes)]))
tp_genes_summary<-NULL
for (i in 1:ncol(tp_genes))
{
  if (grepl("fusion", colnames(tp_genes)[i]))
  {tp_genes_summary<-cbind(tp_genes_summary, as.numeric(as.character(length(tp_genes[tp_genes[,i]!="none",i])/length(tp_genes[,i]))))}
  else if (grepl("mut|MUT", colnames(tp_genes)[i]))
  {tp_genes_summary<-cbind(tp_genes_summary, as.numeric(as.character(length(tp_genes[tp_genes[,i]==1,i])/length(tp_genes[,i]))))}
  else if (grepl("CNA", colnames(tp_genes)[i]))
  {tp_genes_summary<-cbind(tp_genes_summary, as.numeric(as.character(length(tp_genes[tp_genes[,i]=="homdel",i])/length(tp_genes[,i]))))}
  else
    #{tp_genes_summary<-cbind(tp_genes_summary, as.numeric(as.character(length(tp_genes[tp_genes[,i]>0,i])/length(tp_genes[,i]))))}
  {tp_genes_summary<-cbind(tp_genes_summary, as.numeric(as.character(length(tp_genes[tp_genes[,i]==2,i])/length(tp_genes[,i]))))}
}

#tm gene summary
tm_genes<-merge(tm, tm_amp_3[,2:ncol(tp_amp_3)], by.x="PATIENT.ID", by="SAMPLE_ID")
tm_genes<-cbind(tm$ERG_status, tm$ETS_status, tm$SPOP_mut, tm_genes[,seq((ncol(tm_genes)-40),(ncol(tm_genes)),by=1)])
colnames(tm_genes)<-c("ERG_fusion", "ETS_fusion", "SPOP_mut", colnames(tm_genes[,4:ncol(tm_genes)]))
tm_genes_summary<-NULL
for (i in 1:ncol(tm_genes))
{
  if (grepl("fusion", colnames(tm_genes)[i]))
  {tm_genes_summary<-cbind(tm_genes_summary, as.numeric(as.character(length(tm_genes[tm_genes[,i]!="none",i])/length(tm_genes[,i]))))}
  else if (grepl("SPOP_mut", colnames(tm_genes)[i]))
  {tm_genes_summary<-cbind(tm_genes_summary, as.numeric(as.character(length(tm_genes[tm_genes[,i]==1,i])/length(tm_genes[,i]))))}
  else if (grepl("mut|MUT", colnames(tm_genes)[i]))
  {tm_genes_summary<-cbind(tm_genes_summary, sum(table(tm_genes[,i]))/length(tm_genes[,i]))}
  else if (grepl("CNA", colnames(tm_genes)[i]))
  {tm_genes_summary<-cbind(tm_genes_summary, as.numeric(as.character(length(tm_genes[grepl("-2",tm_genes[,i]),i])/length(tm_genes[,i]))))}
  else
    #{tm_genes_summary<-cbind(tm_genes_summary, as.numeric(as.character(length(tm_genes[tm_genes[,i]>0,i])/length(tm_genes[,i]))))}
  {tm_genes_summary<-cbind(tm_genes_summary, as.numeric(as.character(length(tm_genes[tm_genes[,i]==2,i])/length(tm_genes[,i]))))}
}

#plot(tm_genes_summary, tp_genes_summary): mut + del + amp
#pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/TCGA-figures//TCGA+CRPC-sample_alter_correlation-mut+homdel+amp.pdf", width=9, height=9)
pval_amp_AR<-fisher.test(rbind(c(nrow(tp_genes[tp_genes$AR==2,]), 333-nrow(tp_genes[tp_genes$AR==2,])), 
                               c(nrow(tm_genes[tm_genes$AR==2,]), 150-nrow(tm_genes[tm_genes$AR==2,]))))$p.value
pval_amp_MYC<-fisher.test(rbind(c(nrow(tp_genes[tp_genes$MYC==2,]), 333-nrow(tp_genes[tp_genes$MYC==2,])), 
                                c(nrow(tm_genes[tm_genes$MYC==2,]), 150-nrow(tm_genes[tm_genes$MYC==2,]))))$p.value
pval_amp_CCND1<-fisher.test(rbind(c(nrow(tp_genes[tp_genes$CCND1==2,]), 333-nrow(tp_genes[tp_genes$CCND1==2,])), 
                                  c(nrow(tm_genes[tm_genes$CCND1==2,]), 150-nrow(tm_genes[tm_genes$CCND1==2,]))))$p.value
pvals_amp<-data.frame(c(pval_amp_CCND1, pval_amp_MYC, pval_amp_AR))
colnames(pvals_amp)<-"p.value"
pvals_amp$info<-as.factor(c("CCND1_AMP", "MYC_AMP", "AR_AMP"))
pvals_amp$log.pval<-0-log10(as.numeric(as.character(pvals_amp$p.value)))
pvals_amp$NAME<-c("CCND1_AMP*", "MYC_AMP", "AR_AMP*")

pvals_combine<-rbind(pvals_amp, pvals)
summary<-data.frame(pvals_combine[pvals_combine$NAME>0,][nrow(pvals_combine[pvals_combine$NAME>0,]):1,])
summary$TCGA_genes<-t(tp_genes_summary)
summary$CRPC_genes<-t(tm_genes_summary)
cex_val <- log2(summary$log.pval+2)
plot.col <- ifelse(grepl("ERG", summary$NAME), plot_ERG, 
                   ifelse(grepl("ETS", summary$NAME), plot_ETS, 
                          ifelse(grepl("SPOP|mut|MUT", summary$NAME), plot_MUT, 
                                 ifelse(grepl("AMP", summary$NAME), plot_AMP, plot_DEL_homo))))
plot.pch <- ifelse(as.numeric(as.character(summary$p.value))<0.05, 16, 1)
plot(summary$CRPC_genes, summary$TCGA_genes, 
     xlim=c(0,0.6), ylim=c(0,0.6), 
     pch=plot.pch, cex=cex_val, col=plot.col)
gene_label<-summary[summary$NAME=="ERG" | summary$NAME=="ETS" | summary$NAME=="SPOP" ,]
text(gene_label$CRPC_genes, gene_label$TCGA_genes, labels=gene_label$NAME, pos=1, cex=0.75, col="red")
gene_sig<-summary[as.numeric(as.character(summary$p.value))<0.01,]
text(gene_sig$CRPC_genes, gene_sig$TCGA_genes, labels=gene_sig$NAME, pos=1, cex=0.5, col="red")
abline(0,1)
legend("topleft", legend = c("MUT", "DEL_hom", "AMP"), 
       col=c(plot_MUT, plot_DEL_homo, plot_AMP), pch=16, cex=0.75)
#dev.off()

#plot(tm_genes_summary, tp_genes_summary): mut + del + amp, and annotate Figure 1A genes (pval<0.05)
pdf("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/TCGA-figures/TCGA+CRPC-sample_alter_correlation-mut+homdel+amp_select.pdf", width=9, height=9, useDingbats=FALSE)
nrow(summary)
nrow(pvals_TCGA_sort)
cex_val <- log2(summary$log.pval+2)
plot.col <- ifelse(grepl("ERG", summary$NAME), plot_ERG, 
                   ifelse(grepl("ETS", summary$NAME), plot_ETS, 
                          ifelse(grepl("SPOP|mut|MUT", summary$NAME), plot_MUT, 
                                 ifelse(grepl("AMP", summary$NAME), plot_AMP, plot_DEL_homo))))
plot.pch <- ifelse(as.numeric(as.character(summary$p.value))<0.05, 16, 1)
plot(summary$CRPC_genes, summary$TCGA_genes, 
     xlim=c(0,0.6), ylim=c(0,0.6), 
     pch=plot.pch, cex=cex_val, col=plot.col)
#gene_label<-summary[summary$NAME=="ERG" | summary$NAME=="ETS" | summary$NAME=="SPOP" ,]
pvals_TCGA_sort$info_gene<-gsub("_alter", "", pvals_TCGA_sort$info)
#pvals_TCGA_sort$info_gene<-gsub("_CNA|_mut", "", pvals_TCGA_sort$info)
summary$info_gene<-gsub("_CNA|_mut", "", summary$info)
#gene_label<-merge(summary, pvals_TCGA_sort[pvals_TCGA_sort$p.value<0.05,], by="info_gene")
gene_label <- summary[grepl("ERG|SPOP|PTEN_CNA|CHD1|FOXA1|AR|SPOPL|KMT2A", summary$info), ]
#text(gene_label$CRPC_genes, gene_label$TCGA_genes, labels=gene_label$info.x, pos=4, cex=0.75, col="red")
text(gene_label$CRPC_genes, gene_label$TCGA_genes, labels=gene_label$info_gene, pos=4, cex=0.75, col="red")
abline(0,1)
legend("topleft", legend = c("MUT", "DEL_hom+hete", "AMP"), 
       col=c(plot_MUT, plot_DEL_homo, plot_AMP), pch=16, cex=0.75)
dev.off()





