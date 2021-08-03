library(ggplot2)
library(ggalluvial)
library(alluvial)
library(tidyr)
library(dplyr)

#Figure 5B: Alluvial diagrams of Gleason scores, lymph node invasion status and tumor stages from molecular subclasses in retrospective cohort. 
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$other=="Pos",])

svm.pred.info_clas_cohort_table<-data.frame(svm.pred.info_clas_cohort)
nrow(svm.pred.info_clas_cohort_table)
svm.pred.info_clas_cohort_table$Subtype<-ifelse(svm.pred.info_clas_cohort_table$ERG_PTENdel=="Pos", "ERG_PTENdel", 
                                                ifelse(svm.pred.info_clas_cohort_table$ERG_PTENwt=="Pos", "ERG_PTENwt", 
                                                       ifelse(svm.pred.info_clas_cohort_table$SPOP_signature=="Pos", "SPOP_signature",
                                                              ifelse(svm.pred.info_clas_cohort_table$CHD1_signature=="Pos", "CHD1_signature", 
                                                                     ifelse(svm.pred.info_clas_cohort_table$ETS_status=="Pos", "ETS","Other")))))
table(svm.pred.info_clas_cohort_table$Subtype)

#GS<7, 3+4, 4+3, GS>7
svm.pred.info_clas_cohort_table$category<-ifelse(as.numeric(svm.pred.info_clas_cohort_table$pathgs)<7, "GS<7",
                                                 ifelse(as.numeric(svm.pred.info_clas_cohort_table$pathgs_p)==3 & as.numeric(svm.pred.info_clas_cohort_table$pathgs_s)==4, "GS=3+4", 
                                                        ifelse(as.numeric(svm.pred.info_clas_cohort_table$pathgs_p)==4 & as.numeric(svm.pred.info_clas_cohort_table$pathgs_s)==3, "GS=4+3",
                                                        "GS>7")))
svm.pred.info_clas_cohort_table_summary<-data.frame(c(table(svm.pred.info_clas_cohort_table$Subtype, svm.pred.info_clas_cohort_table$category)[,1]/table(svm.pred.info_clas_cohort_table$category)[1], 
                                                      table(svm.pred.info_clas_cohort_table$Subtype, svm.pred.info_clas_cohort_table$category)[,2]/table(svm.pred.info_clas_cohort_table$category)[2], 
                                                      table(svm.pred.info_clas_cohort_table$Subtype, svm.pred.info_clas_cohort_table$category)[,3]/table(svm.pred.info_clas_cohort_table$category)[3],
                                                      table(svm.pred.info_clas_cohort_table$Subtype, svm.pred.info_clas_cohort_table$category)[,4]/table(svm.pred.info_clas_cohort_table$category)[4]))
colnames(svm.pred.info_clas_cohort_table_summary)<-"Frequency"
svm.pred.info_clas_cohort_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_cohort_table$Subtype, svm.pred.info_clas_cohort_table$category)), 4)
svm.pred.info_clas_cohort_table_summary$category<-c(rep("GS<7", 6), rep("GS=3+4", 6), rep("GS=4+3", 6), rep("GS>7", 6))
svm.pred.info_clas_cohort_table_summary$Subject<-rep(c(1:6), 4)
svm.pred.info_clas_cohort_table_summary$Subtype <- as.factor(svm.pred.info_clas_cohort_table_summary$Subtype)
svm.pred.info_clas_cohort_table_summary$Subtype <- factor(svm.pred.info_clas_cohort_table_summary$Subtype, 
                                                          levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-alluvium_GS_4.pdf")
ggplot(svm.pred.info_clas_cohort_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#LNI stage category
svm.pred.info_clas_cohort_table_lni<-svm.pred.info_clas_cohort_table
svm.pred.info_clas_cohort_table_lni$category<-ifelse(svm.pred.info_clas_cohort_table_lni$lni==1, "LNI=1", "LNI=0")
svm.pred.info_clas_cohort_table_summary<-data.frame(c(table(svm.pred.info_clas_cohort_table_lni$Subtype, svm.pred.info_clas_cohort_table_lni$category)[,1]/table(svm.pred.info_clas_cohort_table_lni$category)[1], 
                                                      table(svm.pred.info_clas_cohort_table_lni$Subtype, svm.pred.info_clas_cohort_table_lni$category)[,2]/table(svm.pred.info_clas_cohort_table_lni$category)[2]))
colnames(svm.pred.info_clas_cohort_table_summary)<-"Frequency"
svm.pred.info_clas_cohort_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_cohort_table_lni$Subtype, svm.pred.info_clas_cohort_table_lni$category)), 2)
svm.pred.info_clas_cohort_table_summary$category<-c(rep("LNI=0", 6), rep("LNI=1", 6))
svm.pred.info_clas_cohort_table_summary$Subject<-rep(c(1:6), 2)
svm.pred.info_clas_cohort_table_summary$Subtype <- factor(svm.pred.info_clas_cohort_table_summary$Subtype, 
                                                          levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-alluvium_LNI.pdf")
ggplot(svm.pred.info_clas_cohort_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#Tumor stage category
svm.pred.info_clas_cohort_table_pstage<-svm.pred.info_clas_cohort_table[grepl("T2|T3|T4", svm.pred.info_clas_cohort$pstage),]
svm.pred.info_clas_cohort_table_pstage$category<-ifelse(grepl("T2", svm.pred.info_clas_cohort_table_pstage$pstage), "T2", "T3|T4")
svm.pred.info_clas_cohort_table_summary<-data.frame(c(table(svm.pred.info_clas_cohort_table_pstage$Subtype, svm.pred.info_clas_cohort_table_pstage$category)[,1]/table(svm.pred.info_clas_cohort_table_pstage$category)[1], 
                                                      table(svm.pred.info_clas_cohort_table_pstage$Subtype, svm.pred.info_clas_cohort_table_pstage$category)[,2]/table(svm.pred.info_clas_cohort_table_pstage$category)[2]))
colnames(svm.pred.info_clas_cohort_table_summary)<-"Frequency"
svm.pred.info_clas_cohort_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_cohort_table_pstage$Subtype, svm.pred.info_clas_cohort_table_pstage$category)), 2)
svm.pred.info_clas_cohort_table_summary$category<-c(rep("T2", 6), rep("T3|T4", 6))
svm.pred.info_clas_cohort_table_summary$Subject<-rep(c(1:6), 2)
svm.pred.info_clas_cohort_table_summary$Subtype <- factor(svm.pred.info_clas_cohort_table_summary$Subtype, 
                                                          levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-SPOP-sig212-cost_0.01_CHD1-sig167-cost_svm-clas-Decipher-alluvium_TumorStage.pdf")
ggplot(svm.pred.info_clas_cohort_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()



#Figure 5D: Alluvial diagrams of Gleason scores, lymph node invasion status and tumor stages from molecular subclasses in Decipher prospective cohort.
#Aluvial plot to show the GS, LNI, Tumor stage difference across subclasses
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENdel=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ERG_PTENwt=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$ETS_status=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$SPOP_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_signature=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$CHD1_SPOP_share=="Pos",])
nrow(svm.pred.info_clas_6532[svm.pred.info_clas_6532$other=="Pos",])

svm.pred.info_clas_6532_table<-data.frame(svm.pred.info_clas_6532)
nrow(svm.pred.info_clas_6532_table)
svm.pred.info_clas_6532_table$Subtype<-ifelse(svm.pred.info_clas_6532_table$ERG_PTENdel=="Pos", "ERG_PTENdel", 
                                              ifelse(svm.pred.info_clas_6532_table$ERG_PTENwt=="Pos", "ERG_PTENwt", 
                                                     ifelse(svm.pred.info_clas_6532_table$SPOP_signature=="Pos", "SPOP_signature",
                                                            ifelse(svm.pred.info_clas_6532_table$CHD1_signature=="Pos", "CHD1_signature", 
                                                                   ifelse(svm.pred.info_clas_6532_table$ETS_status=="Pos", "ETS","Other")))))
table(svm.pred.info_clas_6532_table$Subtype)

#GS<7, 3+4, 4+3, GS>7
svm.pred.info_clas_6532_table$category<-ifelse(as.numeric(svm.pred.info_clas_6532_table$pathgs)<7, "GS<7",
                                               ifelse(as.numeric(svm.pred.info_clas_6532_table$pathgs_p)==3 & as.numeric(svm.pred.info_clas_6532_table$pathgs_s)==4, "GS=3+4", 
                                                      ifelse(as.numeric(svm.pred.info_clas_6532_table$pathgs_p)==4 & as.numeric(svm.pred.info_clas_6532_table$pathgs_s)==3, "GS=4+3",
                                                             "GS>7")))
svm.pred.info_clas_6532_table_summary<-data.frame(c(table(svm.pred.info_clas_6532_table$Subtype, svm.pred.info_clas_6532_table$category)[,1]/table(svm.pred.info_clas_6532_table$category)[1], 
                                                    table(svm.pred.info_clas_6532_table$Subtype, svm.pred.info_clas_6532_table$category)[,2]/table(svm.pred.info_clas_6532_table$category)[2], 
                                                    table(svm.pred.info_clas_6532_table$Subtype, svm.pred.info_clas_6532_table$category)[,3]/table(svm.pred.info_clas_6532_table$category)[3],
                                                    table(svm.pred.info_clas_6532_table$Subtype, svm.pred.info_clas_6532_table$category)[,4]/table(svm.pred.info_clas_6532_table$category)[4]))
colnames(svm.pred.info_clas_6532_table_summary)<-"Frequency"
svm.pred.info_clas_6532_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_6532_table$Subtype, svm.pred.info_clas_6532_table$category)), 4)
svm.pred.info_clas_6532_table_summary$category<-c(rep("GS<7", 6), rep("GS=3+4", 6), rep("GS=4+3", 6), rep("GS>7", 6))
svm.pred.info_clas_6532_table_summary$Subject<-rep(c(1:6), 4)
svm.pred.info_clas_6532_table_summary$Subtype <- as.factor(svm.pred.info_clas_6532_table_summary$Subtype)
svm.pred.info_clas_6532_table_summary$Subtype <- factor(svm.pred.info_clas_6532_table_summary$Subtype, 
                                                        levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-alluvium_GS_4.pdf")
ggplot(svm.pred.info_clas_6532_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#LNI stage category
svm.pred.info_clas_6532_table_lni<-svm.pred.info_clas_6532_table
svm.pred.info_clas_6532_table_lni$category<-ifelse(svm.pred.info_clas_6532_table_lni$lni==1, "LNI=1", "LNI=0")
svm.pred.info_clas_6532_table_summary<-data.frame(c(table(svm.pred.info_clas_6532_table_lni$Subtype, svm.pred.info_clas_6532_table_lni$category)[,1]/table(svm.pred.info_clas_6532_table_lni$category)[1], 
                                                    table(svm.pred.info_clas_6532_table_lni$Subtype, svm.pred.info_clas_6532_table_lni$category)[,2]/table(svm.pred.info_clas_6532_table_lni$category)[2]))
colnames(svm.pred.info_clas_6532_table_summary)<-"Frequency"
svm.pred.info_clas_6532_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_6532_table_lni$Subtype, svm.pred.info_clas_6532_table_lni$category)), 2)
svm.pred.info_clas_6532_table_summary$category<-c(rep("LNI=0", 6), rep("LNI=1", 6))
svm.pred.info_clas_6532_table_summary$Subject<-rep(c(1:6), 2)
svm.pred.info_clas_6532_table_summary$Subtype <- factor(svm.pred.info_clas_6532_table_summary$Subtype, 
                                                        levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-alluvium_LNI.pdf")
ggplot(svm.pred.info_clas_6532_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#Tumor stage category
svm.pred.info_clas_6532_table_pstage<-svm.pred.info_clas_6532_table[grepl("T2|T3|T4", svm.pred.info_clas_6532$pstage),]
svm.pred.info_clas_6532_table_pstage$category<-ifelse(grepl("T2", svm.pred.info_clas_6532_table_pstage$pstage), "T2", "T3|T4")
svm.pred.info_clas_6532_table_summary<-data.frame(c(table(svm.pred.info_clas_6532_table_pstage$Subtype, svm.pred.info_clas_6532_table_pstage$category)[,1]/table(svm.pred.info_clas_6532_table_pstage$category)[1], 
                                                    table(svm.pred.info_clas_6532_table_pstage$Subtype, svm.pred.info_clas_6532_table_pstage$category)[,2]/table(svm.pred.info_clas_6532_table_pstage$category)[2]))
colnames(svm.pred.info_clas_6532_table_summary)<-"Frequency"
svm.pred.info_clas_6532_table_summary$Subtype<-rep(rownames(table(svm.pred.info_clas_6532_table_pstage$Subtype, svm.pred.info_clas_6532_table_pstage$category)), 2)
svm.pred.info_clas_6532_table_summary$category<-c(rep("T2", 6), rep("T3|T4", 6))
svm.pred.info_clas_6532_table_summary$Subject<-rep(c(1:6), 2)
svm.pred.info_clas_6532_table_summary$Subtype <- factor(svm.pred.info_clas_6532_table_summary$Subtype, 
                                                        levels = c("SPOP_signature", "CHD1_signature", "ERG_PTENwt", "ERG_PTENdel", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/Decipher_prospective_6532_pred-alluvium_TumorStage.pdf")
ggplot(svm.pred.info_clas_6532_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_ERG_PTENwt, sur_col_ERG_PTENdel, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#Figure S20: Alluvial diagrams of Gleason scores, lymph node invasion status and tumor stages from molecular subclasses in TCGA cohort.
#TCGA sample information downloaded from cBioPortal
dinfo_TCGA_cBio<-read.table("/Users/deli/Desktop/Chris/TCGA-PCA/prad_tcga_clinical_data.tsv", header=T, sep="\t")
#Combine with TCGA paper information 
dinfo_TCGA_stage<-merge(dinfo_TCGA, dinfo_TCGA_cBio, by.x="SAMPLE_ID", by.y="Sample.ID")
table(dinfo_TCGA_stage$American.Joint.Committee.on.Cancer.Tumor.Stage.Code)
dinfo_TCGA_stage$other_status<-ifelse(dinfo_TCGA_stage$SPOP_mut==0 & dinfo_TCGA_stage$CHD1_CNA=="diploid" & dinfo_TCGA_stage$ERG_status=="none", "Pos", "Neg")
dinfo_TCGA_stage$Subtype<-ifelse(dinfo_TCGA_stage$SPOP_mut==1 & dinfo_TCGA_stage$CHD1_CNA=="diploid" 
                                 & (dinfo_TCGA_stage$ERG_status=="none" & dinfo_TCGA_stage$ETV1_status=="none"&dinfo_TCGA_stage$ETV4_status=="none"&dinfo_TCGA_stage$FLI1_status=="none"), "SPOPmut_only", 
                                 ifelse(dinfo_TCGA_stage$SPOP_mut==0 & dinfo_TCGA_stage$CHD1_CNA!="diploid" 
                                                               & (dinfo_TCGA_stage$ERG_status=="none" & dinfo_TCGA_stage$ETV1_status=="none"&dinfo_TCGA_stage$ETV4_status=="none"|dinfo_TCGA_stage$FLI1_status=="none"), "CHD1del_only",
                                 ifelse(dinfo_TCGA_stage$SPOP_mut==1 & dinfo_TCGA_stage$CHD1_CNA!="diploid" & (dinfo_TCGA_stage$ERG_status=="none" & dinfo_TCGA_stage$ETV1_status=="none"&dinfo_TCGA_stage$ETV4_status=="none"&dinfo_TCGA_stage$FLI1_status=="none"), "CHD1del+SPOPmut",
                                        ifelse(dinfo_TCGA_stage$ERG_status!="none" & dinfo_TCGA_stage$CHD1_CNA=="diploid" & dinfo_TCGA_stage$PTEN_CNA!="diploid", "ERG+_PTENdel",
                                        ifelse(dinfo_TCGA_stage$ERG_status!="none" & dinfo_TCGA_stage$CHD1_CNA=="diploid" & dinfo_TCGA_stage$PTEN_CNA=="diploid", "ERG+_PTENwt",
                                               ifelse((dinfo_TCGA_stage$ETV1_status!="none"|dinfo_TCGA_stage$ETV4_status!="none"|dinfo_TCGA_stage$FLI1_status!="none") & dinfo_TCGA_stage$CHD1_CNA=="diploid", "ETS", 
                                               "Other" ))))))
table(dinfo_TCGA_stage$Subtype)
dinfo_TCGA_stage_table<-data.frame(dinfo_TCGA_stage)

#Gleason score category: GS<7, GS=3+4, GS=4+3, GS>7
dinfo_TCGA_stage_table$category<-ifelse(dinfo_TCGA_stage_table$Reviewed_Gleason_category=="3+3", "GleasonScore=3+3",
                                        ifelse(dinfo_TCGA_stage_table$Reviewed_Gleason_category=="3+4", "GleasonScore=3+4",
                                               ifelse(dinfo_TCGA_stage_table$Reviewed_Gleason_category=="4+3", "GleasonScore=4+3", "GleasonScore>=8")))
dinfo_TCGA_stage_table_summary<-data.frame(c(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,1]/table(dinfo_TCGA_stage_table$category)[1], 
                                             table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,2]/table(dinfo_TCGA_stage_table$category)[2], 
                                             table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,3]/table(dinfo_TCGA_stage_table$category)[3],
                                             table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,4]/table(dinfo_TCGA_stage_table$category)[4]))
colnames(dinfo_TCGA_stage_table_summary)<-"Frequency"
dinfo_TCGA_stage_table_summary$Subtype<-rep(rownames(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)), 4)
dinfo_TCGA_stage_table_summary$category<-c(rep("GleasonScore=3+3", 7), rep("GleasonScore=3+4", 7), rep("GleasonScore=4+3", 7), rep("GleasonScore>=8", 7))
dinfo_TCGA_stage_table_summary$Subject<-rep(c(1:7), 4)
dinfo_TCGA_stage_table_summary$Subtype <- as.factor(dinfo_TCGA_stage_table_summary$Subtype)
dinfo_TCGA_stage_table_summary$Subtype <- factor(dinfo_TCGA_stage_table_summary$Subtype, 
                                                 levels = c("SPOPmut_only", "CHD1del_only", "CHD1del+SPOPmut", "ERG+_PTENdel", "ERG+_PTENwt", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-alluvium_GS_4.pdf")
ggplot(dinfo_TCGA_stage_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#Lymph node stage
table(dinfo_TCGA_cBio$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code)
dinfo_TCGA_stage_table$category<-ifelse(grepl("N1", dinfo_TCGA_stage_table$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code), "N1", "N0")
dinfo_TCGA_stage_table_summary<-data.frame(c(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,1]/table(dinfo_TCGA_stage_table$category)[1], 
                                             table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,2]/table(dinfo_TCGA_stage_table$category)[2]))
colnames(dinfo_TCGA_stage_table_summary)<-"Frequency"
dinfo_TCGA_stage_table_summary$Subtype<-rep(rownames(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)), 2)
dinfo_TCGA_stage_table_summary$category<-c(rep("N0", 7), rep("N1", 7))
dinfo_TCGA_stage_table_summary$Subject<-rep(c(1:7), 2)
dinfo_TCGA_stage_table_summary$Subtype <- as.factor(dinfo_TCGA_stage_table_summary$Subtype)
dinfo_TCGA_stage_table_summary$Subtype <- factor(dinfo_TCGA_stage_table_summary$Subtype, 
                                                 levels = c("SPOPmut_only", "CHD1del_only", "CHD1del+SPOPmut", "ERG+_PTENdel", "ERG+_PTENwt", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-alluvium_LymphNode.pdf")
ggplot(dinfo_TCGA_stage_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()


#Tumor stage category
dinfo_TCGA_stage_table$category<-ifelse(grepl("T2", dinfo_TCGA_stage_table$American.Joint.Committee.on.Cancer.Tumor.Stage.Code), "T2", "T3|T4")
dinfo_TCGA_stage_table_summary<-data.frame(c(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,1]/table(dinfo_TCGA_stage_table$category)[1], 
                                             table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)[,2]/table(dinfo_TCGA_stage_table$category)[2]))
colnames(dinfo_TCGA_stage_table_summary)<-"Frequency"
dinfo_TCGA_stage_table_summary$Subtype<-rep(rownames(table(dinfo_TCGA_stage_table$Subtype, dinfo_TCGA_stage_table$category)), 2)
dinfo_TCGA_stage_table_summary$category<-c(rep("T2", 7), rep("T3|T4", 7))
dinfo_TCGA_stage_table_summary$Subject<-rep(c(1:7), 2)
dinfo_TCGA_stage_table_summary$Subtype <- as.factor(dinfo_TCGA_stage_table_summary$Subtype)
dinfo_TCGA_stage_table_summary$Subtype <- factor(dinfo_TCGA_stage_table_summary$Subtype, 
                                                 levels = c("SPOPmut_only", "CHD1del_only", "CHD1del+SPOPmut", "ERG+_PTENdel", "ERG+_PTENwt", "ETS", "Other"))
#pdf("/Users/del2017/Desktop/Chris/Michael/TCGA/TCGA-freeze-333-SPOP-sig212_CHD1del-sig148_PTENdel-sig45/TCGA-freeze-333-alluvium_TumorStage.pdf")
ggplot(dinfo_TCGA_stage_table_summary,
       aes(x = category, stratum = Subtype, alluvium = Subject, y = Frequency,
           fill = Subtype, label = Subtype)) +
  #scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c(sur_col_SPOPo, sur_col_CHD1o, sur_col_CHD1_SPOP, sur_col_ERG_PTENdel, sur_col_ERG_PTENwt, sur_col_ETS, sur_col_other)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")
dev.off()




