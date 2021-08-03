library(MASS)
library(gdata)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggpmisc)

#Figure 2E, S7, S8, S9: Upstream Regulator analysis for tumor samples

#Figure S7. Similar predicted upstream regulators from normal to early events between two tumor lineage models in TCGA cohort.
#IPA Upstream Regulator analysis of DEG from TCGA SPOPmut_CHD1wt/N and ERG_PTENwt/N
IPA_UR_TCGA_SPOPmut.CHD1wt_N <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=19)
IPA_UR_TCGA_ERG.PTENwt_N <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=20)
#IPA_UR_TCGA_SPOPmut.CHD1wt_N <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_SPOPmut_CHD1wt_N_rowSum_10_DESeq2_IPA_UR_output.xls", skip=1)
#IPA_UR_TCGA_ERG.PTENwt_N <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_ERG_PTENwt_N_rowSum_10_DESeq2_IPA_UR_output.xls", skip=1)
dat.merge <- merge(IPA_UR_TCGA_SPOPmut.CHD1wt_N[, c(1:5,6)], IPA_UR_TCGA_ERG.PTENwt_N[, c(1:5,7)], 
                   by="Upstream.Regulator")
#dat.merge <- merge(IPA_UR_TCGA_SPOPmut.CHD1wt_N[, c(1:6)], IPA_UR_TCGA_ERG.PTENwt_N[, c(1:5,7)], 
#                   by="Upstream.Regulator")
dat.merge.sig <- dat.merge[dat.merge$p.value.of.overlap.x<0.05 | dat.merge$p.value.of.overlap.y<0.05, ]
dat.merge.sig <- dat.merge.sig[complete.cases(dat.merge.sig$Activation.z.score.x) 
                               & complete.cases(dat.merge.sig$Activation.z.score.y), ]
dat.merge.sig <- dat.merge.sig[order(dat.merge.sig$p.value.of.overlap.x, dat.merge.sig$p.value.of.overlap.y, decreasing = F), ]
colnames(dat.merge.sig) <- c(colnames(dat.merge.sig)[1], 
                             paste(colnames(IPA_UR_TCGA_SPOPmut.CHD1wt_N[, c(2:6)]), ".SPOP.CHD1wt_N", sep=""), 
                             paste(colnames(IPA_UR_TCGA_ERG.PTENwt_N[, c(2:5,7)]), ".ERG.PTENwt_N", sep=""))
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP_CHD1wt_ERG_PTENwt_N_IPA_UR_output_cor_v4.pdf", width=14, height=10, useDingbats=FALSE)
ggplot(dat.merge.sig, aes(x= Activation.z.score.SPOP.CHD1wt_N, y=Activation.z.score.ERG.PTENwt_N)) +
  geom_point(aes(colour = factor(Molecule.Type.SPOP.CHD1wt_N)), size=4)  +
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  theme_classic() 
dev.off()


#Figure 2E: Divergent predicted upstream regulators by IPA software from “early” to “late” events between two tumor lineage models in TCGA cohort.
#IPA Upstream Regulator analysis of DEG from TCGA SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt
#IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_SPOPmut_CHD1wt_CHD1del_rowSum_10_DESeq2_IPA_UR_output.xls", skip=1)
#IPA_UR_TCGA_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA_T_ERG_PTENwt_PTENdel_rowSum_10_DESeq2_IPA_UR_output.xls", skip=1)
IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=21)
IPA_UR_TCGA_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=22)
#ggplot 
dat.merge <- merge(IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt[, c(1:5, 7)], IPA_UR_TCGA_ERG.PTENdel_PTENwt[, c(1:6)], 
                   by="Upstream.Regulator")
dat.merge.sig <- dat.merge[dat.merge$p.value.of.overlap.x<0.05 | dat.merge$p.value.of.overlap.y<0.05, ]
dat.merge.sig <- dat.merge.sig[complete.cases(dat.merge.sig$Activation.z.score.x) 
                               & complete.cases(dat.merge.sig$Activation.z.score.y), ]
dat.merge.sig <- dat.merge.sig[order(dat.merge.sig$p.value.of.overlap.x, dat.merge.sig$p.value.of.overlap.y, decreasing = F), ]
colnames(dat.merge.sig) <- c(colnames(dat.merge.sig)[1], 
                             paste(colnames(IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt[, c(2:5,7)]), ".SPOP.CHD1del_wt", sep=""), 
                             paste(colnames(IPA_UR_TCGA_ERG.PTENdel_PTENwt[, c(2:6)]), ".ERG.PTENdel_wt", sep=""))
dat.merge.sig_dif <- dat.merge.sig[dat.merge.sig$Activation.z.score.SPOP.CHD1del_wt/dat.merge.sig$Activation.z.score.ERG.PTENdel_wt<0 
                                   & (dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05 | dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05)
                                   & grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.SPOP.CHD1del_wt)
                                   & grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.ERG.PTENdel_wt), ]
dat.merge.sig_dif_label <- dat.merge.sig[dat.merge.sig$Activation.z.score.SPOP.CHD1del_wt/dat.merge.sig$Activation.z.score.ERG.PTENdel_wt<0 
                                         & (dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05 | dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05)
                                         & (grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.SPOP.CHD1del_wt)
                                            | grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.ERG.PTENdel_wt)), ]
#write.csv(dat.merge.sig_dif_label, 
#          file="/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP.CHD1del_CHD1wt_ERG.PTENdel_PTENwt_IPA_UR_output_diff.csv", row.names = F)

#UR label results from Mike
dat.merge.sig_label <- read.csv("/Users/deli/Desktop/Chris/Michael/TCGA/FromMike/Mike_EDIT_7_29_20.csv", check.names = F)
dat.merge.sig_label$Predicted.Activation <- ifelse(dat.merge.sig_label$Predicted.Activation.State.SPOP.CHD1del_wt=="Activated"
                                                   & dat.merge.sig_label$Predicted.Activation.State.ERG.PTENdel_wt=="Inhibited", "SPOP.Activated_ERG.Inhibited", 
                                                   ifelse(dat.merge.sig_label$Predicted.Activation.State.SPOP.CHD1del_wt=="Inhibited"
                                                          & dat.merge.sig_label$Predicted.Activation.State.ERG.PTENdel_wt=="Activated", "SPOP.Inhibited_ERG.Activated", 
                                                          ifelse(dat.merge.sig_label$Predicted.Activation.State.SPOP.CHD1del_wt=="Inhibited", "SPOP.Inhibited.only", 
                                                                 ifelse(dat.merge.sig_label$Predicted.Activation.State.ERG.PTENdel_wt=="Inhibited", "ERG.Inhibited.only",
                                                                        ifelse(dat.merge.sig_label$Predicted.Activation.State.SPOP.CHD1del_wt=="Activated", "SPOP.Activated.only", 
                                                                               ifelse(dat.merge.sig_label$Predicted.Activation.State.ERG.PTENdel_wt=="Activated", "ERG.Activated.only","NA"))))))
dat.merge.sig_label$family.v2 <- ifelse(dat.merge.sig_label$family.v1=="MEK" | dat.merge.sig_label$family.v1=="MEK_inhibitor"
                                        | dat.merge.sig_label$family.v1=="PI3K" | dat.merge.sig_label$family.v1=="PI3K_inhibitor"
                                        | dat.merge.sig_label$family.v1=="RTK" | dat.merge.sig_label$family.v1=="RTK_inhibitor"
                                        | dat.merge.sig_label$family.v1=="inflammation"
                                        | dat.merge.sig_label$family.v1=="TGFB"
                                        | dat.merge.sig_label$family.v1=="transcription_factor",
                                        as.matrix(dat.merge.sig_label$family.v1), "other")
dat.merge.sig_label$family.v2 <- factor(dat.merge.sig_label$family.v2, 
                                        level=c("MEK", "PI3K", "RTK", "inflammation", "TGFB", "transcription_factor", "cell_cycle", 
                                                "MEK_inhibitor", "PI3K_inhibitor", "RTK_inhibitor", "other"))

#UR label results from Mike and set others as grey
dat.merge <- merge(IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt[, c(1:5, 7)], IPA_UR_TCGA_ERG.PTENdel_PTENwt[, c(1:6)], 
                   by="Upstream.Regulator")
dat.merge.sig <- dat.merge[dat.merge$p.value.of.overlap.x<0.05 | dat.merge$p.value.of.overlap.y<0.05, ]
dat.merge.sig <- dat.merge.sig[complete.cases(dat.merge.sig$Activation.z.score.x) 
                               & complete.cases(dat.merge.sig$Activation.z.score.y), ]
dat.merge.sig <- dat.merge.sig[order(dat.merge.sig$p.value.of.overlap.x, dat.merge.sig$p.value.of.overlap.y, decreasing = F), ]
colnames(dat.merge.sig) <- c(colnames(dat.merge.sig)[1], 
                             paste(colnames(IPA_UR_TCGA_SPOPmut.CHD1del_CHD1wt[, c(2:5,7)]), ".SPOP.CHD1del_wt", sep=""), 
                             paste(colnames(IPA_UR_TCGA_ERG.PTENdel_PTENwt[, c(2:6)]), ".ERG.PTENdel_wt", sep=""))
dat.merge.sig_label_all <- merge(dat.merge.sig, dat.merge.sig_label[, c("Upstream.Regulator", "family", "family.v1", "family.v2")], by="Upstream.Regulator", all.x=T)
dat.merge.sig_label_all$family.v3 <- ifelse(complete.cases(dat.merge.sig_label_all$family.v2), as.matrix(dat.merge.sig_label_all$family.v2), "other")
dat.merge.sig_label_all$family.v3 <- factor(dat.merge.sig_label_all$family.v3, 
                                            level=c("MEK", "PI3K", "RTK", "inflammation", "TGFB", "transcription_factor", "cell_cycle", 
                                                    "MEK_inhibitor", "PI3K_inhibitor", "RTK_inhibitor", "other"))
dat.input <- dat.merge.sig_label_all
dat.input.name <- dat.input[dat.input$Upstream.Regulator=="AKT1"
                            | dat.input$Upstream.Regulator=="ERBB2"
                            | dat.input$Upstream.Regulator=="GATA2" 
                            #| dat.input$Upstream.Regulator=="IL6"
                            | dat.input$Upstream.Regulator=="LY294002"
                            | dat.input$Upstream.Regulator=="Map3k7"
                            | dat.input$Upstream.Regulator=="Mek"
                            | dat.input$Upstream.Regulator=="PD98059"
                            | dat.input$Upstream.Regulator=="PI3K (family)"
                            | dat.input$Upstream.Regulator=="TGFB2"
                            | dat.input$Upstream.Regulator=="TRIM24"
                            | dat.input$Upstream.Regulator=="U0126",]
#pdf("/Users/deli/Desktop/Chris/Michael/TCGA/TCGA-T_SPOP.CHD1del_CHD1wt_ERG.PTENdel_PTENwt_IPA_UR_output_sig_cor_v4.pdf", width=11.5, height=10, useDingbats=FALSE)
ggplot(dat.input, aes(x= Activation.z.score.SPOP.CHD1del_wt, y=Activation.z.score.ERG.PTENdel_wt)) +
  geom_point(aes(colour = factor(family.v3)), shape = 16, size=6)  +
  scale_color_manual(values=brewer.pal(12,"Set3")[c(1:8,10,9)]) +
  geom_label_repel(aes(x= Activation.z.score.SPOP.CHD1del_wt, y=Activation.z.score.ERG.PTENdel_wt, 
                       label=Upstream.Regulator, fill=factor(family.v3)), 
                   data = dat.input.name,
                   color = 'black', size = 9, label.size=0.1) + 
  scale_fill_manual(values=brewer.pal(12,"Set3")[c(1:8,10)]) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  theme_classic() 
dev.off()



#Figure S8: Different predicted upstream regulators from early to progression events between two tumor lineage models in Taylor cohort.
#IPA Upstream Regulator analysis of DEG from Taylor SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt
#IPA_UR_Taylor_pred_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_SPOP_CHD1del_CHD1wt__log2FC_0.4_IPA_UA.xls", skip=1)
#IPA_UR_Taylor_pred_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_GSE21034_pred_ERG_PTENdel_PTENwt__log2FC_0.2_IPA_UA.xls", skip=1)
IPA_UR_Taylor_pred_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=23)
IPA_UR_Taylor_pred_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=24)
#ggplot 
dat.merge <- merge(IPA_UR_Taylor_pred_SPOPmut.CHD1del_CHD1wt[, c(1:5,7)], IPA_UR_Taylor_pred_ERG.PTENdel_PTENwt[, c(1:5,7)], 
                   by="Upstream.Regulator")
dat.merge.sig <- dat.merge[dat.merge$p.value.of.overlap.x<0.05 | dat.merge$p.value.of.overlap.y<0.05, ]
dat.merge.sig <- dat.merge.sig[complete.cases(dat.merge.sig$Activation.z.score.x) 
                               & complete.cases(dat.merge.sig$Activation.z.score.y), ]
dat.merge.sig <- dat.merge.sig[order(dat.merge.sig$p.value.of.overlap.x, dat.merge.sig$p.value.of.overlap.y, decreasing = F), ]
colnames(dat.merge.sig) <- c(colnames(dat.merge.sig)[1], 
                             paste(colnames(IPA_UR_Taylor_pred_SPOPmut.CHD1del_CHD1wt[, c(2:5,7)]), ".SPOP.CHD1del_wt", sep=""), 
                             paste(colnames(IPA_UR_Taylor_pred_ERG.PTENdel_PTENwt[, c(2:5,7)]), ".ERG.PTENdel_wt", sep=""))
dat.merge.sig_dif <- dat.merge.sig[dat.merge.sig$Activation.z.score.SPOP.CHD1del_wt/dat.merge.sig$Activation.z.score.ERG.PTENdel_wt<0 
                                   & (dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05 | dat.merge.sig$p.value.of.overlap.ERG.PTENdel_wt<0.05)
                                   & (grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.SPOP.CHD1del_wt)
                                      | grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.ERG.PTENdel_wt)), ]
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/Taylor_pred-T_SPOP.CHD1del_CHD1wt_ERG.PTENdel_PTENwt_IPA_UR_output_cor_v3.pdf", width=14, height=10, useDingbats=FALSE)
ggplot(dat.merge.sig, aes(x= Activation.z.score.SPOP.CHD1del_wt, y=Activation.z.score.ERG.PTENdel_wt)) +
  geom_point(aes(colour = factor(Molecule.Type.SPOP.CHD1del_wt)), size=4)  +
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  theme_classic() 
dev.off()


#Figure S9. Different predicted upstream regulators from early to progression events between two tumor lineage models in ICGC cohort.
#IPA Upstream Regulator analysis of DEG from ICGC SPOP.CHD1del/CHD1wt and ERG.PTENdel/PTENwt
#IPA_UR_ICGC_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_SPOP_CHD1del_SPOP_CHD1wt_DESeq2_IPA_UR_output.xls", skip=1)
#IPA_UR_ICGC_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC_read_count_T_ERG_PTENdel_ERG_PTENwt_DESeq2_IPA_UR_output.xls", skip=1)
IPA_UR_ICGC_SPOPmut.CHD1del_CHD1wt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=25)
IPA_UR_ICGC_ERG.PTENdel_PTENwt <- read.xls("/Users/deli/Dropbox/Deli_LabMeeting/CHD1_clinic/Manuscript/SupplementaryTable.xls", skip=1, sheet=26)
#ggplot 
dat.merge <- merge(IPA_UR_ICGC_SPOPmut.CHD1del_CHD1wt[, c(1:6)], IPA_UR_ICGC_ERG.PTENdel_PTENwt[, c(1:6)], 
                   by="Upstream.Regulator")
dat.merge.sig <- dat.merge[dat.merge$p.value.of.overlap.x<0.05 | dat.merge$p.value.of.overlap.y<0.05, ]
dat.merge.sig <- dat.merge.sig[complete.cases(dat.merge.sig$Activation.z.score.x) 
                               & complete.cases(dat.merge.sig$Activation.z.score.y), ]
dat.merge.sig <- dat.merge.sig[order(dat.merge.sig$p.value.of.overlap.x, dat.merge.sig$p.value.of.overlap.y, decreasing = F), ]
colnames(dat.merge.sig) <- c(colnames(dat.merge.sig)[1], 
                             paste(colnames(IPA_UR_TCGA_SPOPmut.CHD1wt_N[, c(2:6)]), ".SPOP.CHD1del_wt", sep=""), 
                             paste(colnames(IPA_UR_TCGA_ERG.PTENwt_N[, c(2:5,7)]), ".ERG.PTENdel_wt", sep=""))
dat.merge.sig_dif <- dat.merge.sig[dat.merge.sig$Activation.z.score.SPOP.CHD1del_wt/dat.merge.sig$Activation.z.score.ERG.PTENdel_wt<0 
                                   & (dat.merge.sig$p.value.of.overlap.SPOP.CHD1del_wt<0.05 | dat.merge.sig$p.value.of.overlap.ERG.PTENdel_wt<0.05)
                                   & (grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.SPOP.CHD1del_wt)
                                      | grepl("Activated|Inhibited", dat.merge.sig$Predicted.Activation.State.ERG.PTENdel_wt)), ]
pdf("/Users/deli/Desktop/Chris/Michael/TCGA/ICGC-T_SPOP.CHD1del_CHD1wt_ERG.PTENdel_PTENwt_IPA_UR_output_cor_v3.pdf", width=14, height=10, useDingbats=FALSE)
ggplot(dat.merge.sig, aes(x= Activation.z.score.SPOP.CHD1del_wt, y=Activation.z.score.ERG.PTENdel_wt)) +
  geom_point(aes(colour = factor(Molecule.Type.SPOP.CHD1del_wt)), size=4)  +
  geom_smooth(method = "lm", se = FALSE, color="grey", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label =  paste(stat(rr.label), 
                                  stat(p.value.label), sep = "*\", \"*")),
               parse = TRUE) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  theme_classic() 
dev.off()


