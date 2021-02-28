setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure1")
library(data.table);library(tidyverse);library(IOBR);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")

exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
exp=exp[,pheno$samplesID]

GSEA_cluster1=c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                "REACTOME_INTERFERON_SIGNALING",
                "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES",
                "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
                "REACTOME_INNATE_IMMUNE_SYSTEM",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1")
GSEA_cluster2=c("HALLMARK_ADIPOGENESIS",
                "REACTOME_TRIGLYCERIDE_METABOLISM",
                "KEGG_INSULIN_SIGNALING_PATHWAY",
                "HALLMARK_FATTY_ACID_METABOLISM",
                "REACTOME_ASSEMBLY_OF_ACTIVE_LPL_AND_LIPC_LIPASE_COMPLEXES",
                "REACTOME_RA_BIOSYNTHESIS_PATHWAY",
                "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_WHITE_ADIPOCYTE_DIFFERENTIATION",
                "REACTOME_TRIGLYCERIDE_CATABOLISM")
GSEA_cluster3=c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY",
                "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                "REACTOME_CELL_CELL_COMMUNICATION",
                "HALLMARK_APICAL_JUNCTION",
                "REACTOME_CELL_JUNCTION_ORGANIZATION",
                "KEGG_FOCAL_ADHESION",
                "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
                "HALLMARK_COAGULATION",
                "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS")
GSEA_cluster4=c("HALLMARK_MTORC1_SIGNALING",
                "KEGG_CELL_CYCLE",
                "REACTOME_DNA_REPLICATION",
                "REACTOME_CELL_CYCLE_CHECKPOINTS",
                "REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE",
                "HALLMARK_MYC_TARGETS_V2",
                "REACTOME_G2_M_CHECKPOINTS",
                "REACTOME_APC_C_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS",
                "REACTOME_MITOTIC_G2_G2_M_PHASES",
                "REACTOME_DNA_REPLICATION_PRE_INITIATION",
                "REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY",
                "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS")
cluster=list(GSEA_cluster1,GSEA_cluster2,GSEA_cluster3,GSEA_cluster4)
P_table=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure1/LMEM_signature.rds")
row.names(P_table)=P_table$V1
P_table=as.data.frame(P_table[c(GSEA_cluster1,GSEA_cluster2,GSEA_cluster3,GSEA_cluster4),-1]) 

################
###Heatmap######
################
library(ComplexHeatmap);library(circlize);library(RColorBrewer)
#################################
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")
signature_score=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step5/signature_score.rds")
signature_score=scale(signature_score)
signature_score=as.data.frame(t(signature_score[pheno$samplesID,c(GSEA_cluster1,GSEA_cluster2,GSEA_cluster3,GSEA_cluster4)]))

mat1 = as.matrix(signature_score[,pheno$samplesID[pheno$tpt=="Baseline"]])   # baseline samples
baseline=pheno[pheno$samplesID%in%colnames(mat1),]
mat2 = as.matrix(signature_score[,pheno$samplesID[pheno$tpt=="Cycle 2"]])  # cycle2 samples
cycle2=pheno[pheno$samplesID%in%colnames(mat2),]
mat3 = as.matrix(signature_score[,pheno$samplesID[pheno$tpt=="Surgery"]]) # surgery samples
surgery=pheno[pheno$samplesID%in%colnames(mat3),]
cluster=c(rep("cluster1",times=8),rep("cluster2",times=8),rep("cluster3",times=10),rep("cluster4",times=12))

m1 = mat1
m2 = mat2
m3 = mat3

library(RColorBrewer)
colKI <- c(
  "gray" = rgb(128, 128, 128, maxColorValue=256), 
  "plum" = rgb(135, 0, 82, maxColorValue=256), 
  "aqua" = rgb(151, 216, 218, maxColorValue=256), 
  "lavender" = rgb(189, 171, 179, maxColorValue=256), 
  "teal" = rgb(136, 196, 197, maxColorValue=256), 
  "cyklamen" = rgb(212, 9, 99, maxColorValue=256))
colKI

pheno=baseline
ha1=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = 2),labels = c("T1 Baseline (n = 122)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black")),
                      na_col="grey",
                      border=TRUE)

pheno=cycle2
ha2=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = 3),labels = c("T2 Cycle2 (n = 82)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black")),
                      na_col="grey",
                      border=TRUE)


pheno=surgery
ha3=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = 4),labels = c("T3 Surgery (n = 71)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black")),
                      na_col="grey",
                      border=TRUE)

row_ha= rowAnnotation(T3vsT2= anno_barplot(-log10(as.numeric(P_table[,1]))),
                            T3vsT1= anno_barplot(-log10(as.numeric(P_table[,2]))),
                            T2vsT1= anno_barplot(-log10(as.numeric(P_table[,3]))))



col_fun = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
ht_list =Heatmap(cluster, name = "Cluster", show_row_names = FALSE, width = unit(10, "mm"),
          col = structure(names = c("cluster1", "cluster2", "cluster3", "cluster4"), c("#D30963","gold2","plum","#96D7D9")),
          row_split = factor(cluster, levels = c("cluster1","cluster2","cluster3","cluster4")))+ 
  Heatmap(m1, col = col_fun, name = "Scaled pathway/Signature score",
                  clustering_distance_columns = "spearman",
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  show_column_names = FALSE,
                  top_annotation = ha1,
                  row_split = factor(cluster, levels = c("cluster1","cluster2","cluster3","cluster4")), 
                  row_title_gp = gpar(col = "#FFFFFF00"), width = unit(10, "cm"))+
  Heatmap(m2,col=col_fun, show_column_names = FALSE, 
          show_column_dend = FALSE,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(8, "cm"))+
  Heatmap(m3,col=col_fun, show_column_names = FALSE, 
          show_column_dend = FALSE,top_annotation=ha3,
          right_annotation =row_ha,
          show_heatmap_legend = FALSE, width = unit(6, "cm"))

p0=draw(ht_list, row_title = paste0("PROMIX Trial"),
     annotation_legend_side = "left", heatmap_legend_side = "left") #8*20
####################line plot######################
ge=geom_errorbar(aes(ymin=LL,ymax=UL),width=0.2)#for absolute risks
gh=geom_hline(yintercept=1)
tc=function(sz) theme_classic(base_size=sz);
facet=facet_wrap(~cluster,nrow = 1)
gx=xlab("Time Series")
gy=ylab("Scaled GSVA score")
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure1")
library(data.table);library(tidyverse);library(IOBR);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(ggplot2);library(reshape2);library(dplyr)
signature_all=read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step5/GSEA_final.xlsx", col_names=TRUE, sheet=1)
signature_lum=read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step5/GSEA_final.xlsx", col_names=TRUE, sheet=2)
signature_TN=read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step5/GSEA_final.xlsx", col_names=TRUE, sheet=3)
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")
row.names(pheno)=pheno$samplesID
pheno$time[pheno$tpt=="Baseline"]="T1"
pheno$time[pheno$tpt=="Cycle 2"]="T2"
pheno$time[pheno$tpt=="Surgery"]="T3"
score=as.data.frame(scale(readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step5/signature_score.rds")))
score=score[pheno$samplesID,]
score$time=pheno[,"time"]
mydata=as.data.frame(t(group_by(score,time)%>%summarize_each(funs(mean))))
colnames(mydata)=mydata[1,]
mydata$name=row.names(mydata)
mydata=mydata[-1,]
data=mydata[signature_all$ID,]
data_all=as.data.frame(cbind(data,signature_all[,c("cluster")]))
data_all=gather(data_all,time,score,T1:T3,-cluster,-name)
data_all$score=as.numeric(data_all$score)
p1=ggplot(data_all, aes(x=time, y=score,group=name,col="grey"))+ geom_line(color ="grey")+ geom_point(color ="white")+ ylim(-0.8, 0.8)+facet_grid(.~cluster)+gx+gy+tc(14)+ theme(legend.position="none")
###6*12

cluster1=data_all[data_all$cluster=="cluster1",]
mean1=as.data.frame(tapply(as.numeric(cluster1$score),INDEX=as.factor(cluster1$time), FUN=mean))
mean1$time=c("T1","T2","T3")
sd=tapply(as.numeric(cluster1$score),INDEX=as.factor(cluster1$time), FUN=sd)
mean1=cbind(mean1,sd)
colnames(mean1)=c("score","time","sd")
mean1$cluster="cluster1"
mean1$name="cluster1"

cluster2=data_all[data_all$cluster=="cluster2",]
mean2=as.data.frame(tapply(as.numeric(cluster2$score),INDEX=as.factor(cluster2$time), FUN=mean))
mean2$time=c("T1","T2","T3")
sd=tapply(as.numeric(cluster2$score),INDEX=as.factor(cluster2$time), FUN=sd)
mean2=cbind(mean2,sd)
colnames(mean2)=c("score","time","sd")
mean2$cluster="cluster2"
mean2$name="cluster2"

cluster3=data_all[data_all$cluster=="cluster3",]
mean3=as.data.frame(tapply(as.numeric(cluster3$score),INDEX=as.factor(cluster3$time), FUN=mean))
mean3$time=c("T1","T2","T3")
sd=tapply(as.numeric(cluster3$score),INDEX=as.factor(cluster3$time), FUN=sd)
mean3=cbind(mean3,sd)
colnames(mean3)=c("score","time","sd")
mean3$cluster="cluster3"
mean3$name="cluster3"

cluster4=data_all[data_all$cluster=="cluster4",]
mean4=as.data.frame(tapply(as.numeric(cluster4$score),INDEX=as.factor(cluster4$time), FUN=mean))
mean4$time=c("T1","T2","T3")
sd=tapply(as.numeric(cluster4$score),INDEX=as.factor(cluster4$time), FUN=sd)
mean4=cbind(mean4,sd)
colnames(mean4)=c("score","time","sd")
mean4$cluster="cluster4"
mean4$name="cluster4"

final_data=rbind(mean1,mean2,mean3,mean4)
p2=ggplot(final_data, aes(x=time, y=score,group=name,col=cluster))+scale_color_manual(values=c("#D30963","gold2","plum","#96D7D9"))+ylim(-0.8, 0.8)+geom_line()+geom_errorbar(aes(ymin=score-sd, ymax=score+sd), width=.2)+geom_point()+facet_grid(.~cluster)+gx+gy+tc(14)+ theme(legend.position="none")

p0
#######################
#######output##########
#######################
library(ggpubr)
ggarrange(p1, p2 + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
10*12
