library(tidyr);library(tibble);library(dplyr);library(survival);library(tableone)    # data manipulation
library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
library(iClusterPlus);library(ConsensusClusterPlus);library(NbClust);library(factoextra);library(tidyverse)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7")
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
mydata=readRDS("pheno_immune.rds")
row.names(mydata)=mydata$samplesID
ADAPTIVE=c("B_cell_PCA_16704732","Bcell_receptors_score",
           "T cells","CD8 T cells","Treg cells",
           "Tcell_receptors_score","T helper cells",
           "Tcm cells","Tgd cells","Th1 cells","Th17 cells",
           "Th2 cells","MHC.I_19272155","MHC.II_19272155","IgG_19272155")
INNATE=c("Macrophages","TAMsurr_score","TAMsurr_TcClassII_ratio","DAP12_data",
         "aDC","DC","pDC","iDC","CD68","CD8_CD68_ratio","CSF1_response","Eosinophils",
         "Mast cells","Neutrophils","NK cells","NK CD56bright cells","NK CD56dim cells")
Stroma=c("CHANG_CORE_SERUM_RESPONSE_UP","CSR_Activated_15701700","Troester_WoundSig_19887484",
         "Lymph vessels","Angiogenesis")
Cytokine=c("IFN_21978456","IL2_score_21050467","IL4_score_21050467","IL8_21978456","IL12_score_21050467","IL13_score_21050467",
           "IFNG_score_21050467","TGFB_score_21050467")
TAAs=c("EMT1","EMT2","EMT3","MUC1","HER2_Immune_PCA_18006808")
Checkpoint=c("Immune_Checkpoint","PD1_data","PDL1_data","PD1_PDL1_score","CTLA4_data","ICR_INHIB_SCORE")
cell_fraction=c("CD4_memory_T_cells","CD8_T_cells","Naive_B_cells","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells","Follicular_Helper_T_cells","Mast_cells","NK_cells")

ID_T2=pheno$samplesID[pheno$tpt=="Cycle 2"]
TME1=c(ADAPTIVE,INNATE,Stroma,Cytokine,TAAs,Checkpoint)
TME1=mydata[ID_T2,TME1]
TME2=as.matrix(mydata[ID_T2,cell_fraction])
TME=cbind(TME1,TME2)
TME<-apply(TME,2,function(x){(x-mean(x))/sd(x)})
results<-ConsensusClusterPlus(t(TME),maxK=5,reps = 1000, 
                              pItem = 0.8, pFeature = 1, title="TME_cluster",
                              clusterAlg='km', distance="euclidean", seed=123456,
                              plot="pdf", #或"png"
                              corUse="pairwise.complete.obs",writeTable=T)
clusters =as.data.frame(results[[3]][["consensusClass"]])
colnames(clusters)="TME"
clusters$TME_cluster_T2[clusters$TME=="1"]="3"
clusters$TME_cluster_T2[clusters$TME=="2"]="1"
clusters$TME_cluster_T2[clusters$TME=="3"]="2"
clusters$samplesID=ID_T2
table(clusters$TME_cluster_T2)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
data=left_join(data,clusters,by="samplesID")
saveRDS(data,file="pheno_TME_T2.rds")


library(data.table);library(tidyverse);library(IOBR);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_T2.rds")
table(mydata$TME_cluster_T2)
mydata=mydata[!is.na(mydata$TME_cluster_T2),]
row.names(mydata)=mydata$samplesID
mydata=mydata[order(mydata$TME_cluster_T2),]
ADAPTIVE=c("B_cell_PCA_16704732","Bcell_receptors_score",
           "T cells","CD8 T cells","Treg cells",
           "Tcell_receptors_score","T helper cells",
           "Tcm cells","Tgd cells","Th1 cells","Th17 cells",
           "Th2 cells","MHC.I_19272155","MHC.II_19272155","IgG_19272155")
INNATE=c("Macrophages","TAMsurr_score","TAMsurr_TcClassII_ratio","DAP12_data",
         "aDC","DC","pDC","iDC","CD68","CD8_CD68_ratio","CSF1_response","Eosinophils",
         "Mast cells","Neutrophils","NK cells","NK CD56bright cells","NK CD56dim cells")
Stroma=c("CHANG_CORE_SERUM_RESPONSE_UP","CSR_Activated_15701700","Troester_WoundSig_19887484",
         "Lymph vessels","Angiogenesis")
Cytokine=c("IFN_21978456","IL2_score_21050467","IL4_score_21050467","IL8_21978456","IL12_score_21050467","IL13_score_21050467",
           "IFNG_score_21050467","TGFB_score_21050467")
TAAs=c("EMT1","EMT2","EMT3","MUC1","HER2_Immune_PCA_18006808")
Checkpoint=c("Immune_Checkpoint","PD1_data","PDL1_data","PD1_PDL1_score","CTLA4_data","ICR_INHIB_SCORE")

TME=c(ADAPTIVE,INNATE,Stroma,Cytokine,TAAs,Checkpoint)
TME_score=scale(mydata[,TME])
TME_score=as.data.frame(t(TME_score[mydata$samplesID,TME]))

mat1 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster_T2=="1"]])   # baseline samples
TME_cluster1=mydata[mydata$samplesID%in%colnames(mat1),]
mat2 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster_T2=="2"]])   # baseline samples
TME_cluster2=mydata[mydata$samplesID%in%colnames(mat2),]
mat3 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster_T2=="3"]])   # baseline samples
TME_cluster3=mydata[mydata$samplesID%in%colnames(mat3),]

TME_function=c(rep("ADAPTIVE",times=15),
               rep("INNATE",times=17),
               rep("Stroma",times=5),
               rep("Cytokine",times=8),
               rep("TAAs",times=5),
               rep("Checkpoint",times=6))

library(RColorBrewer)
colKI <- c(
  "gray" = rgb(128, 128, 128, maxColorValue=256), 
  "plum" = rgb(135, 0, 82, maxColorValue=256), 
  "aqua" = rgb(151, 216, 218, maxColorValue=256), 
  "lavender" = rgb(189, 171, 179, maxColorValue=256), 
  "teal" = rgb(136, 196, 197, maxColorValue=256), 
  "cyklamen" = rgb(212, 9, 99, maxColorValue=256))

pheno=TME_cluster1
max(pheno$StromalScore)
ha1=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#6699FF"),labels = c("Cold (n = 100)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      PAM50=pheno$PAM50,
                      subtype_immunity=pheno$subtype_immunity,
                      Time=pheno$tpt,
                      #StromalScore=pheno$StromalScore,
                      #ImmuneScore=pheno$ImmuneScore,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black"),
                               PAM50=structure(names = c("LumA", "LumB", "Her2", "Basal", 
                                                         "Normal"), brewer.pal(5,"Set1")),
                               subtype_immunity=c("IFN-γ dominant"="#E41A1C","Inflammatory"="#377EB8",
                                                  "TGF-β dominant"="#4DAF4A","Lymphocyte depleted"="#984EA3",
                                                  "Wound healing"="#FF7F00"),
                               Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      #StromalScore=colorRamp2(c(3,5,8),c("blue","white","red")),
                      #ImmuneScore=colorRamp2(c(3,5,8),c("blue","white","red"))),
                      na_col="grey",
                      border=TRUE)

pheno=TME_cluster2
ha2=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#FFCC00"),labels = c("Warm (n = 118)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      PAM50=pheno$PAM50,
                      subtype_immunity=pheno$subtype_immunity,
                      Time=pheno$tpt,
                      #StromalScore=pheno$StromalScore,
                      #ImmuneScore=pheno$ImmuneScore,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black"),
                               PAM50=structure(names = c("LumA", "LumB", "Her2", "Basal", 
                                                         "Normal"), brewer.pal(5,"Set1")),
                               subtype_immunity=c("IFN-γ dominant"="#E41A1C","Inflammatory"="#377EB8",
                                                  "TGF-β dominant"="#4DAF4A","Lymphocyte depleted"="#984EA3",
                                                  "Wound healing"="#FF7F00"),
                               Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      #StromalScore=colorRamp2(c(3,5,8),c("blue","white","red")),
                      #ImmuneScore=colorRamp2(c(3,5,8),c("blue","white","red"))),
                      na_col="grey",
                      border=TRUE)


pheno=TME_cluster3
ha3=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#CC0000"),labels = c("Hot (n = 57)")),
                      Tumor_purity= anno_barplot(pheno$purity),
                      pCR=pheno$`pCR or not`,
                      death=pheno$status,
                      PAM50=pheno$PAM50,
                      subtype_immunity=pheno$subtype_immunity,
                      Time=pheno$tpt,
                      #StromalScore=pheno$StromalScore,
                      #ImmuneScore=pheno$ImmuneScore,
                      col=list(pCR=c("non-pCR"="white","pCR"="black"),
                               death=c("0"="white","1"="black"),
                               PAM50=structure(names = c("LumA", "LumB", "Her2", "Basal", 
                                                         "Normal"), brewer.pal(5,"Set1")),
                               subtype_immunity=c("IFN-γ dominant"="#E41A1C","Inflammatory"="#377EB8",
                                                  "TGF-β dominant"="#4DAF4A","Lymphocyte depleted"="#984EA3",
                                                  "Wound healing"="#FF7F00"),
                               Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      #StromalScore=colorRamp2(c(3,5,8),c("blue","white","red")),
                      #ImmuneScore=colorRamp2(c(3,5,8),c("blue","white","red"))),
                      na_col="grey",
                      border=TRUE)
table(mydata$TME_cluster)

col_fun = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
ht_list =Heatmap(TME_function, name = "signature", show_row_names = FALSE, width = unit(10, "mm"),
                 col = structure(names = c("ADAPTIVE","INNATE","Stroma","Cytokine","TAAs","Checkpoint"),
                                 brewer.pal(6,"Set2")),
                 row_split = factor(TME_function, levels = c("ADAPTIVE","INNATE","Stroma","Cytokine","TAAs","Checkpoint")))+ 
  Heatmap(mat1, col = col_fun, name = "Scaled pathway/Signature score",
          clustering_distance_columns = "spearman",
          show_row_dend = FALSE, show_column_dend = FALSE,
          show_column_names = FALSE,
          top_annotation = ha1,
          row_split = factor(TME_function, levels = c("ADAPTIVE","INNATE","Stroma","Cytokine","TAAs","Checkpoint")), 
          row_title_gp = gpar(col = "#FFFFFF00"), width = unit(7, "cm"))+
  Heatmap(mat2,col=col_fun, show_column_names = FALSE, 
          show_column_dend = FALSE,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(9, "cm"))+
  Heatmap(mat3,col=col_fun, show_column_names = FALSE, 
          show_column_dend = FALSE,top_annotation=ha3,
          show_heatmap_legend = FALSE, width = unit(4, "cm"))

p0=draw(ht_list, row_title = paste0("EMT subtype"),
        annotation_legend_side = "left", heatmap_legend_side = "left") 








data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_T2.rds")
table(data$TME_cluster_T2,data$TME_cluster)
data=data[order(data$TME_cluster),]
data$pCR[data$pCR=="pCR"]=1
data$pCR[data$pCR=="non-pCR"]=0
data$LN=9
data$LN[data$`Regional nodes BL.y`=="Yes"]=1
data$LN[data$`Regional nodes BL.y`=="No"]=0
data$TS=9
data$TS[data$Tumorsize<50]=1
data$TS[data$Tumorsize>=50]=2
variable=c("TME_cluster","patientsID","DFS.status","DFS.time","IHC_subtype","pCR","LN","TS")
Baseline=data[data$tpt=="Baseline",c(variable)]
colnames(Baseline)=c("TME1","patientsID","DFS.status","DFS.time","IHC_subtype","pCR","LN","TS")
Cycle2=data[data$tpt=="Cycle 2",c(variable,"TME_cluster_T2")]
colnames(Cycle2)=c("TME2","patientsID","DFS.status2","DFS.time2","IHC_subtype2","pCR2","LN2","TS2","TME_cluster_T2")
Surgery=data[data$tpt=="Surgery",variable]
colnames(Surgery)=c("TME3","patientsID","DFS.status3","DFS.time3","IHC_subtype3","pCR3","LN3","TS3")
dat1=full_join(Baseline,Cycle2,by="patientsID")
df=full_join(dat1,Surgery,by="patientsID")

df$DFS.status[is.na(df$DFS.status)]=df$DFS.status2[is.na(df$DFS.status)]
df$DFS.status[is.na(df$DFS.status)]=df$DFS.status3[is.na(df$DFS.status)]
df$DFS.time[is.na(df$DFS.time)]=df$DFS.time2[is.na(df$DFS.time)]
df$DFS.time[is.na(df$DFS.time)]=df$DFS.time3[is.na(df$DFS.time)]
df$IHC_subtype[is.na(df$IHC_subtype)]=df$IHC_subtype2[is.na(df$IHC_subtype)]
df$IHC_subtype[is.na(df$IHC_subtype)]=df$IHC_subtype3[is.na(df$IHC_subtype)]
df$pCR[is.na(df$pCR)]=df$pCR2[is.na(df$pCR)]
df$pCR[is.na(df$pCR)]=df$pCR3[is.na(df$pCR)]
df$LN[is.na(df$LN)]=df$LN2[is.na(df$LN)]
df$LN[is.na(df$LN)]=df$LN3[is.na(df$LN)]
df$TS[is.na(df$TS)]=df$TS2[is.na(df$TS)]
df$TS[is.na(df$TS)]=df$TS3[is.na(df$TS)]
###warm,hot->cold was defined as that tumor(baseline,cycle2) immune class changed from hot or warm to cold
###cold->warm,hot change was defined as that tumor(baseline,cycle2) immune class changed from hot or warm to cold
####Define cold#####
df=df[!is.na(df$TME1),]
df=df[!is.na(df$TME_cluster_T2),]
df$cold=0
df$cold[df$TME1%in%c("1")]=1
table(df$cold)
df$group_cold=0
df$group_cold[df$TME1%in%c("1")&df$TME_cluster_T2%in%c("2","3")]=1
cold=df[df$cold=="1",]
table(cold$group_cold,cold$pCR)
table(cold$group_cold,cold$DFS.status)
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_cold), data=cold) ##final model
ShowRegTable(cox.test)
N_logistic <- glm(as.numeric(pCR) ~ as.factor(group_cold)+as.factor(LN)+as.factor(TS)+as.factor(IHC_subtype), data= cold) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)
df$IHC_subtype

df$hot=0
df$hot[df$TME1%in%c("2","3")]=1
table(df$hot)
df$group_hot=0
df$group_hot[df$TME1%in%c("3","2")&df$TME_cluster_T2%in%c("1")]=1
df$group_hot[df$TME1%in%c("3")&df$TME_cluster_T2%in%c("2")]=1
hot=df[df$hot==1,]
table(hot$group_hot,hot$pCR)
table(hot$group_hot,hot$DFS.status)
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_hot)+as.factor(LN), data=hot) ##final model
ShowRegTable(cox.test)
N_logistic <- glm(as.numeric(pCR) ~ as.factor(group_hot)+as.factor(LN)+as.factor(TS)+as.factor(IHC_subtype), data=hot) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)
df$IHC_subtype

library(readxl);library(tidyverse)
mydata<-read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/HR_OR_TME.xlsx", col_names=TRUE, sheet=1)



##########################################
##########################################
##########################################
df=df[!is.na(df$TME_cluster_T2),]
df=df[!is.na(df$TME3),]
df$cold=0
df$cold[df$TME_cluster_T2%in%c("1")]=1
table(df$cold)
df$group_cold=0
df$group_cold[df$TME_cluster_T2%in%c("1")&df$TME3%in%c("2","3")]=1
cold=df[df$cold=="1",]
table(cold$group_cold,cold$pCR)
table(cold$group_cold,cold$DFS.status)
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_cold), data=cold) ##final model
ShowRegTable(cox.test)
N_logistic <- glm(as.numeric(pCR) ~ as.factor(group_cold)+as.factor(LN)+as.factor(TS)+as.factor(IHC_subtype), data= cold) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)

df$hot=0
df$hot[df$TME_cluster_T2%in%c("2","3")]=1
table(df$hot)
df$group_hot=0
df$group_hot[df$TME_cluster_T2%in%c("3","2")&df$TME3%in%c("1")]=1
df$group_hot[df$TME_cluster_T2%in%c("3")&df$TME3%in%c("2")]=1
hot=df[df$hot==1,]
table(hot$group_hot,hot$pCR)
table(hot$group_hot,hot$DFS.status)
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_hot)+as.factor(LN), data=hot) ##final model
ShowRegTable(cox.test)
N_logistic <- glm(as.numeric(pCR) ~ as.factor(group_hot)+as.factor(LN)+as.factor(TS)+as.factor(IHC_subtype), data=hot) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)
df$IHC_subtype


