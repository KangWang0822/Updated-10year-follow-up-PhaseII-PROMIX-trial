#############################################
###################Figure2.A#################
#############################################
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2")
library(data.table);library(tidyverse);library(IOBR);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
table(mydata$TME_cluster)
row.names(mydata)=mydata$samplesID
mydata=mydata[order(mydata$TME_cluster),]
mydata$TME_cluster
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

mat1 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="1"]])   # baseline samples
TME_cluster1=mydata[mydata$samplesID%in%colnames(mat1),]
mat2 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="2"]])   # baseline samples
TME_cluster2=mydata[mydata$samplesID%in%colnames(mat2),]
mat3 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="3"]])   # baseline samples
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

table(mydata$TME_cluster)
table(pheno$PAM50)
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

#8*15
cell_fraction=c("samplesID","CD4_memory_T_cells","CD8_T_cells","Naive_B_cells","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells","Follicular_Helper_T_cells","Mast_cells","NK_cells")
res<-cell_bar_plot(input = mydata[TME_cluster1$samplesID,cell_fraction], title = "CIBERSORT Cell Fraction",id="samplesID")
res<-cell_bar_plot(input = mydata[TME_cluster2$samplesID,cell_fraction], title = "CIBERSORT Cell Fraction",id="samplesID")
res<-cell_bar_plot(input = mydata[TME_cluster3$samplesID,cell_fraction], title = "CIBERSORT Cell Fraction",id="samplesID")
#4*6
#######################################
##############Figure2.B################
#######################################
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")

mydata$pCR[mydata$pCR=="pCR"]=1
mydata$pCR[mydata$pCR=="non-pCR"]=0

mydata$LN=9
mydata$LN[mydata$`Regional nodes BL.y`=="Yes"]=1
mydata$LN[mydata$`Regional nodes BL.y`=="No"]=0

mydata$TS=9
mydata$TS[mydata$Tumorsize<50]=1
mydata$TS[mydata$Tumorsize>=50]=2

data=mydata[mydata$tpt%in%c("Baseline"),]  #Baseline,Cycle 2,Surgery
table(data$IHC_subtype)


table(data$TME_cluster,data$pCR)
table(data$TME_cluster)
table(data$TME_cluster,data$DFS.status)
table(data$TME_cluster,data$pCR)

cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(TME_cluster)+as.factor(LN)+as.factor(TS)+strata(IHC_subtype), data=data) ##final model
(test.ph <- cox.zph(cox.test))
step(cox.test)
ShowRegTable(cox.test)
cox.test <- coxph(Surv(OS.time,OS.status) ~ as.factor(TME_cluster)+as.factor(TS)+as.factor(LN)+as.factor(IHC_subtype), data=data) ##final model
ShowRegTable(cox.test)

N_logistic <- glm(as.numeric(pCR) ~ as.factor(TME_cluster)+as.factor(TS)+as.factor(LN)+as.factor(IHC_subtype), data= data) #注意Surgery应该是数值型变量##
step(N_logistic)
ShowRegTable(N_logistic)
################################################
#######################Forest###################
################################################
library(readxl);library(tidyverse)
mydata<-read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/HR_OR_TME.xlsx", col_names=TRUE, sheet=1)

ggplot(mydata,aes(HR,group))+
  geom_point(size=5, shape=18, color="#6699CC") +
  geom_errorbarh(aes(xmax = uci, xmin = lci,), height = 0.15,color="#6699CC") +
  geom_vline(xintercept = 1, linetype = "longdash",color="#CCCCCC",size=1) +
  scale_x_continuous(breaks = seq(0,14,1), labels = seq(0,14,1)) +
  labs(x="Multivarite HR/OR (95% CI)", y="") +theme_set(theme_bw())+
  theme(text = element_text(size=12),axis.line.x=element_line(size=1,colour="black"),panel.grid=element_blank(),panel.border=element_blank())
#5*7
ggtitle("Hazard/Odds Ratios for EFS/pCR")
#############################################
##############Figure2.B-C####################
#############################################
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ggplot2);library(sunburstR);library(ggpubr)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")

mydata$TME="Cold"
mydata$TME[mydata$TME_cluster=="2"]="Warm"
mydata$TME[mydata$TME_cluster=="3"]="Hot" 

mydata$burst=paste(mydata$tpt,mydata$TME,mydata$PAM50,sep="-")

library(RColorBrewer)
colKI <- c(
  "gray" = rgb(128, 128, 128, maxColorValue=256), 
  "plum" = rgb(135, 0, 82, maxColorValue=256), 
  "aqua" = rgb(151, 216, 218, maxColorValue=256), 
  "lavender" = rgb(189, 171, 179, maxColorValue=256), 
  "teal" = rgb(136, 196, 197, maxColorValue=256), 
  "cyklamen" = rgb(212, 9, 99, maxColorValue=256))

data=mydata%>%group_by(burst)%>%summarize(n = n())  ## 等价于summarize(n = sum(month))
color=brewer.pal(12,"Paired")
display.brewer.all(n=12)  #Paired
sunburst(data,
  count = TRUE,
  colors =color
)


############IHC_density#############
library(RColorBrewer);library(data.table);library(ggpubr);library(readxl);library(ComplexHeatmap);library(tidyverse);library(circlize)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
mydata=mydata[,c("samplesID","TME_cluster","tpt")]
density<-read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/Promix output_Multiplex_FINAL.xlsx",col_names=TRUE,sheet=2) #withoutNA
mydata=inner_join(density,mydata,by="samplesID")
mydata=as.data.frame(mydata[order(mydata$TME_cluster),])
row.names(mydata)=mydata$samplesID
All=c("CD4","CD8","CD68","T_regs","PDL1_CD4","PDL1_CD8","PDL1_CD68","PDL1_T_regs","PD1_CD4","PD1_CD8","PD1_CD68","PD1_T_regs")
Tumor=c("CD4_T","CD8_T","CD68_T","T_regs_T","PDL1_CD4_T","PDL1_CD8_T","PDL1_CD68_T","PDL1_T_regs_T","PD1_CD4_T","PD1_CD8_T","PD1_CD68_T","PD1_T_regs_T")
Stroma=c("CD4_S","CD8_S","CD68_S","T_regs_S","PDL1_CD4_S","PDL1_CD8_S","PDL1_CD68_S","PDL1_T_regs_S","PD1_CD4_S","PD1_CD8_S","PD1_CD68_S")
IHC=c(All,Tumor,Stroma)
TME_score=mydata[mydata$samplesID,IHC]
samplesID=row.names(TME_score)
TME_score=as.data.frame(lapply(TME_score,as.numeric))
row.names(TME_score)=samplesID
TME_score=log(TME_score+1)
#TME_score=scale(TME_score[,IHC])
TME_score=as.data.frame(t(TME_score))
#table(density$TME_cluster) #32 37 12 
mat1 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="1"]])   # baseline samples
TME_cluster1=mydata[mydata$samplesID%in%colnames(mat1),]
mat2 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="2"]])   # baseline samples
TME_cluster2=mydata[mydata$samplesID%in%colnames(mat2),]
mat3 = as.matrix(TME_score[,mydata$samplesID[mydata$TME_cluster=="3"]])   # baseline samples
TME_cluster3=mydata[mydata$samplesID%in%colnames(mat3),]

TME_function=c(rep("All",times=12),
               rep("Tumor",times=12),
               rep("Stroma",times=11))

pheno=TME_cluster1
ha1=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#6699FF"),labels = c("Cold (n=32)")),
                      Time=pheno$tpt,
                      col=list(Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      na_col="grey",
                      border=TRUE)

pheno=TME_cluster2
ha2=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#FFCC00"),labels = c("Warm (n=35)")),
                      Time=pheno$tpt,
                      col=list(Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      na_col="grey",
                      border=TRUE)


pheno=TME_cluster3
ha3=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#CC0000"),labels = c("Hot (n=12)")),
                      Time=pheno$tpt,
                      col=list(Time=c("Baseline"="floralwhite","Cycle 2"="gray","Surgery"="black")),
                      na_col="grey",
                      border=TRUE)
table(mydata$TME_cluster)

col_fun=colorRamp2(c(0, 1, 9), c("navy", "white", "firebrick3"))
ht_list=Heatmap(TME_function, name = "signature", show_row_names = FALSE, width = unit(10, "mm"),
                 col = structure(names = c("All","Tumor","Stroma"),
                                 brewer.pal(6,"Set2")),
                 row_split = factor(TME_function, levels = c("All","Tumor","Stroma")))+ 
  Heatmap(mat1, col = col_fun, name = "Scaled log(number of positive cells/mm2)",
          clustering_distance_columns = "spearman",
          show_row_dend = FALSE, show_column_dend = T,
          show_column_names = T,
          top_annotation = ha1,
          row_split = factor(TME_function, levels = c("All","Tumor","Stroma")), 
          row_title_gp = gpar(col = "#FFFFFF00"), width = unit(10, "cm"))+
  Heatmap(mat2,col=col_fun, show_column_names = T, 
          show_column_dend = T,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(10, "cm"))+
  Heatmap(mat3,col=col_fun, show_column_names = T, 
          show_column_dend = T,top_annotation=ha3,
          show_heatmap_legend = FALSE, width = unit(5, "cm"))

p0=draw(ht_list, row_title = paste0("EMT subtype"),
        annotation_legend_side = "left", heatmap_legend_side = "left") 

#(i.e. 225_op(cold), 504_0(warm), 316_0(hot))
##################################Box plot###############################
#########################################################################
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
mydata=mydata[,c("TME_cluster","StromalScore","ImmuneScore")]
mydata$StromalScore=log(mydata$StromalScore)
mydata$ImmuneScore=log(mydata$ImmuneScore)
mydata$StromalScore[is.na(mydata$StromalScore)]=mean(mydata$StromalScore)
mydata$ImmuneScore[is.na(mydata$ImmuneScore)]=mean(mydata$ImmuneScore)
data<-mydata%>% 
  ## 基因表达数据gather,gather的范围应调整
  gather(key="gene",value="Expression",StromalScore:ImmuneScore) %>% 
  ##
  dplyr::select(TME_cluster,gene,Expression,everything()) 
head(data)  ## 每个基因作为一个变量的宽数据
data$TME_cluster=as.factor(data$TME_cluster)

ggplot(data, aes(x=gene, y=Expression,fill=TME_cluster)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#6699FF","#FFCC00","#CC0000"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("")+ ylim(0,10)
#5*12################################
#################Barplot for percentage#############
####################################################
library(ggplot2);library(plyr)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
mydata$TIL=NA
mydata$TIL[mydata$tpt=="Baseline"]=mydata$TILs_baseline[mydata$tpt=="Baseline"]
mydata$TIL[mydata$tpt=="Cycle 2"]=mydata$TILs_cycle2[mydata$tpt=="Cycle 2"]
mydata$TIL[mydata$tpt=="Surgery"]=mydata$TILs_surgery[mydata$tpt=="Surgery"]
table(mydata$TIL,mydata$TME_cluster)
chisq.test(mydata$TIL,mydata$TME_cluster) #0.006
mydata=mydata[,c("TME_cluster","TIL")]
mydata$TIL=as.factor(mydata$TIL)
mydata=mydata[!is.na(mydata$TIL),]
mydata=mydata[order(mydata$TIL),]
mydata$value=1

ggplot(mydata, aes(x=TME_cluster,  y=value, fill=TIL))+
  geom_bar(stat="identity", position="fill")+scale_fill_brewer(palette="RdPu")+
  labs(x = "TME_cluster",y = "Percentage(%)")+theme_bw()+theme(element_blank())
#5*5
#################################################################################
#####################Biomarkers for cold warm hot tumor##########################
library(Blasso) 
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
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
TME=c("samplesID","TME_cluster",ADAPTIVE,INNATE,Stroma,Cytokine,TAAs,Checkpoint)
TME_score=mydata[,TME]
dat_clin=mydata[,c("samplesID","TME_cluster","DFS.status","DFS.time")]

colnames(TME_score)[1]="ID"
colnames(dat_clin)=c("ID","TME_cluster","status","time")

table(TME_score$TME_cluster)
output=matrix(nrow=0,ncol=2)

score=TME_score[TME_score$TME_cluster=="2",-2]
clin=dat_clin[dat_clin$TME_cluster=="2",-2]

res<-best_predictor_cox(target_data = clin, 
                        features = score, 
                        status = "status",
                        time = "time",
                        nfolds = 5,
                        permutation = 1000)
output=rbind(output,res$res)
write.csv(output,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/surv_freq.csv")
#########################################################################
#########################################################################
#visualization#
library(data.table);library(ggpubr);library(readxl)
mydata<-read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/surv_freq.xlsx",col_names=TRUE,sheet=2)
  ggdotchart(mydata, x="name", y="Freq", color = "group",
             palette = c("#6699FF", "#FFCC00", "#CC0000"),
             sorting = "descending", add = "segments", 
             add.params = list(color = "lightgray", size = 1),
             rotate = TRUE,group = "group", dot.size = 6,
             label = mydata$Freq, font.label = list(color="white",
             size=9, vjust=0.5), ggtheme = theme_pubr())
#7*7

#########################################################################
#########################################################################
#########################################################################
  ############Box plot#############
