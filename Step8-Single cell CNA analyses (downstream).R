#################################################################
####################Extinction scRNA DGE#########################
#########################sig score calculation###################
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/diploid_prediction")
library(Seurat);library(stringr);library(IOBR)
#P126=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P01_KTN126_SNRS_TPM.txt",head=TRUE)
#P129=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P02_KTN129_SNRS_TPM.txt",head=TRUE)
P206=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P06_KTN206_SNRS_TPM.txt",head=TRUE)
P302=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P09_KTN302_SNRS_TPM.txt",head=TRUE)
mydata=cbind(P206,P302)
table(as.data.frame(str_split(colnames(mydata),"_",simplify = T))[,c(1,2)])
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
table(plate)
plate[plate==2]="P"
plate[plate==0]="B"
plate=as.data.frame(plate)
plate$n=1:nrow(plate)
colnames(plate)=c("time","n")
plate=plate[plate$time%in%c("B","P"),]

mydata=mydata[,plate$n]
mydata=log((mydata/10)+1)
saveRDS(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA_extinction.rds")
data=mydata

sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = data,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)
sig_kegg<-calculate_sig_score(pdata           = NULL,
                              eset            = data,
                              signature       = kegg,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_kegg=sig_kegg[,-c(1,2)]
sig_reactome<-calculate_sig_score(pdata           = NULL,
                                  eset            = data,
                                  signature       = reactome,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)
sig_reactome=sig_reactome[,-c(1,2)]
sig_collect=calculate_sig_score(pdata      = NULL,
                                eset            = data,
                                signature       = signature_collection,
                                method          = "ssgsea",
                                mini_gene_count = 2)
sig_collect=sig_collect[,-c(1,2)]

signature_score=cbind(sig_hallmark,sig_kegg,sig_reactome,sig_collect)
row.names(signature_score)=signature_score[,1]
signature_score=signature_score[,-c(1,2)]
signature_score=t(signature_score)
saveRDS(signature_score,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_extinction.rds")

library(Seurat);library(stringr);library(IOBR)
P102=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P10_KTN102_SNRS_TPM.txt",head=TRUE)
P132=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P11_KTN132_SNRS_TPM.txt",head=TRUE)
P152=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P14_KTN152_SNRS_TPM.txt",head=TRUE)
#P615=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P15_KTN615_SNRS_TPM.txt",head=TRUE)
mydata=cbind(P102,P132,P152)
table(as.data.frame(str_split(colnames(mydata),"_",simplify = T))[,c(1,2)])
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
table(plate)
plate[plate==0]="B"
plate=as.data.frame(plate)
plate$n=1:nrow(plate)
colnames(plate)=c("time","n")
plate=plate[plate$time%in%c("B","P"),]

mydata=mydata[,plate$n]
mydata=log((mydata/10)+1)
saveRDS(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA_persistence.rds")
data=mydata

sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = data,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)
sig_kegg<-calculate_sig_score(pdata           = NULL,
                              eset            = data,
                              signature       = kegg,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_kegg=sig_kegg[,-c(1,2)]
sig_reactome<-calculate_sig_score(pdata           = NULL,
                                  eset            = data,
                                  signature       = reactome,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)
sig_reactome=sig_reactome[,-c(1,2)]
sig_collect=calculate_sig_score(pdata      = NULL,
                                eset            = data,
                                signature       = signature_collection,
                                method          = "ssgsea",
                                mini_gene_count = 2)
sig_collect=sig_collect[,-c(1,2)]

signature_score=cbind(sig_hallmark,sig_kegg,sig_reactome,sig_collect)
row.names(signature_score)=signature_score[,1]
signature_score=signature_score[,-c(1,2)]
signature_score_=t(signature_score)
saveRDS(signature_score,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_persistence.rds")

#################################################################################
###########################differential ssGSEA score#############################
#################################################################################
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA_extinction.rds")
signature_score=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_extinction.rds")
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
signature_score=as.data.frame(signature_score)
signames=rownames(signature_score) 
signature_score=as.data.frame(lapply(signature_score,as.numeric))
rownames(signature_score)=signames
plate[plate==2]="P"
plate[plate==0]="B"
design<-model.matrix(~0+as.factor(plate))
colnames(design)<-c("baseline","op")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(signature_score, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
extinction = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,number=50000,sort.by = 'logFC')
extinction$absolute_diff=2^(extinction$logFC)*extinction$AveExpr
extinction_list=extinction[extinction$absolute_diff>0.2|extinction$absolute_diff<(-0.2),]
extinction_list$GSVA=row.names(extinction_list)

mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA_persistence.rds")
signature_score=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_persistence.rds")
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
signature_score=as.data.frame(signature_score)
signames=rownames(signature_score) 
signature_score=as.data.frame(lapply(signature_score,as.numeric))
rownames(signature_score)=signames
table(plate)
design<-model.matrix(~0+as.factor(plate))
colnames(design)<-c("baseline","op")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(signature_score, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
persistence = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,number=50000,sort.by = 'logFC')
persistence$absolute_diff=2^(persistence$logFC)*persistence$AveExpr
persistence_list=persistence[persistence$absolute_diff>0.2|persistence$absolute_diff<(-0.2),]
persistence_list$GSVA=row.names(persistence_list)
###输出2个genelist##
library(rJava);library(xlsxjars);library(xlsx);library(openxlsx)
mydata=list("extinction"=extinction_list,"persistence"=persistence_list)
write.xlsx(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/Diff_GSVA_scRNA.xlsx")





sig_IOBR=as.data.frame(sig_IOBR)
row.names(sig_IOBR)=sig_IOBR$ID
sig_IOBR=as.data.frame(sig_IOBR[,-c(1,2)])
sig_IOBR=as.data.frame(t(scale(sig_IOBR)))    
names=c("Fatty_Acid_Biosynthesis","Tyrosine_Metabolism","Tryptophan_Metabolism","Citric_Acid_Cycle","Cell_cycle","EMT1")
sig=sig_IOBR[names,]
row.names(sig)=c("FattyAcid","Tyrosine","Tryptophan","TAC","Cellcycle","EMT1")
#reorder dataframes
order=mydata[order(row.names(mydata)), ]
#single cells that did not express either GAPDH or ACTB
kept_cells=(order["GAPDH",])>0|(order["ACTB",]>0)
order=order[,kept_cells]
#remove genes in less than 30 cells
rm_genes= which(rowSums(order>0)<ncol(order)*0.3)
r2 =order[-rm_genes,]

genelis=c(row.names(r2))        ######加入代谢特征######
mydata=rbind(r2,sig)
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
obj= CreateSeuratObject(counts=mydata)
#find variable features to reduce run time
obj=FindVariableFeatures(obj)
#regress out batch
obj_nobatch = ScaleData(obj)
obj_nobatch@meta.data[,"protocol"]=plate
obj = RunPCA(obj_nobatch, features = genelis)
obj = FindNeighbors(obj)
obj = FindClusters(obj, resolution = 0.15)
obj = RunTSNE(obj)
DimPlot (obj , reduction = "tsne", label = TRUE)
DimPlot(object =obj, reduction = "tsne", group.by = "protocol", pt.size = 0.5)
FeaturePlot(obj,shape.by="protocol",features = c("EMT1"))











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

cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group)+as.factor(LN)+as.factor(TS)+strata(IHC_subtype), data=data) ##final model
(test.ph <- cox.zph(cox.test))
step(cox.test)
ShowRegTable(cox.test)
cox.test <- coxph(Surv(OS.time,OS.status) ~ as.factor(TME_cluster)+as.factor(TS)+as.factor(LN)+as.factor(IHC_subtype), data=data) ##final model
ShowRegTable(cox.test)

N_logistic <- glm(as.numeric(pCR) ~ as.factor(group)+as.factor(TS)+as.factor(LN)+as.factor(IHC_subtype), data= data) #注意Surgery应该是数值型变量##
step(N_logistic)
ShowRegTable(N_logistic)


library(tidyverse);library(SEERaBomb);library(ggsci)#load packages          
library(survival);library(survminer);library(bbmle)
gp=geom_point();gl=geom_line();jco=scale_color_jco()
dat <- data.frame(group = factor(c("Infliximab","Etanercept","Adalimumab","Golimumab","Certolizumab","Abatacept","Rituximab","Tocilizumab","Total"), 
                                 
                                 levels=c("Total","Certolizumab","Golimumab","Adalimumab","Etanercept","Infliximab","Tocilizumab","Rituximab","Abatacept")),
                  cen = c(1.88,2.67,2.35,2.14,5.08,1.68,2.11,2.01,2.16),
                  low = c(1.01,1.44,1.52,1.59,3.46,1.47,1.64,1.57,1.83),
                  high = c(3.51,4.94,3.65,2.89,7.48,1.90,2.72,2.57,2.55))



p_eff <- ggplot(dat,aes(cen,group)) +
  geom_point(size=5, shape=18) +
  geom_errorbarh(aes(xmax = high, xmin = low), height = 0.15) +
  geom_vline(xintercept = 1, linetype = "longdash") +
  scale_x_continuous(breaks = seq(0,14,1), labels = seq(0,14,1)) +
  labs(x="Favours Experimental", y="") +
  theme(text = element_text(size=20)) +
  ggtitle("Risk Ratio for ACR20 response at 6 months")

p_eff

dat2 <- data.frame(group = 
                     factor(c("Infliximab","Etanercept","Adalimumab","Golimumab","Certolizumab","Total"), 
                            
                            levels=c("Total","Certolizumab","Golimumab","Adalimumab","Etanercept","Infliximab")),
                   cen = c(3.22,0.71,1.59,0.98,2.72,1.26),
                   low = c(1.76,0.54,1.13,0.46,1.23,0.93),
                   high = c(5.91,0.92,2.23,2.08,6.01,1.71))

ggplot(dat2,aes(HR,group))+
  geom_point(size=5, shape=18, color="#6699CC") +
  geom_errorbarh(aes(xmax = uci, xmin = lci,), height = 0.15,color="#6699CC") +
  geom_vline(xintercept = 1, linetype = "longdash",color="#CCCCCC",size=1) +
  scale_x_continuous(breaks = seq(0,14,1), labels = seq(0,14,1)) +
  labs(x="Multivarite HR/OR (95% CI)", y="") +
  theme(text = element_text(size=20),panel.grid=element_blank(),panel.border=element_blank(),axis.line.x=element_line(size=1,colour="black"))
  ggtitle("Hazard/Odds Ratios for EFS/pCR")











mydata$group_metabolism[mydata$Metabolism_cluster=="3"]=1
mydata$group_metabolism[mydata$Metabolism_cluster=="1"]=2
mydata$group_metabolism[mydata$Metabolism_cluster=="2"]=3

mydata$pCR[mydata$pCR=="pCR"]=1
mydata$pCR[mydata$pCR=="non-pCR"]=0

mydata$LN=9
mydata$LN[mydata$`Regional nodes BL.y`=="Yes"]=1
mydata$LN[mydata$`Regional nodes BL.y`=="No"]=0

mydata$TS=9
mydata$TS[mydata$Tumorsize<50]=1
mydata$TS[mydata$Tumorsize>=50]=2

data=mydata[mydata$tpt=="Baseline",]  #Baseline,Cycle 2,Surgery

mydata$IHC_subtype
table(data$integrate_cluster,data$DFS.status)
table(data$TME_cluster,data$DFS.status)
#####################Metabolism###################
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_metabolism)+as.factor(TS)+as.factor(LN), data=data) ##final model
(test.ph <- cox.zph(cox.test))
step(cox.test)
ShowRegTable(cox.test)
cox.test <- coxph(Surv(OS.time,OS.status) ~ as.factor(group_metabolism)+as.factor(TS)+as.factor(LN), data=data) ##final model
ShowRegTable(cox.test)


N_logistic <- glm(as.numeric(pCR) ~ as.factor(group_metabolism)+as.factor(TS)+as.factor(LN), data= data) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)
###########################TME#############################
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(TME_cluster), data=data) ##final model
(test.ph <- cox.zph(cox.test))
step(cox.test)
ShowRegTable(cox.test)
cox.test <- coxph(Surv(OS.time,OS.status) ~ as.factor(TME_cluster)+as.factor(LN), data=data) ##final model
ShowRegTable(cox.test)

N_logistic <- glm(as.numeric(pCR) ~ as.factor(TME_cluster)+as.factor(TS)+as.factor(LN), data= data) #注意Surgery应该是数值型变量##
ShowRegTable(N_logistic)
########################################################

coxreg_1<- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_metabolism)+as.factor(TME_cluster)+as.factor(LN), data=data) ##DII_density_with_supp
coxreg_2 <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(group_metabolism)+as.factor(TME_cluster)+as.factor(Metabolism_cluster)*as.factor(TME_cluster)+as.factor(LN), data=data) ##DII_density_with_supp
anova(coxreg_2, coxreg_1)

fit=survfit(Surv(DFS.time,DFS.status)~Metabolism_cluster, data = data) 
ggsurvplot(fit,data=data, 
           title="Metabolism_cluster",
           ylab="Event-free Survival(EFS)", #"Overall Survival"
           xlab="Follow-up time(month)",
           pval = TRUE, 
           pval.method=T,
           conf.int="F",
           risk.table = TRUE,
           risk.table.col = "Metabolism_cluster",
           legend.title = "Metabolism_cluster",
           risk.table.y.text=F,
           xlim = c(0,120),
           break.time.by = 40,
           legend = "right",
           palette = "jco")

fit=survfit(Surv(OS.time,OS.status)~Metabolism_cluster, data = data) 
ggsurvplot(fit,data=data, 
           title="Metabolism_cluster",
           ylab="Overall survival(OS)", #"Overall Survival"
           xlab="Follow-up time(month)",
           pval = TRUE, 
           pval.method=T,
           conf.int="F",
           risk.table = TRUE,
           risk.table.col = "Metabolism_cluster",
           legend.title = "Metabolism_cluster",
           risk.table.y.text=F,
           xlim = c(0,120),
           break.time.by = 40,
           legend = "right",
           palette = "jco")

fit=survfit(Surv(DFS.time,DFS.status)~TME_cluster, data = data) 
ggsurvplot(fit,data=data, 
           title="TME_cluster",
           ylab="Event-free Survival(EFS)", #"Overall Survival"
           xlab="Follow-up time(month)",
           pval = TRUE, 
           pval.method=T,
           conf.int="F",
           risk.table = TRUE,
           risk.table.col = "TME_cluster",
           legend.title = "TME_cluster",
           risk.table.y.text=F,
           xlim = c(0,120),
           break.time.by = 40,
           legend = "right",
           palette = "jco")

fit=survfit(Surv(OS.time,OS.status)~TME_cluster, data = data) 
ggsurvplot(fit,data=data, 
           title="TME_cluster",
           ylab="Overall survival(OS)", #"Overall Survival"
           xlab="Follow-up time(month)",
           pval = TRUE, 
           pval.method=T,
           conf.int="F",
           risk.table = TRUE,
           risk.table.col = "TME_cluster",
           legend.title = "TME_cluster",
           risk.table.y.text=F,
           xlim = c(0,120),
           break.time.by = 40,
           legend = "right",
           palette = "jco")










cox.test <- coxph(Surv(DFS.time,DFS.status) ~ TK1_Change+strata(ER_status)+KI67+as.factor(LN)+as.factor(TS), data=TK1_C) ##DII_density_with_supp
cox.test <- coxph(Surv(DFS.time,DFS.status) ~ as.factor(TK1_Change)+as.factor(LN)+as.factor(TS)+as.factor(KI67)+strata(ER_status), data=TK1_C) ##DII_density_with_supp

step(cox.test)
ShowRegTable(cox.test)

mydata$

  




library(copykat)

P126=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/scRNA/Extinct_P01_KTN126_SNRS_TPM.txt",head=TRUE|FALSE)
copykat.test <- copykat(rawmat=P126, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, distance="euclidean")

P129=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/scRNA/Extinct_P02_KTN129_SNRS_TPM.txt",head=TRUE|FALSE)
copykat.test <- copykat(rawmat=P129, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, distance="euclidean")

