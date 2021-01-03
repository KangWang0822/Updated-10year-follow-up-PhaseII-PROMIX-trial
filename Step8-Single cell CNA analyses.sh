######
###############################
###############################
library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
pheno<- readRDS("tumor_purity.rds")
table(pheno$PAM50,pheno$IHC_subtype)
pheno$group="baseline"
pheno$group[pheno$tpt=="Cycle 2"]="cycle"
pheno$group[pheno$tpt=="Surgery"]="op"
#########################################
#############Pair-wise DGE###############
#########################################
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html
library(limma);library(Glimma);library(edgeR)
##op vs. baseline##
mat=exp[,pheno$samplesID[pheno$group%in%c("baseline","op")]]
clin_1=pheno[pheno$group%in%c("baseline","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.factor(clin_1$IHC_subtype)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("baseline","op","subtype","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
opVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
##op vs. cycle##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","op")]]
clin_1=pheno[pheno$group%in%c("cycle","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.factor(clin_1$IHC_subtype)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("cycle","op","subtype","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-cycle, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
opVScycle = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
##op vs. cycle##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","baseline")]]
clin_1=pheno[pheno$group%in%c("cycle","baseline"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.factor(clin_1$IHC_subtype)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("baseline","cycle","subtype","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cycle-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
cycleVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')

genelist=Reduce(union,list(cycleVSbase=row.names(cycleVSbase), 
                      opVSbase= row.names(opVSbase), 
                      opVScycle =row.names(opVScycle))) 

############################
############################
#############Lum############
############################
library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
pheno<- readRDS("tumor_purity.rds")
pheno$group="baseline"
pheno$group[pheno$tpt=="Cycle 2"]="cycle"
pheno$group[pheno$tpt=="Surgery"]="op"

#Lum#
exp=exp[,pheno$samplesID[pheno$IHC_subtype%in%c("Lum")]]
pheno=pheno[pheno$IHC_subtype%in%c("Lum"),]
##op vs. baseline##
mat=exp[,pheno$samplesID[pheno$group%in%c("baseline","op")]]
clin_1=pheno[pheno$group%in%c("baseline","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("baseline","op","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
opVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
##op vs. cycle##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","op")]]
clin_1=pheno[pheno$group%in%c("cycle","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("cycle","op","purity")
design
