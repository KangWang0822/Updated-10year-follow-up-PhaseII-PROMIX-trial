###############################
###############################
library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
pheno<- readRDS("tumor_purity.rds")
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
opVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
opVSbase_total= topTable(efit,adjust.method="BH",genelist=final_genelist,number=50000,sort.by = 'logFC')
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
opVScycle = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
opVScycle_total= topTable(efit,adjust.method="BH",genelist=final_genelist,number=50000,sort.by = 'logFC')
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
cycleVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
cycleVSbase_total = topTable(efit,adjust.method="BH",genelist=final_genelist,number=50000,sort.by = 'logFC')


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
opVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
##op vs. cycle##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","op")]]
clin_1=pheno[pheno$group%in%c("cycle","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("cycle","op","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-cycle, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
opVScycle = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
##cycle vs. baseline##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","baseline")]]
clin_1=pheno[pheno$group%in%c("cycle","baseline"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("baseline","cycle","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cycle-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
cycleVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')

genelist_lum=Reduce(union,list(cycleVSbase=row.names(cycleVSbase), 
                               opVSbase= row.names(opVSbase), 
                               opVScycle =row.names(opVScycle)))

##################
##################
#####TNBC#########
##################
##################
library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
pheno<- readRDS("tumor_purity.rds")
pheno$group="baseline"
pheno$group[pheno$tpt=="Cycle 2"]="cycle"
pheno$group[pheno$tpt=="Surgery"]="op"

#TNBC#
exp=exp[,pheno$samplesID[pheno$IHC_subtype%in%c("TNBC")]]
pheno=pheno[pheno$IHC_subtype%in%c("TNBC"),]
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
opVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
##op vs. cycle##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","op")]]
clin_1=pheno[pheno$group%in%c("cycle","op"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("cycle","op","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=op-cycle, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
opVScycle = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
##cycle vs. baseline##
mat=exp[,pheno$samplesID[pheno$group%in%c("cycle","baseline")]]
clin_1=pheno[pheno$group%in%c("cycle","baseline"),]
design<-model.matrix(~0+as.factor(clin_1$group)+as.numeric(clin_1$TumorPurity))
colnames(design)<-c("baseline","cycle","purity")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cycle-baseline, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(mat, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
cycleVSbase = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')

genelist_TN=Reduce(union,list(cycleVSbase=row.names(cycleVSbase), 
                              opVSbase= row.names(opVSbase), 
                              opVScycle =row.names(opVScycle)))

###输出4个genelist 三组比较的并集###
library(rJava);library(xlsxjars);library(xlsx);library(openxlsx)
final_genelist=Reduce(union,list(cycleVSbase=genelist, ####得到后再次做差异分析！！！！genelist=final_genelist########
                                 opVSbase=genelist_lum, 
                                opVScycle=genelist_TN))    
mydata=list("DGE_total"=final_genelist,"DGE_all"=genelist,"DGE_lum"=genelist_lum,"DGE_TNBC"=genelist_TN)
write.xlsx(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/gene_list.xlsx")

############################################
#######DGE list based on final_genelist#####
############################################
library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
exp=exp[final_genelist,]
pheno<- readRDS("tumor_purity.rds")
pheno$group="baseline"
pheno$group[pheno$tpt=="Cycle 2"]="cycle"
pheno$group[pheno$tpt=="Surgery"]="op"
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
opVSbase_total= topTable(efit,adjust.method="BH",number=50000,sort.by = 'logFC')
opVSbase_total$group="opVSbase"
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
opVScycle_total= topTable(efit,adjust.method="BH",number=50000,sort.by = 'logFC')
opVScycle_total$group="opVScycle"
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
cycleVSbase_total = topTable(efit,adjust.method="BH",number=50000,sort.by = 'logFC')
cycleVSbase_total$group="cycleVSbase"

mydata_total=rbind(opVSbase_total,opVScycle_total,cycleVSbase_total)
mydata_total$gene=row.names(mydata_total)
mydata=list("DGE"=mydata_total)
write.xlsx(mydata, file = "F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/Pair-wise-DGE.xlsx")
##############################################
######If RNA-seq voom conversion##############
##############################################
##########PROMIX-count#############
#####################################
library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2")
gse=readRDS(file="GSE87455_annotation.rds")
id_probe <- gse@gpls$GPL10558@dataTable@table
dim(id_probe)
head(id_probe)
id_probe[1:4,1:15]
View(head(id_probe))## you need to check this , which column do you need
probe2gene <- id_probe[,c("ID","Symbol")]

tmp=as.data.frame(read_tsv("GSE87455_non-normalized.txt"))
promix_exprs <- as.matrix(tmp[, grep("AVG_Signal", names(tmp), fixed = TRUE)])
rownames(promix_exprs) <- as.character(tmp[, "PROBE_ID"])
colnames(promix_exprs) <- gsub(".AVG_Signal", "", colnames(promix_exprs), fixed = TRUE)
## Assemble ExpressionSet
promix <- ExpressionSet(promix_exprs)
annotation(promix) <- "illuminaHumanv4"
## Check
stopifnot(validObject(promix))
## Clean up
remove(list = c("tmp"))
exprSet=promix_exprs

library(dplyr)
exprSet <- exprSet[rownames(exprSet) %in% probe2gene$ID,]
dim(exprSet)
exprSet[1:5,1:5]

#ids过滤探针
probe2gene <- probe2gene[match(rownames(exprSet),probe2gene$ID),]
dim(probe2gene)
probe2gene[1:5,1:2]

idcombine <- function(exprSet, probe2gene){
  tmp <- by(exprSet,
            probe2gene$Symbol,
            function(x) rownames(x)[which.max(rowMeans(x))])
  probes <- as.character(tmp)
  print(dim(exprSet))
  exprSet <- exprSet[rownames(exprSet) %in% probes,]
  
  print(dim(exprSet))
  rownames(exprSet) <- probe2gene[match(rownames(exprSet), probe2gene$ID),2]
  return(exprSet)
}

new_exprSet <- idcombine(exprSet,probe2gene)
new_exprSet[1:4,1:6]
dim(new_exprSet)

rownames(new_exprSet)
saveRDS(new_exprSet,fil="promix_exprSet_count.rds")
######################
######################
par(mfrow=c(1,2))
v <- voom(exp, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
##############################################
##############################################
##############################################
