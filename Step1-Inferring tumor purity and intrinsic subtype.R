library(readr);library(IOBR);library(ggplot2);library(readxl);library(tidyverse);library(data.table);library(readr)
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step1")
library(estimate)
in.file <- system.file("extdata", "sample_input.txt", package="estimate")
filterCommonGenes("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step1/exp.txt",output.f="tumor_purity_10395genes.gct",id="GeneSymbol")
estimateScore("tumor_purity_10395genes.gct", "tumor_purity.gct")
#####################################
######combination clin data##########
#####################################
mydata=as.data.frame(fread("F:/Ph.D projects/0.Cohorts/PROMIX/as-provided/2020-12-01-KW/PROMIX Master File 4 Dec 2020_ki67imputed.csv"))
sub=as.data.frame(fread("F:/Ph.D projects/0.Cohorts/PROMIX/as-provided/2020-12-01-KW/2020-12-01-updated-KW.csv"))
sub=sub[,c("patientsID","subtype","ER","PR","grade","Tumorsize","Regional nodes BL")]
mydata=left_join(mydata,sub,by="patientsID")
mydata$IHC_subtype="Lum"
mydata$IHC_subtype[mydata$subtype=="TNBC"]="TNBC"
table(mydata$`PAM50 subtype_op`)

tumor_purity=as.data.frame(fread("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step1/tumor_purity.csv"))
samplesID=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/samplesID.csv"))
tumor_purity=left_join(tumor_purity,samplesID,by="samplesID")

tumor_purity=left_join(tumor_purity,mydata,by="patientsID")
tumor_purity$TumorPurity[(tumor_purity$`pCR or not`=="pCR")&(tumor_purity$tpt=="Surgery")]=0
#####################################
######Intrinsic——subtype############
###################################
library('clusterProfiler');library(genefu)
samplesID=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/samplesID.csv"))
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
annot<- bitr(exp$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
annot<- annot[!duplicated(annot$SYMBOL),]
colnames(annot)=c("probe","EntrezGene.ID")
rownames(annot)=annot$probe
row.names(exp)=exp[,1]
exp=exp[,-1]
exp=as.matrix(t(exp[annot$probe,]))
exp=scale(exp, center = TRUE, scale = TRUE)
## load SSP fitted in pam50####
data(pam50)
data(pam50.scale) #scale(x, center = TRUE, scale = TRUE)
data(pam50.robust)

Baseline=exp[samplesID$samplesID[samplesID$tpt=="Baseline"],]
Baseline=scale(Baseline, center = TRUE, scale = TRUE)
dim(Baseline)
pam50_baseline<- intrinsic.cluster.predict(sbt.model=pam50.scale,
                                          data=Baseline, annot=annot, do.mapping=T,
                                          do.prediction.strength=FALSE, verbose=TRUE)
table(pam50_baseline$subtype)

Cycle2=exp[samplesID$samplesID[samplesID$tpt=="Cycle 2"],]
Cycle2=scale(Cycle2, center = TRUE, scale = TRUE)
pam50_Cycle2<- intrinsic.cluster.predict(sbt.model=pam50.scale,
                                           data=Cycle2, annot=annot, do.mapping=T,
                                           do.prediction.strength=FALSE, verbose=TRUE)
table(pam50_Cycle2$subtype)

Surgery=exp[samplesID$samplesID[samplesID$tpt=="Surgery"],]
Surgery=scale(Surgery, center = TRUE, scale = TRUE)
pam50_Surgery<- intrinsic.cluster.predict(sbt.model=pam50.scale,
                                           data=Surgery, annot=annot, do.mapping=T,
                                           do.prediction.strength=FALSE, verbose=TRUE)
table(pam50_Surgery$subtype)

pam50_baseline=as.data.frame(pam50_baseline$subtype) 
colnames(pam50_baseline)="PAM50"

pam50_Cycle2=as.data.frame(pam50_Cycle2$subtype) 
colnames(pam50_Cycle2)="PAM50"

pam50_Surgery=as.data.frame(pam50_Surgery$subtype) 
colnames(pam50_Surgery)="PAM50"

pam50=rbind(pam50_baseline,pam50_Cycle2,pam50_Surgery)
pam50$samplesID=rownames(pam50)

tumor_purity=left_join(tumor_purity,pam50,by="samplesID") 
table(tumor_purity$PAM50)
tumor_purity$PAM50[(tumor_purity$`pCR or not`=="pCR")&(tumor_purity$tpt=="Surgery")]="Normal"

saveRDS(tumor_purity,fil="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
#####################################
###############Ending################
#####################################
