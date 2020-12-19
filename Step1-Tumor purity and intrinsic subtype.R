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

tumor_purity=as.data.frame(fread("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step1/tumor_purity.csv"))
samplesID=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/samplesID.csv"))
tumor_purity=left_join(tumor_purity,samplesID,by="samplesID")

tumor_purity=left_join(tumor_purity,mydata,by="patientsID")
saveRDS(tumor_purity,fil="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
#####################################
###############Ending################
#####################################
#####################################
