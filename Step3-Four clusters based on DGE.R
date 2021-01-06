setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3")
library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
genelist=read_excel("gene_list.xlsx", col_names=TRUE, sheet=1)
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
pheno<- readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
row.names(exp)=exp[,1]
exp=exp[,-1]
exp <- as.matrix(exp[genelist$genelist,])
##生成聚类矩阵，并保存至setp3文件夹内
library(ConsensusClusterPlus)
exp=exp[,pheno$samplesID]
results<-ConsensusClusterPlus(exp,maxK=8,reps = 1000, 
                                  pItem = 0.8, pFeature = 1, title="output_cluster",
                                  clusterAlg='km', distance="euclidean", seed=123456,
                                  plot="pdf", #或"png"
                                  corUse="pairwise.complete.obs",writeTable=T)
clu=as.data.frame(results[[4]][["consensusClass"]])
colnames(clu)="cluster4"
clu$samplesID=row.names(clu)
pheno<- readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
pheno=left_join(pheno,clu,by="samplesID")
saveRDS(pheno,fil="tumor_pheno.rds")





pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")
pheno_lum=pheno[pheno$IHC_subtype=="Lum",]
exp_lum=exp[,pheno_lum$samplesID]
results_lum<-ConsensusClusterPlus(exp_lum,maxK=8,reps = 1000, 
                            pItem = 0.8, pFeature = 1, title="output_Lumcluster",
                            clusterAlg='km', distance="euclidean", seed=123456,
                            plot="pdf", #或"png"
                            corUse="pairwise.complete.obs",writeTable=T)
clu_lum=as.data.frame(results_lum[[4]][["consensusClass"]])
colnames(clu_lum)="cluster4"
clu_lum$samplesID=row.names(clu_lum)


pheno_TNBC=pheno[pheno$IHC_subtype=="TNBC",]
exp_TNBC=exp[,pheno_TNBC$samplesID]
results_TNBC<-ConsensusClusterPlus(exp_TNBC,maxK=8,reps = 1000, 
                              pItem = 0.8, pFeature = 1, title="output_TNBCcluster",
                              clusterAlg='km', distance="euclidean", seed=123456,
                              plot="pdf", #或"png"
                              corUse="pairwise.complete.obs",writeTable=T)
clu_TNBC=as.data.frame(results_TNBC[[4]][["consensusClass"]])
colnames(clu_TNBC)="cluster4"
clu_TNBC$samplesID=row.names(clu_TNBC)

cluster=rbind(clu_lum,clu_TNBC)

pheno<- readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
pheno=left_join(pheno,cluster,by="samplesID")
saveRDS(pheno,fil="tumor_pheno.rds")
#library(NMF)
#res <- nmf(exp,3, nrun=100)
#consensusmap(res)
rm(list=ls())




