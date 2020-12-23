setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3")
library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
genelist=read_excel("gene_list.xlsx", col_names=TRUE, sheet=1)
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
exp <- as.matrix(exp[genelist$genelist,])
##生成聚类矩阵，并保存至setp3文件夹内
library(ConsensusClusterPlus)
results<-ConsensusClusterPlus(exp,maxK=3,reps = 1000, 
                            pItem = 0.8, pFeature = 1, title="output_cluster",
                            clusterAlg='km', distance="euclidean", seed=123456,
                            plot="pdf", #或"png"
                            corUse="pairwise.complete.obs",writeTable=T)
clu3=as.data.frame(results[[3]][["consensusClass"]])
colnames(clu3)="cluster3"
clu3$samplesID=row.names(clu3)
clu3$cluster=as.factor(clu3$cluster3)
pheno<- readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step2/tumor_purity.rds")
pheno=left_join(pheno,clu3,by="samplesID")
saveRDS(pheno,fil="tumor_pheno.rds")
