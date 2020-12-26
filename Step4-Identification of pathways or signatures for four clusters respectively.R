###############################################
############寻找cluster的特征富集##############
###############################################
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4")
library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")
#genelist=read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/gene_list.xlsx", col_names=TRUE, sheet=1)
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
exp=exp[,pheno$samplesID]

library(msigdbr);library(clusterProfiler)
hallmark=msigdbr(species = "Homo sapiens", category = c("H"))%>%select(gs_name,gene_symbol)
kegg=msigdbr(species = "Homo sapiens", category = c("C2"),subcategory = "CP:KEGG")%>%select(gs_name,gene_symbol)
reactome=msigdbr(species = "Homo sapiens", category = c("C2"),subcategory = "CP:REACTOME")%>%select(gs_name,gene_symbol)
GSEA=rbind(hallmark,kegg,reactome)
immun=msigdbr(species="Homo sapiens",category=c("C7"))%>%select(gs_name,gene_symbol)

#exp <- as.data.frame(exp[genelist$genelist,])


#table(pheno$group)
#table(pheno$cluster4)

#cg=names(tail(sort(apply(exp,1,sd)),5000))  #########变化最大得5000gene
#exp=exp[cg,]
#write.csv(exp,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4/GSEA/exp.csv")
#write.csv(pheno,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4/GSEA/pheno.csv")
########Cluster1 features##################immune
library(limma);library(Glimma);library(edgeR);library(clusterProfiler)
pheno$group="cluster"
pheno$group[pheno$cluster4%in%c("2","3","4")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster1=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH"))
GSEA_cluster1_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                     pAdjustMethod = "BH"))
########Cluster2 features##################代谢
pheno$group="cluster"
pheno$group[pheno$cluster4%in%c("1","3","4")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster2=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH"))
GSEA_cluster2_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))
########Cluster3 features##################EMT
pheno$group="cluster"
pheno$group[pheno$cluster4%in%c("1","2","4")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster3=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH"))
GSEA_cluster3_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))
########Cluster4 features##################cell-cycle
pheno$group="cluster"
pheno$group[pheno$cluster4%in%c("1","2","3")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster4=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH"))
GSEA_cluster4_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))

library(rJava);library(xlsxjars);library(xlsx);library(openxlsx)
GSEA_results=Reduce(union,list(cluster1=GSEA_cluster1, ####得到后再次做差异分析！！！！genelist=final_genelist########
                                 cluster2=GSEA_cluster2, 
                                 cluster3=GSEA_cluster3,
                                 cluster4=GSEA_cluster4))  
################################
############Output##############
################################
mydata=list("cluster1"=GSEA_cluster1,"cluster2"=GSEA_cluster2,"cluster3"=GSEA_cluster3,"cluster4"=GSEA_cluster4)
write.xlsx(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4/GSEA_results.xlsx")

mydata=list("cluster1"=GSEA_cluster1_immun,"cluster2"=GSEA_cluster2_immun,"cluster3"=GSEA_cluster3_immun,"cluster4"=GSEA_cluster4_immun)
write.xlsx(mydata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4/GSEA_immune_results.xlsx")


##########################################################
###################3 cluster version######################
##########################################################
########Cluster2 features##################EMT
library(limma);library(Glimma);library(edgeR);library(clusterProfiler)
pheno$group="cluster"
pheno$group[pheno$cluster3%in%c("2","3")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster1=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.01,
                                     pAdjustMethod = "BH"))
GSEA_cluster1_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))
########Cluster2 features##################代谢
pheno$group="cluster"
pheno$group[pheno$cluster3%in%c("1","3")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster2=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.01,
                                     pAdjustMethod = "BH"))
GSEA_cluster2_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))
########Cluster3 features##################EMT
pheno$group="cluster"
pheno$group[pheno$cluster3%in%c("1","2")]="control"
design<-model.matrix(~0+as.factor(pheno$group)+as.numeric(pheno$TumorPurity)+as.factor(pheno$IHC_subtype))
colnames(design)<-c("cluster","control","purity","subtype")
design

contr.matrix <- makeContrasts(
  op.vs.baseline=cluster-control, 
  levels=colnames(design))
contr.matrix

vfit <- lmFit(exp, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.58,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster3=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,pvalueCutoff = 0.01,
                                     pAdjustMethod = "BH"))
GSEA_cluster3_immun=as.data.frame(enricher(gene=genelist,TERM2GENE=immun,pvalueCutoff = 0.01,
                                           pAdjustMethod = "BH"))

###############################
###############################














