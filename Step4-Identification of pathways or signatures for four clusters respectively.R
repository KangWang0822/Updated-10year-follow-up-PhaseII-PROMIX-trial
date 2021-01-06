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
###########################################
#################All#######################
###########################################
########Cluster1 features##################cell cycle/immune
library(limma);library(Glimma);library(edgeR);library(clusterProfiler)
table(pheno$cluster4)
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
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.01,lfc=0.5,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster1=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,minGSSize = 5,pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH"))
GSEA_cluster1$cluster="cluster1"
GSEA_cluster1$group="All"
########Cluster2 features##################metabolism
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
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster2=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,minGSSize = 5,pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH"))

GSEA_cluster2$cluster="cluster2"
GSEA_cluster2$group="All"
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
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster3=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,minGSSize = 5,pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH"))
GSEA_cluster3$cluster="cluster3"
GSEA_cluster3$group="All"
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
DGE = topTable(efit,adjust.method="BH",coef=1,p.value=0.05,lfc=0.5,number=50000,sort.by = 'logFC')
genelist=DGE%>%filter(logFC>0)
genelist=row.names(genelist)
GSEA_cluster4=as.data.frame(enricher(gene=genelist,TERM2GENE=GSEA,minGSSize = 5,pvalueCutoff = 0.05,
                                        pAdjustMethod = "BH"))

GSEA_cluster4$cluster="cluster4"
GSEA_cluster4$group="All"

GSEA_ALL=rbind(GSEA_cluster1,GSEA_cluster2,GSEA_cluster3,GSEA_cluster4)
################################
############Output##############
################################
library(rJava);library(xlsxjars);library(xlsx);library(openxlsx)
#mydata=list("GSEA_ALL"=GSEA_ALL,"GSEA_Lum"=GSEA_Lum,"GSEA_TNBC"=GSEA_TNBC)
write.xlsx(GSEA_ALL,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step4/GSEA_results.xlsx")
