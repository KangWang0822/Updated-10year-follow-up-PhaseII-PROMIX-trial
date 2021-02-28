#######Icluster for immune phenotype#######
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6")
library(IOBR);library(ISOpureR);library(readr);library(yaml);library(ImmuneSubtypeClassifier);library(data.table)
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
res0 <- ImmuneSubtypeClassifier::callEnsemble(X = exp, geneids = 'symbol')
res0$subtype_immunity[res0$BestCall==1]="Wound healing"
res0$subtype_immunity[res0$BestCall==2]="IFN-γ dominant"
res0$subtype_immunity[res0$BestCall==3]="Inflammatory"
res0$subtype_immunity[res0$BestCall==4]="Lymphocyte depleted"
res0$subtype_immunity[res0$BestCall==6]="TGF-β dominant"
res=res0[,c("SampleIDs","subtype_immunity")]
colnames(res)=c("samplesID","subtype_immunity")
#################################################
sig_IOBR<-calculate_sig_score(pdata           = NULL,
                             eset            = exp,
                             signature       = signature_collection,
                             method          = "ssgsea",
                             mini_gene_count = 1)

colnames(sig_IOBR)[1]="samplesID"
#timer<-deconvo_tme(eset=exp,method="timer",group_list=rep("stad",dim(exp)[2]))
cibersort<-deconvo_tme(eset=exp,method ="cibersort",arrays = FALSE,perm=200)
colnames(cibersort)
cibersort$CD4_memory_T_cells=cibersort$T_cells_CD4_memory_resting_CIBERSORT+cibersort$T_cells_CD4_memory_activated_CIBERSORT+cibersort$T_cells_CD4_naive_CIBERSORT
cibersort$CD8_T_cells=cibersort$T_cells_CD8_CIBERSORT
cibersort$Naive_B_cells=cibersort$B_cells_naive_CIBERSORT
cibersort$Macrophages_M0=cibersort$Macrophages_M0_CIBERSORT
cibersort$Macrophages_M1=cibersort$Macrophages_M1_CIBERSORT
cibersort$Macrophages_M2=cibersort$Macrophages_M2_CIBERSORT
cibersort$Dendritic_cells=cibersort$Dendritic_cells_resting_CIBERSORT+cibersort$Dendritic_cells_activated_CIBERSORT
cibersort$Follicular_Helper_T_cells=cibersort$T_cells_follicular_helper_CIBERSORT
cibersort$Mast_cells=cibersort$Mast_cells_resting_CIBERSORT+cibersort$Mast_cells_activated_CIBERSORT
cibersort$NK_cells=cibersort$NK_cells_resting_CIBERSORT+cibersort$NK_cells_activated_CIBERSORT
cibersort$total_cells=cibersort$CD4_memory_T_cells+cibersort$CD8_T_cells+cibersort$Naive_B_cells+cibersort$Macrophages_M0+cibersort$Macrophages_M1+cibersort$Macrophages_M2+cibersort$Dendritic_cells+cibersort$Follicular_Helper_T_cells+cibersort$Mast_cells+cibersort$NK_cells

cibersort$CD4_memory_T_cells=cibersort$CD4_memory_T_cells/cibersort$total_cells
cibersort$CD8_T_cells=cibersort$CD8_T_cells/cibersort$total_cells
cibersort$Naive_B_cells=cibersort$Naive_B_cells/cibersort$total_cells
cibersort$Macrophages_M0=cibersort$Macrophages_M0/cibersort$total_cells
cibersort$Macrophages_M1=cibersort$Macrophages_M1/cibersort$total_cells
cibersort$Macrophages_M2=cibersort$Macrophages_M2/cibersort$total_cells
cibersort$Dendritic_cells=cibersort$Dendritic_cells/cibersort$total_cells
cibersort$Follicular_Helper_T_cells=cibersort$Follicular_Helper_T_cells/cibersort$total_cells
cibersort$Mast_cells=cibersort$Mast_cells/cibersort$total_cells
cibersort$NK_cells=cibersort$NK_cells/cibersort$total_cells

cell_fraction=c("ID","CD4_memory_T_cells","CD8_T_cells","Naive_B_cells","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells","Follicular_Helper_T_cells","Mast_cells","NK_cells")
cibersort=cibersort[,cell_fraction]
colnames(cibersort)[1]="samplesID"

#################################################
##########subtype_immunity_score#################
#################################################
exp=as.data.frame(fread("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/promix_non_normal.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]

load('comparative_immuneSigs_geneLists4.rda')

sigs1_2_geneIDs2<-as.character(na.omit(sigs1_2_eg2[[1]]$probe))
for(i in 2:length(sigs1_2_eg2)){
  sigs1_2_geneIDs2<-c(sigs1_2_geneIDs2,as.character(na.omit(sigs1_2_eg2[[i]]$probe)))
}
sigs1_2_geneIDs2<-unique(sigs1_2_geneIDs2)   ## 2652 unique
length(sigs1_2_geneIDs2)

# let's work on a subset of the data, and compare to the reported scores #

# first to remove some duplicate symbols '?' and 'SLC35E2'
geneSymbols=row.names(exp)
jdx <- which(geneSymbols %in% sigs1_2_geneIDs2)
gs2 <- geneSymbols[jdx]

# then make a smaller data set and give it row names
datSubset <- as.data.frame(exp[jdx,]) #  sample(x = 1:ncol(dat), size = 1000, replace = F)]) 

# just to have a look
dim(datSubset)
datSubset[1:5,1:5]
# apply the transforms #
# first the log2
datSubsetTransformed <- apply(datSubset, 2, function(a) log2(a+1))

dim(datSubsetTransformed)

# then the median scale
datSubsetTransformed <- t(apply(datSubsetTransformed, 1, function(a) a - median(a, na.rm=T)))

dim(datSubsetTransformed)

source('ImmuneSigs68_function.R')
scores<-ImmuneSigs_function(datSubsetTransformed, 
                              sigs1_2_eg2,
                              sigs12_weighted_means,
                              sigs12_module_weights,
                              sigs1_2_names2,
                              sigs1_2_type2)
	
scores=t(as.data.frame(scores[c("LIexpression_score","CSF1_response","Module3_IFN_score","TGFB_score_21050467","CHANG_CORE_SERUM_RESPONSE_UP"),]))
scores=as.data.frame(scores)
colnames(scores)=c("LIexpression_score_immunity","CSF1_response_immunity","Module3_IFN_score_immunity","TGFB_score_21050467_immunity","CHANG_CORE_SERUM_RESPONSE_UP_immunity")
scores$samplesID=as.character(row.names(scores))  

#################################
#################################
#################################
library(IOBR);library(reshape2)
##准备gmt file####
#mydata=as.data.frame(fread("sig_immunity.csv"))
#mydata<-dcast(mydata,Gene~SetName)
#dim(mydata)
#for(i in 1:3015){
#  for(j in 2:110){
#    if(!is.na(mydata[i,j]))
#      mydata[i,j]=mydata[i,1]
#  }
#}
#mydata=mydata[-which(is.na(mydata$Gene)),-1]
#gmt=matrix(nrow=109,ncol=1000)
#for(i in 1:109){gmt[i,1:length(mydata[!is.na(mydata[,i]),i])]=mydata[!is.na(mydata[,i]),i]}
#row.names(gmt)=as.character(colnames(mydata))
#write.csv(gmt,file="sig_immunity_gmt.csv")

library(GSVA);library(GSEABase)
geneSet=getGmt("sig_immunity.gmt")
score_immunity=gsva(as.matrix(exp),geneSet,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)#ssGSEA计算
score_immunity=t(as.data.frame(score_immunity))
score_immunity=as.data.frame(score_immunity)
score_immunity$samplesID=as.character(row.names(score_immunity))  

pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")

mydata<-pheno%>%inner_join(.,res,by="samplesID")%>%
  inner_join(.,scores,by="samplesID")%>%
  inner_join(.,cibersort,by="samplesID")%>%
  inner_join(.,score_immunity,by="samplesID")%>%
  inner_join(.,sig_IOBR,by="samplesID")
saveRDS(mydata,file="pheno_immune.rds")
####################################
#########TME_subtype################
####################################
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6")
library(iClusterPlus);library(ConsensusClusterPlus);library(NbClust);library(factoextra);library(tidyverse)
########check missing score########
#for(i in 1:length(Checkpoint)) {
#a=Checkpoint[i]
#print(mydata[,a])
#print(i)
#}
#Checkpoint[9]
##1.Immune cells##
mydata=readRDS("pheno_immune.rds")
row.names(mydata)=mydata$samplesID
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
cell_fraction=c("CD4_memory_T_cells","CD8_T_cells","Naive_B_cells","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells","Follicular_Helper_T_cells","Mast_cells","NK_cells")

#res<-cell_bar_plot(input = cibersort[,cell_fraction], title = "CIBERSORT Cell Fraction")
TME1=c(ADAPTIVE,INNATE,Stroma,Cytokine,TAAs,Checkpoint)
TME1=scale(mydata[,TME1])
TME2=as.matrix(mydata[,cell_fraction])
#Number of data types n.lambda
#1 any
#2 8, 13, 21, 34, 55, 89, 144, 233, 377, 610
#3 35, 101, 135, 185, 266, 418, 597, 828, 1010
#4 307, 562, 701, 1019, 2129, 3001, 4001, 5003, 6007
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Immune_cluster")
TME=cbind(TME1,TME2)
TME<-apply(TME,2,function(x){(x-mean(x))/sd(x)})
#数据集群性评估，使用get_clust_tendency()计算Hopkins统计量
res=get_clust_tendency(TME,40,graph=TRUE)
res$hopkins_stat #Hopkins统计量的值<0.5,表明数据是高度可聚合的
res$plot
res_nbclust<-NbClust(TME,distance = "euclidean",min.nc = 1,max.nc = 5,method = "complete",index = "all")
table(res_nbclust$Best.partition)

#TME<-apply(TME,2,function(x){(x-mean(x))/sd(x)})
results<-ConsensusClusterPlus(t(TME),maxK=5,reps = 1000, 
                              pItem = 0.8, pFeature = 1, title="TME_cluster",
                              clusterAlg='km', distance="euclidean", seed=123456,
                              plot="pdf", #或"png"
                              corUse="pairwise.complete.obs",writeTable=T)
clusters =as.data.frame(results[[2]][["consensusClass"]])
colnames(clusters)="TME_cluster"
clusters$samplesID=row.names(clusters)
table(clusters$TME_cluster)


set.seed(123) # 随机数种子保证每次结果一致性
before_run <- Sys.Date()
for(k in 1:5){
  file_name <- paste("cv.fit.k",k,".Rdata",sep="")
  if ( ! file.exists(file_name)){
    cv.fit  <-  tune.iClusterPlus(cpus=1,dt1=TME1, dt2=TME2,
                                  type=c("gaussian","gaussian"),
                                  K=k,n.lambda=89,
                                  scale.lambda=c(1,0.05),maxiter=20)
    save(cv.fit, file=file_name) #保存到Rdata文件里，便于中断后继续
  }
}
after_run <- Sys.Date()
paste0("Running time is ", after_run-before_run)

################
#读取保存的Rdata
################
output <- alist() # alist里的参数惰性求值
files <- grep("cv.fit",dir()) 
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]] <- cv.fit
}
# 提取之前设置的K和Lambda
nLambda <-  nrow(output[[1]]$lambda)
nK <-  length(output)
# 获取每个拟合结果中的BIC值
BIC <-  getBIC(output)
# 获取每个拟合结果中的
devR <-  getDevR(output)
# 寻找每个K中最小BIC对应lambda ID
minBICid  <-  apply(BIC,2,which.min)
# 计算最小BIC的
devRatMinBIC <-  rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] <-  devR[minBICid[i],i]
}
plot(1:(nK+1),c(0,devRatMinBIC),
     type="b",
     xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")

clusters <- getClusters(output)
rownames(clusters) <- rownames(mydata)
colnames(clusters) <- paste("K=",2:(length(output)+1),sep="")
head(clusters)

k <- 2
best.cluster <- clusters[,k]
best.fit <- output[[k]]$fit[[which.min(BIC[,k])]]

features <-  alist()
features[[1]] <-  colnames(TME1)
features[[2]] <-  colnames(TME2)
sigfeatures=alist()
for(i in 1:2){
  rowsum <- apply(abs(best.fit$beta[[i]]),1, sum)
  upper <- quantile(rowsum,prob=0.75)
  sigfeatures[[i]] <- (features[[i]])[which(rowsum>upper)]
}
names(sigfeatures) <- c("signature","cell_fraction")
#输出每个数据集的前几个特征
sigfeatures$signature
sigfeatures$cell_fraction

clusters=as.data.frame(clusters)
colnames(clusters)[2]="TME_cluster"
clusters$samplesID=row.names(clusters)
clusters=as.data.frame(clusters[,c("samplesID","TME_cluster")])
table(clusters$TME_cluster)

clusters=clusters[,c("samplesID","TME_cluster")]
####################################################################
#####################Pathological variables#########################
########################CD163,FOXP3,cellularity#####################
####################################################################
library(readxl);library(tidyverse);library(openxlsx)
Path=as.data.frame(readxl::read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Pathological data/FINAL SOS2_promix_tils_CD163_Cellularity_JL.xlsx"))
FOXP3=as.data.frame(readxl::read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Pathological data/170719 TIL CD163 FOXP3 biopsier promix 1-129.xlsx",sheet = 3)%>%select("glassseccd","glassnum","FOXP3"))
Path$glassnum=as.factor(Path$glassnum)
Path$glassseccd=as.factor(Path$glassseccd)
FOXP3$glassnum=as.factor(FOXP3$glassnum)
FOXP3$glassseccd=as.factor(FOXP3$glassseccd)
FOXP3_a=left_join(Path[Path$glassseccd=="a",],FOXP3[FOXP3$glassseccd=="a",],by="glassnum")%>%select("subjid","tptcd","tpt","TIL","CD163","Cellularity","glassseccd.x","glassnum","FOXP3")
FOXP3_b=left_join(Path[Path$glassseccd=="b",],FOXP3[FOXP3$glassseccd=="b",],by="glassnum")%>%select("subjid","tptcd","tpt","TIL","CD163","Cellularity","glassseccd.x","glassnum","FOXP3")
FOXP3_c=left_join(Path[Path$glassseccd=="c",],FOXP3[FOXP3$glassseccd=="c",],by="glassnum")%>%select("subjid","tptcd","tpt","TIL","CD163","Cellularity","glassseccd.x","glassnum","FOXP3")
FOXP3_d=left_join(Path[Path$glassseccd=="d",],FOXP3[FOXP3$glassseccd=="d",],by="glassnum")%>%select("subjid","tptcd","tpt","TIL","CD163","Cellularity","glassseccd.x","glassnum","FOXP3")
Path=rbind(FOXP3_a,FOXP3_b,FOXP3_c,FOXP3_d)%>%select("subjid","tpt","TIL","CD163","Cellularity","glassseccd.x","glassnum","FOXP3")
colnames(Path)=c("patientsID","tpt","TIL","CD163","Cellularity","glassseccd","glassnum","FOXP3")
write.xlsx(Path,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Pathological data/PROMIX_TIL_CD163_FOXP3.xlsx")


###################################################################
#########################Metabolic subtype#########################
###################################################################
##Step1 deconvolution TC and TAC#######
#######################################
setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster")
library(ISOpureR);library(clusterProfiler);library(readr);library(data.table);library(enrichplot);library(reshape2);library(tidyverse)
gmt=read.gmt("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster/metabolite.gmt")
uni_gene=unique(gmt$gene)
load("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/promix_non_normal.Rdata")
exp=new_exprSet
exp$gene=row.names(exp)
genelist=Reduce(intersect,list(exp$gene,uni_gene))
exp=exp[genelist,]
exp=exp[order(exp$gene),]
exp=exp[-1783,-276]
pheno=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step3/tumor_pheno.rds")
normal=exp[,pheno$samplesID[pheno$tpt=="Surgery"&pheno$`pCR or not`=="pCR"]]
tumor=exp[,!(colnames(exp)%in%colnames(normal))]
min(normal)
set.seed(123);
ISOpureS1model=ISOpure.step1.CPE(as.matrix(tumor),as.matrix(normal))
ISOpureS1model$alphapurities
saveRDS(ISOpureS1model,file="ISOpureS1model.rds")
# For reproducible results, set the random seed
set.seed(456);
# Run ISOpureR Step 2 - Patient Profile Estimation
ISOpureS2model=ISOpure.step2.PPE(as.matrix(tumor),as.matrix(normal),ISOpureS1model)
ISOpureS2model$cc_cancerprofiles
saveRDS(ISOpureS2model,file="ISOpureS2model.rds")



setwd("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster")
#####output FDR#####
gmt=read.gmt("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster/metabolite.gmt")
uni_gene=unique(gmt$gene)
exp=as.data.frame(fread("F:/Other projects/1.ribosome/data processing/datasets/PROMIX/new_exprSet.csv"))
row.names(exp)=exp[,1]
exp=exp[,-1]
cg=names(tail(sort(apply(exp,1,sd)),10000))  #########变化最大得2000gene

genelist=Reduce(intersect,list(row.names(exp),uni_gene))
genelist=Reduce(union,list(genelist,cg))
exp=as.data.frame(t(exp[genelist,]))
exp=as.data.frame(scale(exp))
Metabolite_group=c(rep(NA,4))
for(i in 1:275){
a=i
gene=t(exp[i,])
genelist=gene[,1]
names(genelist)=row.names(gene)
genelist=sort(genelist,decreasing=T)
GSEA=GSEA(genelist,TERM2GENE=gmt,minGSSize=1,maxGSSize=1000,pvalueCutoff=1,by="fgsea",nPerm=10000)
table=GSEA@result%>%select("ID","NES","p.adjust")
table$samplesID=row.names(exp)[i]
Metabolite_group=rbind(Metabolite_group,table)
print(paste0("Finish:",a))}
saveRDS(Metabolite_group[-1,],file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster/Metabolism_cluster10000.RDS")
#gseaplot2(GSEA,2) #可视化GSEA
####################################
####################################
####################################
Metabolite=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/Metabolism_cluster/Metabolism_cluster.RDS")

Metabolite$Amino_acid_group[Metabolite$ID=="Amino_acid"]="Neutral"
Metabolite$Amino_acid_group[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Amino_acid_group[Metabolite$ID=="Amino_acid"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Amino_acid_group[Metabolite$ID=="Amino_acid"])

Metabolite$Lipid_group[Metabolite$ID=="Lipid"]="Neutral"
Metabolite$Lipid_group[Metabolite$ID=="Lipid"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Lipid_group[Metabolite$ID=="Lipid"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Lipid_group[Metabolite$ID=="Lipid"])

Metabolite$Carbohydrate_group[Metabolite$ID=="Carbohydrate"]="Neutral"
Metabolite$Carbohydrate_group[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Carbohydrate_group[Metabolite$ID=="Carbohydrate"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Carbohydrate_group[Metabolite$ID=="Carbohydrate"])

Metabolite$TCA_group[Metabolite$ID=="TCA_cycle"]="Neutral"
Metabolite$TCA_group[Metabolite$ID=="TCA_cycle"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$TCA_group[Metabolite$ID=="TCA_cycle"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$TCA_group[Metabolite$ID=="TCA_cycle"])

Metabolite$Energy_group[Metabolite$ID=="Energy"]="Neutral"
Metabolite$Energy_group[Metabolite$ID=="Energy"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Energy_group[Metabolite$ID=="Energy"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Energy_group[Metabolite$ID=="Energy"])

Metabolite$Nucleotide_group[Metabolite$ID=="Nucleotide"]="Neutral"
Metabolite$Nucleotide_group[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Nucleotide_group[Metabolite$ID=="Nucleotide"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Nucleotide_group[Metabolite$ID=="Nucleotide"])

Metabolite$Vitamin_cofactor_group[Metabolite$ID=="Vitamin"]="Neutral"
Metabolite$Vitamin_cofactor_group[Metabolite$ID=="Vitamin"&Metabolite$p.adjust<0.3&Metabolite$NES>0]="Upregulated"
Metabolite$Vitamin_cofactor_group[Metabolite$ID=="Vitamin"&Metabolite$p.adjust<0.3&Metabolite$NES<0]="Downregulated"
table(Metabolite$Vitamin_cofactor_group[Metabolite$ID=="Vitamin"])

variable=c("samplesID","NES","Amino_acid_group","Lipid_group","Carbohydrate_group","TCA_group","Energy_group","Nucleotide_group","Vitamin_cofactor_group")

Metabolite_group=Metabolite[Metabolite$ID=="Lipid",c("samplesID","Lipid_group")]%>%left_join(Metabolite[Metabolite$ID=="Amino_acid",c("samplesID","Amino_acid_group")],by="samplesID")%>%
  left_join(Metabolite[Metabolite$ID=="Carbohydrate",c("samplesID","Carbohydrate_group")],by="samplesID")%>%left_join(Metabolite[Metabolite$ID=="TCA_cycle",c("samplesID","TCA_group")],by="samplesID")%>%
  left_join(Metabolite[Metabolite$ID=="Energy",c("samplesID","Energy_group")],by="samplesID")%>%left_join(Metabolite[Metabolite$ID=="Nucleotide",c("samplesID","Nucleotide_group")],by="samplesID")%>%
  left_join(Metabolite[Metabolite$ID=="Vitamin",c("samplesID","Vitamin_cofactor_group")],by="samplesID")
table(Metabolite_group$Amino_acid_group)
Metabolite_group$Amino_acid_group[is.na(Metabolite_group$Amino_acid_group)]="Neutral"

###################################################################
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step6/pheno_immune.rds")
row.names(mydata)=mydata$samplesID
finaldata=left_join(mydata,clusters,by="samplesID")%>%left_join(Metabolite_group,by="samplesID")
saveRDS(finaldata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure2/pheno_TME_MET.rds")
saveRDS(finaldata,file="F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")



