########################################################
########################################################
########################################################
library(maftools)
install.packages('devtools')
library(devtools)
install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
library(clonevol)
data(aml1)
x <- aml1$variants




set.seed(33)
dat = data.frame(Subject = 1:10, 
                 Months = sample(4:20, 10, replace=TRUE),
                 Treated=sample(0:1, 10, replace=TRUE),
                 Stage = sample(1:4, 10, replace=TRUE),
                 Continued=sample(0:1, 10, replace=TRUE))

dat = dat %>%
  group_by(Subject) %>%
  mutate(Complete=sample(c(4:(max(Months)-1),NA), 1, 
                         prob=c(rep(1, length(4:(max(Months)-1))),5), replace=TRUE),
         Partial=sample(c(4:(max(Months)-1),NA), 1, 
                        prob=c(rep(1, length(4:(max(Months)-1))),5), replace=TRUE),
         Durable=sample(c(-0.5,NA), 1, replace=TRUE))

# Order Subjects by Months
dat$Subject = factor(dat$Subject, levels=dat$Subject[order(dat$Months)])

# Melt part of data frame for adding points to bars
dat.m = melt(dat %>% select(Subject, Months, Complete, Partial, Durable),
             id.var=c("Subject","Months"))


library(copynumber)
data(lymphoma)
head(lymphoma)
sub.lymphoma <- subsetData(data=lymphoma,sample=1:3)
lymph.wins <- winsorize(data=sub.lymphoma,verbose=FALSE)
single.seg <- pcf(data=lymph.wins,gamma=12,verbose=FALSE)
single.seg <- pcf(data=lymph.wins,gamma=12,verbose=FALSE)
plotGenome(data=sub.lymphoma,segments=single.seg,sample=1,cex=3)

#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 6 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:6)

#Winsorize data:
wins.data <- winsorize(data=sub.lymphoma,return.outliers=TRUE)

#Use pcf to find segments:        
uni.segments <- pcf(data=wins.data,gamma=12)

#Use multipcf to find segments as well:
multi.segments <- multipcf(data=wins.data,gamma=12)

#Plot data and pcf-segments over entire genome for all six samples (one page
#for each sample):
plotGenome(data=sub.lymphoma,segments=uni.segments)

#Let each sample define its own range, and adjust range to fit all observations:
plotGenome(data=sub.lymphoma,segments=uni.segments,equalRange=FALSE,q=0)
BiocManager::install("timescape")
library(timescape)
library(Rsamtools)
library(GenomicAlignments)
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGAlignments(ex1_file)

subject <- granges(galn)[1]
subject
## Note the absence of query no. 9 (i.e. 'galn[9]') in this result:
as.matrix(findOverlaps(galn, subject))

## This is because, by default, findOverlaps()/countOverlaps() are
## strand specific:
galn[8:10]
countOverlaps(galn[8:10], subject)
countOverlaps(galn[8:10], subject, ignore.strand=TRUE)

## Count alignments in 'galn' that DO overlap with 'subject' vs those
## that do NOT:
table(overlapsAny(galn, subject))
## Extract those that DO:
subsetByOverlaps(galn, subject)

## GAlignmentsList
galist <- GAlignmentsList(galn[8:10], galn[3000:3002])
gr <- GRanges(c("seq1", "seq1", "seq2"), 
              IRanges(c(15, 18, 1233), width=1),
              strand=c("-", "+", "+"))

countOverlaps(galist, gr)
countOverlaps(galist, gr, ignore.strand=TRUE)
findOverlaps(galist, gr)
findOverlaps(galist, gr, ignore.strand=TRUE)
########################################################
########################################################
########################################################
######################scRNA extinction##################
library(Seurat);library(stringr);library(rafalib);library(clustree);library(tidyverse)
library(RColorBrewer)
colKI <- c("gray" = rgb(128, 128, 128, maxColorValue=256), 
  "plum" = rgb(135, 0, 82, maxColorValue=256), 
  "aqua" = rgb(151, 216, 218, maxColorValue=256), 
  "lavender" = rgb(189, 171, 179, maxColorValue=256), 
  "teal" = rgb(136, 196, 197, maxColorValue=256), 
  "cyklamen" = rgb(212, 9, 99, maxColorValue=256))
#P126=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/scRNA/Extinct_P01_KTN126_SNRS_TPM.txt",head=TRUE|FALSE)
#P129=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/scRNA/Extinct_P02_KTN129_SNRS_TPM.txt",head=TRUE|FALSE)
P206=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P06_KTN206_SNRS_TPM.txt",head=TRUE|FALSE)
P302=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P09_KTN302_SNRS_TPM.txt",head=TRUE|FALSE)
mydata=P302 #mydata=P302
gsva=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_extinction.rds")
sig_extinction=c("CellCycle_Reg","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                  "HALLMARK_ADIPOGENESIS","HALLMARK_GLYCOLYSIS","Citric_Acid_Cycle",
                 "REACTOME_THE_FATTY_ACID_CYCLING_MODEL","KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
                 "HALLMARK_G2M_CHECKPOINT","Purine_Biosynthesis")
gsva=as.data.frame(gsva[sig_extinction,colnames(mydata)])
row.names(gsva)=c("CellCycle","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","TAC",
                  "FATTYACIDCYCLING","UNSATURATEDFATTYACIDS","G2M","Purine")
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
table(plate)
rm_genes= which(rowSums(mydata>0)<ncol(mydata)*0.2)
mydata=mydata[-rm_genes,]
all.genes=rownames(mydata)
mydata=rbind(mydata,gsva)
obj_206<-CreateSeuratObject(counts=mydata, project = "P206", min.cells = 5)
obj_206<- NormalizeData(obj_206,verbose = FALSE)
obj_206<- FindVariableFeatures(obj_206, selection.method = "vst")
obj_206<- ScaleData(obj_206,features =all.genes)
obj_206<- RunPCA(obj_206, features = all.genes)
# Examine and visualize PCA results a few different ways
print(obj_206[["pca"]], dims = 1:5, nfeatures = 5)
obj_206<- RunTSNE(obj_206, reduction = "pca", dims = 1:20)
obj_206<- FindNeighbors(obj_206,reduction = "pca", dims = 1:20)
obj_206<- FindClusters(obj_206,resolution = 0.5)
obj_206$type=plate
colKI
DimPlot(obj_206,reduction="tsne",group.by ="type",cols =c("#D30963","#808080"),pt.size = 1) #5*5
FeaturePlot(obj_206, features = c("CellCycle","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","TAC",
            "FATTYACIDCYCLING","UNSATURATEDFATTYACIDS","G2M","Purine"),min.cutoff = "q10", pt.size = 1)
#15*10
###################################################
##################scRNA persistence#################
###################################################
library(Seurat);library(stringr);library(rafalib);library(clustree);library(tidyverse)
library(RColorBrewer)

P152=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P14_KTN152_SNRS_TPM.txt",head=TRUE|FALSE)
P132=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P11_KTN132_SNRS_TPM.txt",head=TRUE|FALSE)
P102=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Persist_P10_KTN102_SNRS_TPM.txt",head=TRUE|FALSE)

mydata=P132 #mydata=P302
gsva=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_persistence.rds")
sig_extinction=c("HALLMARK_MYC_TARGETS_V1","HALLMARK_G2M_CHECKPOINT","REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS",
                 "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_ADIPOGENESIS",
                 "HALLMARK_GLYCOLYSIS","Purine_Biosynthesis","REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE_",
                 "Pyrimidine_Biosynthesis","KEGG_CYSTEINE_AND_METHIONINE_METABOLISM")
gsva=as.data.frame(gsva[sig_extinction,colnames(mydata)])
row.names(gsva)=c("MYC","G2M","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","Purine",
                  "TCA","Pyrimidine","METHIONINE")
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
table(plate)
rm_genes= which(rowSums(mydata>0)<ncol(mydata)*0.2)
mydata=mydata[-rm_genes,]
all.genes=rownames(mydata)
mydata=rbind(mydata,gsva)
obj_206<-CreateSeuratObject(counts=mydata, project = "P206", min.cells = 5)
obj_206<- NormalizeData(obj_206,verbose = FALSE)
obj_206<- FindVariableFeatures(obj_206, selection.method = "vst")
obj_206<- ScaleData(obj_206,features =all.genes)
obj_206<- RunPCA(obj_206, features = all.genes)
# Examine and visualize PCA results a few different ways
print(obj_206[["pca"]], dims = 1:5, nfeatures = 5)
obj_206<- RunTSNE(obj_206, reduction = "pca", dims = 1:20)
obj_206<- FindNeighbors(obj_206,reduction = "pca", dims = 1:20)
obj_206<- FindClusters(obj_206,resolution = 0.5)
obj_206$type=plate
DimPlot(obj_206,reduction="tsne",group.by ="type",cols =c("#D30963","#808080"),pt.size = 1) #5*5
FeaturePlot(obj_206, features = c("MYC","G2M","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","Purine",
                                  "TCA","Pyrimidine","METHIONINE"),min.cutoff = "q10", pt.size = 1)
#15*10






#############integrated extinction#################
#P126=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P01_KTN126_SNRS_TPM.txt",head=TRUE|FALSE)
#P129=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P02_KTN129_SNRS_TPM.txt",head=TRUE|FALSE)
#obj_126<-CreateSeuratObject(counts=P126, project = "P126", min.cells = 5)
#obj_126<- NormalizeData(obj_126,verbose = FALSE)
#obj_126<- FindVariableFeatures(obj_126, selection.method = "vst", nfeatures = 2000)
#obj_129<-CreateSeuratObject(counts=P129, project = "P129", min.cells = 5)
#obj_129<- NormalizeData(obj_129,verbose = FALSE)
#obj_129<- FindVariableFeatures(obj_129, selection.method = "vst", nfeatures = 2000)
P206=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P06_KTN206_SNRS_TPM.txt",head=TRUE|FALSE)
P302=read.table("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/scRNA/Extinct_P09_KTN302_SNRS_TPM.txt",head=TRUE|FALSE)
mydata=cbind(P206,P302)
rm_genes= which(rowSums(mydata>0)<ncol(mydata)*0.2)
mydata=mydata[-rm_genes,]
all.genes=rownames(mydata)
P206=P206[all.genes,]
P302=P302[all.genes,]
plate=str_split(colnames(mydata),"_",simplify = T)[,2]
table(plate)
plate[plate=="P"]="2"
gsva=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step8/sigscore_scRNA_extinction.rds")
sig_extinction=c("CellCycle_Reg","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                 "HALLMARK_ADIPOGENESIS","HALLMARK_GLYCOLYSIS","Citric_Acid_Cycle",
                 "REACTOME_THE_FATTY_ACID_CYCLING_MODEL","KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
                 "HALLMARK_G2M_CHECKPOINT","Purine_Biosynthesis")
gsva=as.data.frame(gsva[sig_extinction,colnames(mydata)])
row.names(gsva)=c("CellCycle","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","TAC",
                  "FATTYACIDCYCLING","UNSATURATEDFATTYACIDS","G2M","Purine")
P206=rbind(P206,gsva[,colnames(P206)])
P302=rbind(P302,gsva[,colnames(P302)])

obj_206<-CreateSeuratObject(counts=P206,project="P206",min.cells=5)
obj_206<- NormalizeData(obj_206,verbose= FALSE)
obj_206<- FindVariableFeatures(obj_206,selection.method="vst",nfeatures=2000)

obj_302<-CreateSeuratObject(counts=P302,project="P302",min.cells=5)
obj_302<- NormalizeData(obj_302,verbose = FALSE)
obj_302<- FindVariableFeatures(obj_302, selection.method = "vst", nfeatures = 2000)

tumor.anchors=FindIntegrationAnchors(object.list = list(obj_206,obj_302), dims = 1:20)
tumor.combined <- IntegrateData(anchorset = tumor.anchors, dims = 1:20)

DefaultAssay(tumor.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
tumor.combined<- ScaleData(tumor.combined, verbose = FALSE)
tumor.combined<- RunPCA(tumor.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
tumor.combined<- RunTSNE(tumor.combined, reduction = "pca", dims = 1:20)
tumor.combined<- FindNeighbors(tumor.combined,reduction = "pca", dims = 1:20)
tumor.combined<- FindClusters(tumor.combined,resolution = 0.5)
tumor.combined$type=plate
tumor.combined <- RunUMAP(tumor.combined, dims = 1:10)
# Visualization
DimPlot(tumor.combined, reduction = "tsne", group.by ="type",cols =c("#D30963","#808080"),pt.size = 1)
tumor.combined$orig.ident=tumor.combined$type
table(tumor.combined$orig.ident)
####################Dot######################
B_obj=subset(tumor.combined,idents="0")
avg.B <- log1p(AverageExpression(B_obj, verbose = FALSE)$RNA)
avg.B$gene <- rownames(avg.B)
colnames(avg.B)[1]="B"

P_obj=subset(tumor.combined,idents="2")
avg.P <- log1p(AverageExpression(P_obj, verbose = FALSE)$RNA)
avg.P$gene <- rownames(avg.P)
colnames(avg.P)[1]="P"

dot=left_join(avg.B,avg.P,by="gene")
  
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(dot, aes(B,P)) + geom_point() + ggtitle("extinction")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
################Violin########################
plots <- VlnPlot(tumor.combined, features = c("CellCycle","EMT","FATTYACID","CHOLESTEROL","ADIPOGENESIS","GLYCOLYSIS","TAC",
      "FATTYACIDCYCLING","UNSATURATEDFATTYACIDS","G2M","Purine"), group.by = "type",pt.size = 0,cols =c("#D30963","#808080"),combine = F)
CombinePlots(plots = plots, ncol = 11)
##################Heatmap#####################
# find all markers of cluster 1
tumor.combined$orig.ident=tumor.combined$type
table(tumor.combined$orig.ident)
cluster0.markers <- FindMarkers(tumor.combined,ident.1 =0,min.pct = 0.25)
head(cluster0.markers, n = 5)
cluster2.markers <- FindMarkers(tumor.combined,ident.1 =2,min.pct = 0.25)
head(cluster2.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(tumor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#############################################
#############################################
BiocManager::install("deepSNV")
library(deepSNV)
regions <- data.frame(chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 2074, stop=3585)
data(HIVmix) # Attach the data instead, as it could fail in routine checks without internet connection.
show(HIVmix)
control(HIVmix)[100:110,]
plot(HIVmix)
SNVs <- summary(HIVmix, sig.level=0.05, adjust.method="BH")
head(SNVs)



library(ISOpureR)
path.to.data <- paste0(file.path(system.file(package = "ISOpureR"), 'extdata', 'Beer'));
load(file.path(path.to.data , 'beer.tumordata.250.transcripts.30.patients.RData'));
load(file.path(path.to.data , 'beer.normaldata.250.transcripts.RData'));
beer.tumordata=log(beer.tumordata+1)
beer.normaldata=log(beer.normaldata+1)
ISOpureS1model <- ISOpure.step1.CPE(beer.tumordata,beer.normaldata)
ISOpureS1model$alphapurities
ISOpureS2model <- ISOpure.step2.PPE(
  beer.tumordata,
  beer.normaldata,
  ISOpureS1model
);
cancer=ISOpureS2model$cc_cancerprofiles
TAC=ISOpure.calculate.tac(beer.tumordata,ISOpureS2model$cc_cancerprofiles,ISOpureS2model$alphapurities)
# Check what the cancer profiles look like
str(ISOpureS2model$cc_cancerprofiles);
## num [1:250, 1:30] 203 1464 153 1021 194 ...
# Look at the first entries in the cancer profile for a particular patient,
# say patient 3
head(ISOpureS2model$cc_cancerprofiles[ ,3]);
