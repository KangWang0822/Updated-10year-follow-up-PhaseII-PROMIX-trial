library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
#Lipid Amino_acid TCA Nucleotide

chisq.test(mydata$TME_cluster,mydata$Vitamin_cofactor_group)

library(ggpubr);library(ggplot2);library(ggpie)
mydata$
#Lipid_group Amino_acid_group Carbohydrate_group TCA_group Energy_group Nucleotide_group Vitamin_cofactor_group
  ggpie(mydata,Vitamin_cofactor_group,TME_cluster, 
        nrow=1,                    # number of rows
        border.color="white",      # border color
        border.width=1.5,          # border width
        label.color="black")+
  scale_fill_manual(values=c("#99CCFF","#CCCCCC","#FF6699"))
#6*6
########################Barplot for percentage##########################
library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
data$value=1
#Lipid_group Amino_acid_group TCA_group Energy_group Carbohydrate_group Nucleotide_group Vitamin_cofactor_group
metabolite="Vitamin_cofactor_group"
mydata=data[data$tpt=="Baseline",]
chisq.test(mydata$TME_cluster,mydata[,metabolite])
ggplot(mydata, aes(x=TME_cluster,  y=value, fill=mydata[,metabolite]))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=c("#99CCFF","#CCCCCC","#FF6699"))+
  labs(x = "TME_cluster",y = "Percentage(%)")+theme_bw()+theme(element_blank())
mydata=data[data$tpt=="Cycle 2",]
chisq.test(mydata$TME_cluster,mydata[,metabolite])
ggplot(mydata, aes(x=TME_cluster,  y=value, fill=mydata[,metabolite]))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=c("#99CCFF","#CCCCCC","#FF6699"))+
  labs(x = "TME_cluster",y = "Percentage(%)")+theme_bw()+theme(element_blank())
mydata=data[data$tpt=="Surgery",]
chisq.test(mydata$TME_cluster,mydata[,metabolite])
ggplot(mydata, aes(x=TME_cluster,  y=value, fill=mydata[,metabolite]))+
  geom_bar(stat="identity",position="fill")+scale_fill_manual(values=c("#99CCFF","#CCCCCC","#FF6699"))+
  labs(x = "TME_cluster",y = "Percentage(%)")+theme_bw()+theme(element_blank())
#6*8 10%
############################################################
########################Network#############################
######################prepare OR_table######################
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
mydata=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")

mydata$pCR[mydata$pCR=="pCR"]=1
mydata$pCR[mydata$pCR=="non-pCR"]=0

mydata$LN=9
mydata$LN[mydata$`Regional nodes BL.y`=="Yes"]=1
mydata$LN[mydata$`Regional nodes BL.y`=="No"]=0

mydata$TS=9
mydata$TS[mydata$Tumorsize<50]=1
mydata$TS[mydata$Tumorsize>=50]=2
data=mydata
#data=mydata[mydata$tpt%in%c("Baseline"),]  #Baseline,Cycle 2,Surgery
table(data$IHC_subtype)

#TME_cluster Lipid_group Amino_acid_group TCA_group Energy_group Carbohydrate_group Nucleotide_group Vitamin_cofactor_group
logit=glm(as.numeric(pCR)~as.factor(Vitamin_cofactor_group)+as.factor(TS)+as.factor(LN)+as.factor(IHC_subtype), data= data) #注意Surgery应该是数值型变量##
ShowRegTable(logit)
#########################################
#################Visualization##########
library(readxl);library(reshape2);library(corrplot);library(plyr);library(igraph);library(RColorBrewer)
bb=read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure3/Network/figure3.B/OR_input.xlsx", sheet = 1)
bb$weight <- abs(log10(bb$pvalue))
#用HR标圆心点的颜色
bb$weight_HR <- (as.numeric(bb$OR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
poscol <- "#FB9A99" #正相关用红色连线
negcol <- "#C6DBEF" #负相关用蓝色连线
mycol <-brewer.pal(8,'Dark2')#cluster的颜色，如果有更多类，就给更多的颜色

ID=c("TME_cluster","Lipid_group","Amino_acid_group","TCA_group","Energy_group","Carbohydrate_group","Nucleotide_group","Vitamin_cofactor_group")
input_data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
input_data=input_data[,ID]
input_data[input_data=="Downregulated"]=-1
input_data[input_data=="Neutral"]=0
input_data[input_data=="Upregulated"]=1
input_data$Lipid_group=as.numeric(input_data$Lipid_group)
input_data$Amino_acid_group=as.numeric(input_data$Amino_acid_group)
input_data$TCA_group=as.numeric(input_data$TCA_group)
input_data$Energy_group=as.numeric(input_data$Energy_group)
input_data$Carbohydrate_group=as.numeric(input_data$Carbohydrate_group)
input_data$Nucleotide_group=as.numeric(input_data$Nucleotide_group)
input_data$Vitamin_cofactor_group=as.numeric(input_data$Vitamin_cofactor_group)

corr <- cor(input_data, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 3, #hclust聚为4类，根据数据的具体情况调整
         mar = c(4,0,4,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)

cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(input_data) 

#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)
corpvlue <- corpvlue[!corpvlue$cor==1,]
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)
#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)
#保存到文件，便于查看
#write.csv(corpvlue,"F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure3/Network/figure3.B/output_links.csv")

cellcluster <- as.data.frame(t(input_data))
#cellcluster[1:5,1:5]

hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus8 <- cutree(hc, 8)) #分4类

A <- as.character(rownames(as.data.frame(subset(clus8,clus8==1))))
B <- as.character(rownames(as.data.frame(subset(clus8,clus8==2))))
C <- as.character(rownames(as.data.frame(subset(clus8,clus8==3))))
D <- as.character(rownames(as.data.frame(subset(clus8,clus8==4))))
E <- as.character(rownames(as.data.frame(subset(clus8,clus8==5))))
f <- as.character(rownames(as.data.frame(subset(clus8,clus8==6))))
G <- as.character(rownames(as.data.frame(subset(clus8,clus8==7))))
H <- as.character(rownames(as.data.frame(subset(clus8,clus8==8))))
cls <- list(A,B,C,D,E,f,G,H)

nodes <- as.data.frame(unlist(cls))
nodes$type <- c(rep("A",1),rep("B",1),rep("C",1),rep("D",1),rep("E",1),rep("f",1),rep("G",1),rep("H",1))
names(nodes) <- c("media","type.label")
nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes
summary(nodes$media %in% bb$ID) #检查细胞名是否一致
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge

nodes$Fraction <- abs(nodes$weight_HR)
nodes$id <- paste("S", 1:8, sep = "")
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")

corpvlue$from <- revalue(corpvlue$from,c('TME_cluster'='S1','Lipid_group'='S2','Amino_acid_group'='S3','TCA_group'='S4','Energy_group'='S5','Carbohydrate_group'='S6','Nucleotide_group'='S7','Vitamin_cofactor_group'='S8'))

corpvlue$to <- revalue(corpvlue$to,c('TME_cluster'='S1','Lipid_group'='S2','Amino_acid_group'='S3','TCA_group'='S4','Energy_group'='S5','Carbohydrate_group'='S6','Nucleotide_group'='S7','Vitamin_cofactor_group'='S8'))
(links <- corpvlue)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2],"C"=mycol[3],"D"=mycol[4],"E"=mycol[5],"f"=mycol[6],"G"=mycol[7],"H"=mycol[8]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*10 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重
pdf("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure3/Network/figure3.B/network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
#legend("topright", #图例的位置
#       c("Cell cluster-A", "Cell cluster-B", "Cell cluster-C", "Cell cluster-D"),
#       pch=21, col="black", pt.bg=mycol, pt.cex=3,
#       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.1, 0.05, 0.01)
s <- sqrt(abs(log10(f)))*10
legend("topright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "Mutivariate Logit regression test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)
dev.off()

#############################################
#############Correlation heatmap#############
#############################################
library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
row.names(data)=data$samplesID
library(IOBR)
names(signature_metabolism)
names(signature_tumor)
names(signature_collection)
TME=c("B_cell_PCA_16704732","CD8 T cells","CD4_memory_T_cells","Treg cells",
      "T helper cells","Macrophages_M1","Macrophages_M2","Dendritic_cells",
      "NK_cells","PD1_data","PDL1_data","Lymph vessels","Angiogenesis")

Amino_acid=c("Homocysteine_Biosynthesis","Methionine_Cycle","Glycine_Serine_and_Threonine_Metabolism",
             "Phenylalanine_Metabolism","Arginine_Biosynthesis","Cysteine_and_Methionine_Metabolism","Histidine_Metabolism",
             "Tyrosine_Metabolism","Tryptophan_Metabolism","Beta_Alanine_Metabolism",
             "Taurine_and_Hypotaurine_Metabolism","Glutathione_Metabolism")
Lipid=c("Fatty_Acid_Biosynthesis","Steroid_Biosynthesis","Steroid_Hormone_Biosynthesis","Cholesterol_Biosynthesis","Biosynthesis_of_Unsaturated_Fatty_Acids","Glycerolipid_Metabolism","Cardiolipin_Metabolism","Ether_Lipid_Metabolism","Linoleic_Acid_Metabolism","Glycerophospholipid_Metabolism")
TCA=c("Citric_Acid_Cycle","Pyruvate_Metabolism","Pantothenate_and_CoA_Biosynthesis","Oxidative_Phosphorylation")
Metabolite=c(Lipid,Amino_acid,TCA)

TME_data=data[,TME]
TME_data=scale(TME_data)
Metabolite_data=(data[,Metabolite])
Metabolite_data=as.data.frame(t(scale(Metabolite_data)))

immuscore <- function(Metabolite){
  y <- as.numeric(Metabolite_data[Metabolite,])
  colnames <- colnames(TME_data)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(TME_data[,x]),y,type="spearman")
    data.frame(factor1=Metabolite,factor2=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
immuscore("Homocysteine_Biosynthesis")

data<- do.call(rbind,lapply(Metabolite,immuscore))
head(data)
#保存到文件
#write.csv(data,"F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure3/Network/figure3.C/correlation.csv", quote = F, row.names = F)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]
library(ggplot2);library(dplyr)
ggplot(data, aes(x=factor(factor2,levels = c(unique(factor2))),
                 y=factor(factor1,levels = c(unique(factor1))))) + 
  geom_tile(aes(fill = cor),size=1)+
  scale_fill_gradient2(low = "#99CCFF",mid = "white",high = "#FF6699")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legend
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
#7*7 80%
###########################################################
#########################Upset figure######################
###########################################################
library(UpSetR);library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
data$TME=0
data$TME[data$TME_cluster%in%c("1","2")]=1
data$Lipid=0
data$Lipid[data$Lipid_group=="Upregulated"]=1
data$Amino=0
data$Amino[data$Amino_acid_group=="Upregulated"]=1
data$Carbohydrate=0
data$Carbohydrate[data$Carbohydrate_group=="Upregulated"]=1
data$TCA=0
data$TCA[data$TCA_group=="Upregulated"]=1
data$Energy=0
data$Energy[data$Energy_group=="Upregulated"]=1
data$Nucleotide=0
data$Nucleotide[data$Nucleotide_group=="Upregulated"]=1
data$Vitamin=0
data$Vitamin[data$Vitamin_cofactor_group=="Upregulated"]=1

upset(data, sets = c("Lipid_group","Amino_acid_group","Carbohydrate_group"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on") ######效果不好
################################################################

