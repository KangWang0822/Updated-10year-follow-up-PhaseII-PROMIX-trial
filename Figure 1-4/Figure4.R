library(networkD3);library(tidyr);library(tibble);library(dplyr)    # data manipulation
library(UpSetR);library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
data=data[order(data$TME_cluster),]
variable=c("TME_cluster","patientsID","pCR or not","DFS.status")
Baseline=data[data$tpt=="Baseline",variable]
colnames(Baseline)=c("TME1","patientsID","pCR1","DFS1")
Cycle2=data[data$tpt=="Cycle 2",variable]
colnames(Cycle2)=c("TME2","patientsID","pCR2","DFS2")
Surgery=data[data$tpt=="Surgery",variable]
colnames(Surgery)=c("TME3","patientsID","pCR3","DFS3")
dat1=full_join(Baseline,Cycle2,by="patientsID")
df=full_join(dat1,Surgery,by="patientsID")

df$TME1=as.character(df$TME1)
df$TME2=as.character(df$TME2)
df$TME3=as.character(df$TME3)
df=df[,c("TME1","TME2","TME3")]

################################################
###################T1->T2 T3####################
################################################
df1=df[(!is.na(df$TME2))&(!is.na(df$TME3)),]
table(df1$TME3)
# put your df in two columns, and preserve the ordering in many levels (columns) with paste0
links <- data.frame(source = c(paste0(df1$TME1,'_1'),paste0(df1$TME2,'_2')),
                    target   = c(paste0(df1$TME2,'_2'),paste0(df1$TME3,'_3')))
table(links$source)
links$color[links$source%in%c("1_1","1_2")]="cold"
links$color[links$source%in%c("2_1","2_2")]="warm"
links$color[links$source%in%c("3_1","3_2")]="hot"


my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "my_unique_group"]) .range(["#6699FF","#FFCC00","#CC0000"])'

# now convert as character
links$source <- as.character(links$source)
links$target<- as.character(links$target)
table(links$target)
links=links[links$target!="NA_3",]
nodes <- data.frame(name = unique(c(links$source,links$target)))
nodes$color[nodes$name%in%c("1_1","1_2","1_3")]="cold"
nodes$color[nodes$name%in%c("2_1","2_2","2_3")]="warm"
nodes$color[nodes$name%in%c("3_1","3_2","3_3")]="hot"
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
links$value <- 1 # add also a value

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',colourScale=my_color,
              Target = 'target', Value = 'value', NodeID = 'name',LinkGroup ="color",NodeGroup="color",
              height ="600",width ="800",nodeWidth ="50")

table(df1$TME3)
df1=df[(!is.na(df$TME1))&(!is.na(df$TME2)),]
########################################################
########################################################
######################Forest plot#######################
########################################################
library(networkD3);library(tidyr);library(tibble);library(dplyr) # data manipulation
library(UpSetR);library(ggplot2);library(plyr);library(gridExtra);library(ggpubr)
data=readRDS("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/step7/pheno_TME_MET.rds")
data=data[order(data$TME_cluster),]
variable=c("TME_cluster","patientsID","pCR or not","DFS.status")
Baseline=data[data$tpt=="Baseline",variable]
colnames(Baseline)=c("TME1","patientsID","pCR1","DFS1")
Cycle2=data[data$tpt=="Cycle 2",variable]
colnames(Cycle2)=c("TME2","patientsID","pCR2","DFS2")
Surgery=data[data$tpt=="Surgery",variable]
colnames(Surgery)=c("TME3","patientsID","pCR3","DFS3")
dat1=full_join(Baseline,Cycle2,by="patientsID")
df=full_join(dat1,Surgery,by="patientsID")

df$TME1=as.character(df$TME1)
df$TME2=as.character(df$TME2)
df$TME3=as.character(df$TME3)
df=df[,c("TME2","TME3")]

####T2->T3#####
df1=df[(!is.na(df$TME2))&(!is.na(df$TME3)),]
# put your df in two columns, and preserve the ordering in many levels (columns) with paste0
links=data.frame(source = c(paste0(df1$TME2,'_2')),
                    target   = c(paste0(df1$TME3,'_3')))
table(links$source)
links$color[links$source%in%c("1_2")]="cold"
links$color[links$source%in%c("2_2")]="warm"
links$color[links$source%in%c("3_2")]="hot"

my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "my_unique_group"]) .range(["#6699FF","#FFCC00","#CC0000"])'

# now convert as character
links$source <- as.character(links$source)
links$target<- as.character(links$target)
table(links$target)
nodes <- data.frame(name = unique(c(links$source,links$target)))
nodes$color[nodes$name%in%c("1_2","1_3")]="cold"
nodes$color[nodes$name%in%c("2_2","2_3")]="warm"
nodes$color[nodes$name%in%c("3_2","3_3")]="hot"
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
links$value <- 1 # add also a value

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',colourScale=my_color,
              Target = 'target', Value = 'value', NodeID = 'name',LinkGroup ="color",NodeGroup="color",
              height ="600",width ="800",nodeWidth ="50")

table(df1$TME2,df1$TME3)
df1=df[(!is.na(df$TME2))&(!is.na(df$TME3)),]
##############################################
################Forest plot###################
##############################################
library(readxl);library(tidyverse);library(ggplot2)
mydata<-read_excel("F:/Ph.D projects/1.PROMIX_10y_followup/dataprocessing/Figure4/HRandOR_forest.xlsx", col_names=TRUE, sheet=1)
ggplot(mydata,aes(OR,group))+
  geom_point(size=5, shape=18, color="#6699CC") +
  geom_errorbarh(aes(xmax = uci, xmin = lci,), height = 0.15,color="#6699CC") +
  geom_vline(xintercept = 1, linetype = "longdash",color="#CCCCCC",size=1) +
  scale_x_continuous(breaks = seq(0,14,1), labels = seq(0,14,1)) +
  labs(x="Multivarite HR/OR (95% CI)", y="") +theme_set(theme_bw())+
  theme(text = element_text(size=12),axis.line.x=element_line(size=1,colour="black"),panel.grid=element_blank(),panel.border=element_blank())
#5*7
ggtitle("Hazard/Odds Ratios for EFS/pCR")

mydata$group
