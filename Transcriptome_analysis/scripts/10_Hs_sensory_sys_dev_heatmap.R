#This script generates 3 heatmaps - the first is for the full gene set using all the hydractinia genes found in this gene set. 
#The second is using only significant DEGs identified for this gene set. 
#The Third heatmap is of the genes that are expressed but not significantly differentially expressed

#Input files you will need for the first part (full gene set): 
  #Salmon output (the mapping dir we used before in edgeR)
  #libraries_to_stages.txt file that we used before in edgeR
  #The R output file from step 8 (contains gene_acc, Homo_seqid, OG, gene symbol, Actinula_seqid)
  #all_symbols_for_OGs output from 9.D_get_all_symbols_for_OGs.py

#Input files you will need for the second part (sig DEGs in gene set):
  #Salmon output (the mapping dir we used before in edgeR)
  #libraries_to_stages.txt file that we used before in edgeR
  #The R output file from step 8 (contains gene_acc, Homo_seqid, OG, gene symbol, Actinula_seqid)
  #all_symbols_for_OGs output from 9.D_get_all_symbols_for_OGs.py
  #The output string from 9.C_prep_sig_DEGs_for_heatmap.py 


library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("tximport")


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

setwd("~/Desktop/R /Hydractinia_transcriptomics_3-19-20")
system('ls')

dir <- getwd()
list.files()

#### Part 1: Make heatmap of all Hydractinia Genes in Gene Set ####

#library info file
devstages<-read.table("libraries_to_stages.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4

txi.salmon<- tximport(c(dev1files,dev2files,dev3files,dev4files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2,3,3,3,4,4,4))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#set up geneID map (hashtable, dictionary, lookup, etc)
homohydra.map<- read.table("Reduced_Hs_sensory_sys_dev_R_output.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/sens_sys_dev-all_symbols_for_OGs.txt")
names(all_symbols_for_OGs)[2] <- c("tot_symbol")
str(all_symbols_for_OGs)
all_symbols_for_OGs$OG<- as.character(all_symbols_for_OGs$OG)
all_symbols_for_OGs$tot_symbol<- as.character(all_symbols_for_OGs$tot_symbol)
str(all_symbols_for_OGs)

nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:5],all_symbols_for_OGs,by="OG")
str(homohydra.map)

#make new column for naming the transcripts
library(tidyr)
homohydra.map <-unite(homohydra.map, symbol_seqid, tot_symbol:Hydractinia_seqid, sep= "    ", remove=F, na.rm=FALSE )
head(homohydra.map)

#Keep only unique hydractinia headers
#library(dplyr)
length(homohydra.map$Hydractinia_seqid)
length(unique(homohydra.map$Hydractinia_seqid))

homohydra.map<- homohydra.map[!duplicated(homohydra.map$Hydractinia_seqid), ]
length(homohydra.map$Hydractinia_seqid)

keep<-unique(homohydra.map$Hydractinia_seqid); length(keep) 


#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(homohydra.map$symbol)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))
#summary(TFcts) #to figure out range

library(gplots)

#pdf("BestTFhits all and averaged.pdf",height=11,width=8)  
#using row as color break (normal)
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette(120),  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
#                                       change to symbol col  	
#specifying color break
heatmap.2(TFcts,  scale="none", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2,3,3,3,4,4,4));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break (normal)
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);

#heatmaps specifiying color break
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", notecex = 1 , labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(3,11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);



##### Clear objects from workspace #####



#### Part 2: Make heatmap of HydractiniaSig DEGS in Gene Set ####

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("~/Desktop/R /Hydractinia_transcriptomics_3-19-20")
system('ls')

dir <- getwd()
list.files()

#library info file
devstages<-read.table("libraries_to_stages.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4

txi.salmon<- tximport(c(dev1files,dev2files,dev3files,dev4files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2,3,3,3,4,4,4))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("Hs_planula_t.72185","Hs_planula_t.62380","Hs_planula_t.84960","Hs_planula_t.47977","Hs_planula_t.92541","Hs_planula_t.92543","Hs_planula_t.83120","Hs_planula_t.38478","Hs_planula_t.47382","Hs_planula_t.67842","Hs_planula_t.70828","Hs_planula_t.47814","Hs_planula_t.24836","Hs_planula_t.90763","Hs_planula_t.83711","Hs_planula_t.83712","Hs_planula_t.24723","Hs_planula_t.78038","Hs_planula_t.72185","Hs_planula_t.47554","Hs_planula_t.48675","Hs_planula_t.62380","Hs_planula_t.25179","Hs_planula_t.78332","Hs_planula_t.84960","Hs_planula_t.88759","Hs_planula_t.35598","Hs_planula_t.47977","Hs_planula_t.82643","Hs_planula_t.30195","Hs_planula_t.92543","Hs_planula_t.24912","Hs_planula_t.90552","Hs_planula_t.27004","Hs_planula_t.38528","Hs_planula_t.92454","Hs_planula_t.83120","Hs_planula_t.66373","Hs_planula_t.81533","Hs_planula_t.36699","Hs_planula_t.38478","Hs_planula_t.86065","Hs_planula_t.26803","Hs_planula_t.47734","Hs_planula_t.29846","Hs_planula_t.47477","Hs_planula_t.38054","Hs_planula_t.47256","Hs_planula_t.67841","Hs_planula_t.67842","Hs_planula_t.72997","Hs_planula_t.67607","Hs_planula_t.59560","Hs_planula_t.70828","Hs_planula_t.47814","Hs_planula_t.24836","Hs_planula_t.90763","Hs_planula_t.87581","Hs_planula_t.83711","Hs_planula_t.83712","Hs_planula_t.78038","Hs_planula_t.72185","Hs_planula_t.47554","Hs_planula_t.48675","Hs_planula_t.62380","Hs_planula_t.23286","Hs_planula_t.25179","Hs_planula_t.49134","Hs_planula_t.78332","Hs_planula_t.88759","Hs_planula_t.35598","Hs_planula_t.47977","Hs_planula_t.36723","Hs_planula_t.83005","Hs_planula_t.92543","Hs_planula_t.27004","Hs_planula_t.35455","Hs_planula_t.36173","Hs_planula_t.38528","Hs_planula_t.92454","Hs_planula_t.83120","Hs_planula_t.66373","Hs_planula_t.81533","Hs_planula_t.38478","Hs_planula_t.86065","Hs_planula_t.47734","Hs_planula_t.29846","Hs_planula_t.47477","Hs_planula_t.38054","Hs_planula_t.49570","Hs_planula_t.27144","Hs_planula_t.47256","Hs_planula_t.67841","Hs_planula_t.67842","Hs_planula_t.42186","Hs_planula_t.72997","Hs_planula_t.79337","Hs_planula_t.67607","Hs_planula_t.70828","Hs_planula_t.47814","Hs_planula_t.26123","Hs_planula_t.24836","Hs_planula_t.25243","Hs_planula_t.90763","Hs_planula_t.83711","Hs_planula_t.83712","Hs_planula_t.26403","Hs_planula_t.88569","Hs_planula_t.30195","Hs_planula_t.92541","Hs_planula_t.24883","Hs_planula_t.48878","Hs_planula_t.85060","Hs_planula_t.83120","Hs_planula_t.66373","Hs_planula_t.24883","Hs_planula_t.35461","Hs_planula_t.47477","Hs_planula_t.27144","Hs_planula_t.79337","Hs_planula_t.28757","Hs_planula_t.24724","Hs_planula_t.83120","Hs_planula_t.49187","Hs_planula_t.35461","Hs_planula_t.47518","Hs_planula_t.79337","Hs_planula_t.47814","Hs_planula_t.90763"
) #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Hydractinia_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Hydractinia_seqid<- as.character(sig_degs_df$Hydractinia_seqid)
str(sig_degs_df)

homohydra.map<- read.table("Reduced_Hs_sensory_sys_dev_R_output.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("~/Desktop/R /Hydractinia_transcriptomics_3-19-20/DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/sens_sys_dev-all_symbols_for_OGs.txt")
names(all_symbols_for_OGs)[2] <- c("tot_symbol")
str(all_symbols_for_OGs)
all_symbols_for_OGs$OG<- as.character(all_symbols_for_OGs$OG)
all_symbols_for_OGs$tot_symbol<- as.character(all_symbols_for_OGs$tot_symbol)
str(all_symbols_for_OGs)

nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:5],all_symbols_for_OGs,by="OG")
str(homohydra.map)

#make new column for naming the transcripts
library(tidyr)
homohydra.map <-unite(homohydra.map, symbol_seqid, tot_symbol:Hydractinia_seqid, sep= "    ", remove=F, na.rm=FALSE )

library(dplyr)
#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:7],sig_degs_df,by="Hydractinia_seqid")
nrow(sig_set)
length(unique(sig_set$Hydractinia_seqid))

#Keep only unique hydractinia headers
#library(dplyr)
unique_sig_set<- sig_set[!duplicated(sig_set$Hydractinia_seqid), ]
length(unique_sig_set$Hydractinia_seqid)

homohydra.map<-unique_sig_set
keep<-unique(homohydra.map$Hydractinia_seqid); length(keep) #1 col with 32 enteries 


#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(homohydra.map$symbol)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c("#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c("#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#pdf("BestTFhits all and averaged.pdf",height=11,width=8)  
#using row as color break
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette(120),  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
#                                       change to symbol col  	
#specifying color break
#heatmap.2(TFcts,  scale="none", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2,3,3,3,4,4,4));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 14),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);

#heatmaps specifying color break
#heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(3,11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);


##### Clear objects from workspace #####



#### Part 3: Make heatmap of Hydractinia Genes NOT Sig DEG (just expressed) #### 

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("/Users/Sydney/Desktop/R/Hydractinia_transcriptomics_3-19-20")
system('ls')

dir <- getwd()
list.files()

#4 groups

#library info file
devstages<-read.table("libraries_to_stages.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4

txi.salmon<- tximport(c(dev1files,dev2files,dev3files,dev4files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2,3,3,3,4,4,4))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#set up geneID map (hashtable, dictionary, lookup, etc) 
homohydra.map<- read.table("Reduced_Hs_sensory_sys_dev_R_output.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 
length(unique(homohydra.map$symbol)) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/sens_sys_dev-all_symbols_for_OGs.txt")
names(all_symbols_for_OGs)[2] <- c("tot_symbol")
str(all_symbols_for_OGs)
all_symbols_for_OGs$OG<- as.character(all_symbols_for_OGs$OG)
all_symbols_for_OGs$tot_symbol<- as.character(all_symbols_for_OGs$tot_symbol)
str(all_symbols_for_OGs)

nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:5],all_symbols_for_OGs,by="OG")
str(homohydra.map)
## make object of genes that are expressed but not sig de
not_sig_degs<-c("KDM5B","TBC1D20","CRYGB","MITF","CRYGA","FZR1","WNT2","HDAC1","LPCAT1","SMAD3","CHRDL1","ZDHHC16","ARHGAP35","PTF1A","LRP6","BAX","PROM1","CYP1A1","HDAC2","DLL4","P2RY12","SDHAF4","MYOM1","NF2","ATP2B4","AGTPBP1","SMARCA4","GATA3","FBN1","SMOC1","TENT2","WDPCP","TWIST1","NF1","ATF6","WNT5A","OPA1","ACTL6A","WNT16","TULP1","SHH","EPHA2","RING1","EPHB1","BMPR1B","PROX1","AHI1","TTLL5","B9D1","TBC1D32","ADAMTS18","PPP1R13L","DLL1","WNT2B","SDK2","RET","TMEM231","TRAF3IP1","RAB18","LHX1","HES5","UNC45B","EPHB2","MEIS1","CLCN2","FOXC1","PKNOX1","SDK1","NPHP4","TOPORS","IFT140","ATP2B1","OSR2","WNT5B","PRDM1","ATOH7","FGF9","PHACTR4","MFN2","CRB1","TH","NPHP1","FBN2","SLC25A25","BMP7","GRHL2","HIF1A","HCN1","TUB","GPD2","TENM3","ABCB5","SLC44A4","HSF4","CRYGD","BBS7","MAN2A1","TGFBR1","WNT7B","SIX6","ARL6","FGF10","BCL2","RXRA","WNT7A","ABI2","MED1","IFT122","CRYAB","TULP3","IMPDH2","DSCAM","IHH","DIO3","ARID1A","RDH13","WDR19","CDKN1B","RP1","SH3PXD2B","BMP6","RDH10","CYP1B1","LIMK2","JAG1","MYOM2","BAK1","FZD5","SMG9","CTNS","NEUROD1","CELF4","LRP5","PYGO2","INHBA","NIPBL","ACVRL1","KDM2B","MEIS2","PRKCI","CEP290","PSEN1","RPGRIP1L","GRHL3","NAGLU","ACVR2B","BCL11B","SOS1","MTNR1B","BBS4","ISL1","MYO7A","SIX3","CALB1","SRF","PFDN5"				
)

#make into a data frame and change name of column 
not_sig_degs_df<-data.frame(not_sig_degs)
str(not_sig_degs_df)
names(not_sig_degs_df)[1] <- c("symbol")


#join with homohydra map to get genes that are expressed but not sig DE expressed
nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:6],not_sig_degs_df,by="symbol")
str(homohydra.map)
length(unique(homohydra.map$symbol)) 


#make new column for naming the transcripts
library(tidyr)
homohydra.map <-unite(homohydra.map, symbol_seqid, tot_symbol:Hydractinia_seqid, sep= "   ", remove=F, na.rm=FALSE )
head(homohydra.map)

#Keep only unique hydractinia headers
library(dplyr)
length(homohydra.map$Hydractinia_seqid) 
length(unique(homohydra.map$Hydractinia_seqid)) 
length(unique(homohydra.map$symbol)) 
length(unique(homohydra.map$OG)) 
length(unique(homohydra.map$tot_symbol)) 

#homohydra.map<- homohydra.map[!duplicated(homohydra.map$Hydractinia_seqid), ]
#length(homohydra.map$Hydractinia_seqid)

keep<-unique(homohydra.map$Hydractinia_seqid); length(keep) 
length(unique(homohydra.map$tot_symbol)) 
length(unique(homohydra.map$symbol)) 

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) 
sort(homohydra.map$symbol)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
sort(rownames(TFcts)) 

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))
#summary(TFcts) #to figure out range

library(gplots)

#pdf("BestTFhits all and averaged.pdf",height=11,width=8)  
#using row as color break (normal)
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette(120),  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
#                                       change to symbol col  	
#specifying color break
heatmap.2(TFcts,  scale="none", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2,3,3,3,4,4,4));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break (normal)
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);

#heatmaps specifiying color break
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);




#dev.off()



#layout key (lmat)
# 1-Heatmap,
# 2-Row dendrogram,
# 3-Column dendrogram,
# 4-Key

