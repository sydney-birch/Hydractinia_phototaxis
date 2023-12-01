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

setwd("/Users/Sydney/Desktop/R/Hydractinia_transcriptomics_3-19-20")
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
homohydra.map<- read.table("Reduced_Hs_Sensory_Percep_Light_Stim_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 
length(unique(homohydra.map$symbol)) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("sens_percep_light_stim-all_symbols_for_OGs.txt") #Reduced
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
homohydra.map <-unite(homohydra.map, symbol_seqid, tot_symbol:Hydractinia_seqid, sep= "   ", remove=F, na.rm=FALSE )
head(homohydra.map)

#Keep only unique hydractinia headers
library(dplyr)
length(homohydra.map$Hydractinia_seqid)
length(unique(homohydra.map$Hydractinia_seqid))

homohydra.map<- homohydra.map[!duplicated(homohydra.map$Hydractinia_seqid), ]
length(homohydra.map$Hydractinia_seqid)

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

#heatmaps specifiying color break
#heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);






##### Clear objects from workspace #####



#### Part 2: Make heatmap of HydractiniaSig DEGS in Gene Set ####

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("/Users/Sydney/Desktop/R/Hydractinia_transcriptomics_3-19-20")
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


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formatted in this way
sig_degs<-c("Hs_planula_t.15325","Hs_planula_t.35573","Hs_planula_t.48076","Hs_planula_t.62380","Hs_planula_t.84960","Hs_planula_t.92541","Hs_planula_t.92543","Hs_planula_t.71193","Hs_planula_t.88622","Hs_planula_t.25801","Hs_planula_t.70828","Hs_planula_t.25244","Hs_planula_t.24836","Hs_planula_t.90763","Hs_planula_t.36207","Hs_planula_t.24723","Hs_planula_t.25097","Hs_planula_t.70841","Hs_planula_t.47554","Hs_planula_t.48675","Hs_planula_t.62380","Hs_planula_t.25179","Hs_planula_t.78332","Hs_planula_t.84960","Hs_planula_t.88759","Hs_planula_t.92543","Hs_planula_t.35971","Hs_planula_t.71823","Hs_planula_t.38528","Hs_planula_t.71193","Hs_planula_t.28270","Hs_planula_t.81533","Hs_planula_t.37383","Hs_planula_t.47734","Hs_planula_t.72997","Hs_planula_t.91957","Hs_planula_t.25801","Hs_planula_t.74420","Hs_planula_t.70828","Hs_planula_t.47720","Hs_planula_t.25244","Hs_planula_t.84026","Hs_planula_t.24836","Hs_planula_t.90763","Hs_planula_t.24676","Hs_planula_t.36207","Hs_planula_t.25097","Hs_planula_t.47554","Hs_planula_t.48675","Hs_planula_t.62380","Hs_planula_t.23286","Hs_planula_t.25179","Hs_planula_t.49134","Hs_planula_t.78332","Hs_planula_t.88759","Hs_planula_t.36723","Hs_planula_t.92543","Hs_planula_t.50251","Hs_planula_t.25278","Hs_planula_t.35971","Hs_planula_t.36173","Hs_planula_t.71823","Hs_planula_t.38528","Hs_planula_t.71193","Hs_planula_t.81533","Hs_planula_t.37383","Hs_planula_t.47734","Hs_planula_t.35740","Hs_planula_t.90667","Hs_planula_t.72997","Hs_planula_t.91957","Hs_planula_t.25801","Hs_planula_t.74420","Hs_planula_t.70828","Hs_planula_t.47720","Hs_planula_t.25244","Hs_planula_t.84026","Hs_planula_t.24836","Hs_planula_t.25243","Hs_planula_t.90763","Hs_planula_t.24676","Hs_planula_t.26325","Hs_planula_t.48076","Hs_planula_t.88569","Hs_planula_t.92541","Hs_planula_t.24676","Hs_planula_t.74636","Hs_planula_t.48878","Hs_planula_t.48076","Hs_planula_t.24676","Hs_planula_t.49187","Hs_planula_t.90763") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Hydractinia_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Hydractinia_seqid<- as.character(sig_degs_df$Hydractinia_seqid)
str(sig_degs_df)

homohydra.map<- read.table("Reduced_Hs_Sensory_Percep_Light_Stim_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("sens_percep_light_stim-all_symbols_for_OGs.txt")
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
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
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
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
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
homohydra.map<- read.table("Reduced_Hs_Sensory_Percep_Light_Stim_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 
length(unique(homohydra.map$symbol)) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("sens_percep_light_stim-all_symbols_for_OGs.txt") #Reduced
names(all_symbols_for_OGs)[2] <- c("tot_symbol")
str(all_symbols_for_OGs)
all_symbols_for_OGs$OG<- as.character(all_symbols_for_OGs$OG)
all_symbols_for_OGs$tot_symbol<- as.character(all_symbols_for_OGs$tot_symbol)
str(all_symbols_for_OGs)

nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:5],all_symbols_for_OGs,by="OG")
str(homohydra.map)
## make object of genes that are expressed but not sig de
not_sig_degs<-c("RDH10","CHRNB2","CRYGB","CNNM4","ARL6","MKKS","GABRR2","BBS1","SFRP5","CLDN19","HSF4","TIMP3","OAT","TH","OPA1","HPS1","RABGGTB","USH2A","KIFC3","DNAJC19","GJC1","CYP1B1","ZIC2","PDC","BBS5","UNC119","RPGR","RDH5","SIX3","AIPL1","DRAM2","OPA3","RP1","LUM","EYA3","RABGGTA","CNGA3","TULP1","CRYAA","SIX6","BBS2","RDH8","TULP2","RDH12","CRYGA","MYO5A","EYS","BBS4","RDH11","CRYGC","ABLIM1","PPT1","CNGB3","PDCL","GJD2","MYO7A","PDE6D","GLRA1","PITPNA","EYA4","WHRN","RP2","ATF6","SLC24A2","REEP6","CACNA2D4","NXNL2","RLBP1","BBS7","CRYZ")

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
#heatmap.2(TFcts,  scale="none", labRow = homohydra.map$symbol_seqid,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);


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

