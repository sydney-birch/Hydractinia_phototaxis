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


library(tximport); library(readr); library(edgeR)

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
homohydra.map<- read.table("Reduced_Hs_TF_activity_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/tf_activity-all_symbols_for_OGs.txt")
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
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(3, 10),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
#heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);

#heatmaps specifiying color break
#heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="none", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette2, breaks=col_breaks, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);



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
sig_degs<-c("Hs_planula_t.47705","Hs_planula_t.70805","Hs_planula_t.48339","Hs_planula_t.25178","Hs_planula_t.27430","Hs_planula_t.69265","Hs_planula_t.48371","Hs_planula_t.78088","Hs_planula_t.49923","Hs_planula_t.47814","Hs_planula_t.80358","Hs_planula_t.50013","Hs_planula_t.78195","Hs_planula_t.35505","Hs_planula_t.84109","Hs_planula_t.29942","Hs_planula_t.47705","Hs_planula_t.81533","Hs_planula_t.80291","Hs_planula_t.92607","Hs_planula_t.26803","Hs_planula_t.70805","Hs_planula_t.25178","Hs_planula_t.25953","Hs_planula_t.25766","Hs_planula_t.27430","Hs_planula_t.47531","Hs_planula_t.69265","Hs_planula_t.48371","Hs_planula_t.49854","Hs_planula_t.26576","Hs_planula_t.29026","Hs_planula_t.39746","Hs_planula_t.92085","Hs_planula_t.48444","Hs_planula_t.30796","Hs_planula_t.78088","Hs_planula_t.40269","Hs_planula_t.92611","Hs_planula_t.28487","Hs_planula_t.82707","Hs_planula_t.25996","Hs_planula_t.49385","Hs_planula_t.36887","Hs_planula_t.59560","Hs_planula_t.75906","Hs_planula_t.48810","Hs_planula_t.82724","Hs_planula_t.28669","Hs_planula_t.49923","Hs_planula_t.91872","Hs_planula_t.78510","Hs_planula_t.47814","Hs_planula_t.35632","Hs_planula_t.89953","Hs_planula_t.26313","Hs_planula_t.80358","Hs_planula_t.73874","Hs_planula_t.38966","Hs_planula_t.50013","Hs_planula_t.26485","Hs_planula_t.78195","Hs_planula_t.35505","Hs_planula_t.29942","Hs_planula_t.90558","Hs_planula_t.28219","Hs_planula_t.27398","Hs_planula_t.47705","Hs_planula_t.81533","Hs_planula_t.80291","Hs_planula_t.92607","Hs_planula_t.70805","Hs_planula_t.88704","Hs_planula_t.36186","Hs_planula_t.24738","Hs_planula_t.49951","Hs_planula_t.37204","Hs_planula_t.24909","Hs_planula_t.47531","Hs_planula_t.69265","Hs_planula_t.48371","Hs_planula_t.91026","Hs_planula_t.42186","Hs_planula_t.26576","Hs_planula_t.39686","Hs_planula_t.24927","Hs_planula_t.91768","Hs_planula_t.48444","Hs_planula_t.78088","Hs_planula_t.40269","Hs_planula_t.82707","Hs_planula_t.25889","Hs_planula_t.75906","Hs_planula_t.66399","Hs_planula_t.48810","Hs_planula_t.28669","Hs_planula_t.91872","Hs_planula_t.78510","Hs_planula_t.47814","Hs_planula_t.35632","Hs_planula_t.89953","Hs_planula_t.26313","Hs_planula_t.80358","Hs_planula_t.73874","Hs_planula_t.35843","Hs_planula_t.48240","Hs_planula_t.74424","Hs_planula_t.48593","Hs_planula_t.48386","Hs_planula_t.83583","Hs_planula_t.38966","Hs_planula_t.50013","Hs_planula_t.26485","Hs_planula_t.78195","Hs_planula_t.24738","Hs_planula_t.91768","Hs_planula_t.40269","Hs_planula_t.36887","Hs_planula_t.75211","Hs_planula_t.81810","Hs_planula_t.35632","Hs_planula_t.78195","Hs_planula_t.24738","Hs_planula_t.47531","Hs_planula_t.69264","Hs_planula_t.91768","Hs_planula_t.48444","Hs_planula_t.35632","Hs_planula_t.78195","Hs_planula_t.47531","Hs_planula_t.39746","Hs_planula_t.48275","Hs_planula_t.30796","Hs_planula_t.36887","Hs_planula_t.38820","Hs_planula_t.25889","Hs_planula_t.81810","Hs_planula_t.24627","Hs_planula_t.47814","Hs_planula_t.83583","Hs_planula_t.78195"
) #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Hydractinia_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Hydractinia_seqid<- as.character(sig_degs_df$Hydractinia_seqid)
str(sig_degs_df)

homohydra.map<- read.table("Reduced_Hs_TF_activity_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("~/Desktop/R /Hydractinia_transcriptomics_3-19-20/DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/tf_activity-all_symbols_for_OGs.txt")
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
heatmap.2(t(dev.means2), labRow = homohydra.map$symbol_seqid, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL,  margins=c(2, 14),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
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
homohydra.map<- read.table("Reduced_Hs_TF_activity_R_output_genes.txt",header=T, stringsAsFactors=F) #Reduced

str(homohydra.map) #should be data frame
length(homohydra.map) 
length(unique(homohydra.map$symbol)) 

#import collapsed symbols for OG - add to dataset
all_symbols_for_OGs <- read.delim("DEGs 4 groups 0.0001/Updated_Reduced_sig_DEG_Summaries_3-26-21/tf_activity-all_symbols_for_OGs.txt")
names(all_symbols_for_OGs)[2] <- c("tot_symbol")
str(all_symbols_for_OGs)
all_symbols_for_OGs$OG<- as.character(all_symbols_for_OGs$OG)
all_symbols_for_OGs$tot_symbol<- as.character(all_symbols_for_OGs$tot_symbol)
str(all_symbols_for_OGs)

nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:5],all_symbols_for_OGs,by="OG")
str(homohydra.map)
## make object of genes that are expressed but not sig de
## there were two many items in my string so I had to break up the string and merge into one data frame
not_sig_degs <- c("TWIST1","OLIG2","TSHZ1","OTX1","EN1","ZRANB2","IRF1","MAFG","NFAT5","DLX3","POU4F2","TP53","POU2F3","BARX1","NR2F2","HES2","NOC3L","ZIC4","STAT5A","GRHL1","CXXC1","PCBP3","HOXB2","IRX1","FOXA2","SOX7","TBX4","SUPT6H","DBX1","LITAF","FOXO4","NKX2-3","ZNF639","MLXIP","ZBED5","VAX2","BCL11B","CEBPG","MEIS2","THAP12","PAX8","SIM2","E2F6","ZBED1","THAP11","IRF2BP2","ZNF593","MIER3","TLX3","IRF2BPL","FOXJ1","SOX1","DLX5","MLX","OVOL3","LYL1","TERF1","RUNX1T1","GSC2","HIRA","OSR1","ZNF197","TP73","TLX2","FOXI2","ARID3C","TCF15","SOX10","ZZZ3","SP7","ATOH8","MTA3","FIGLA","GATA2","TRPS1","MECOM","FOXO6","IRX2","FOXS1","RELB","PRDM1","POGZ","MAEL","SIRT1","ZNF770","PRDM11","PRDM6","NFXL1","HMX1","YEATS2","NKX1-1","SMAD2","TBX20","ZBTB41","BLZF1","NKX6-2","NFIA","UNKL","SOX3","GRHL3","ARID5B","HSF2","CREM","ARID2","ISL1","HMX2","HLTF","GLMP","PITX2","ARID4A","ZBTB47","MITF","GLIS3","ZIC3","KDM5A","NR2F6","PHB","LHX5","STAT6","BARHL2","SOX14","TCF21","YBX1","BHLHE23","TBX10","KDM5D","SMAD1","PRDM14","EAF2","ESRRA","IRF8","IRF2BP1","NHLH2","SOX9","MXD1","EZH2","TBR1","PROP1","SHOX2","RAX2","CBFA2T2","CC2D1A","EN2","KLF14","NEUROG3","ZIC5","FUBP1","SRY","SP8","MSC","CRX","ZNF395","KDM3A","ATOH1","SMAD9","FOXD4L5","ZNF28","ERG","ZNF518A","TAF1","BHLHA9","MTF1","PURG","PFDN1","RXRB","NPAS3","CAMTA2","CARF","HOXB5","LMX1B","ARNT","HOXD4","KLF5","NKX3-1","USF3","TEF","NEUROG1","ZBED6","KLF10","YBX2","MSX2","PAXBP1","TAL2","ZIC2","NFE2L3","FUBP3","PHF5A","ASCL2","MXI1","TFAP4","THAP1","ATRX","ARID1B","FOXK1","NKX2-8","FOXA3","FOXK2","CREB3L1","MYBL2","ZC3H15","AFF1","FOXL2","SIX3","HOXA4","JARID2","ATOH7","FOXD4","ESRRG","LMO4","FOXC2","SP9","KDM5C","ARID1A","YBX3","VSX1","ZBTB18","CDX2","YY1","FOXL1","MLLT3","MAFF","SIX6","AEBP2","LZTR1","LHX3","LHX4","KLF1","PLAGL1","SLC2A4RG","NONO","TXK","YY2","JUN","SOX2","FOXN2","OSR2","CRAMP1","GSX1","FOXF1","PHOX2A","PAX3","E2F3","IRF4","CBFA2T3","DR1","HSF1","ZFP42","SOX17","LMO2","PLSCR1","DLX6","SNAPC4","NFE2L2","FERD3L","LMX1A","E2F1","ATF1","TCF25","RXRA","KLF4","BARHL1","CASZ1","NFKB2","TADA2B","ASCL3","GABPB1","POU3F2","AHR","MEIS3","ZNF367","BHLHE22","RERE","CREB5","NKX1-2","STK16","BRD8","IRF2","BOLA3","WNT5A","MXD4","HES1","ZBED4","CREB1","TSHZ2","TAF13","GBX1","TBX6","DPF2","E2F5","PSPC1","SMAD5","NFX1","FOXN3","TLE4","HOXC5","NFIC","NFYB","NKX2-2","TFE3","POU3F4","BARX2","TFCP2","ZNF385A","PAX5","ESRRB","ZBTB37","IRX3","NFE2","ZNFX1","TFEC","MAFA","ISL2","HSF4","SOX18","YEATS4","ZNF280B","MNT","SOHLH2","THAP2","DBX2","TAF1L","OLIG1","TEAD3","PITX3","IRX5","ARID4B","UBP1","PRDM9","POU3F3","FOXD4L4","PHTF1","ATF7","TP63","LHX1","FOXD2","SP5","ZNF280A","ZNF622","PRRX1","HMX3","EGR1","ZBTB8B","BMPR1A","NSD2","KLF7","TEAD4","POU4F3","SCX","PRDM16","PLAG1","SHOX","ETV3","TARDBP","ARID3B","NEUROD1","HOXA5","TFEB","MIER2","ZNF516","AHCTF1","GLIS2","ASCL5","SFPQ","ZNF579","NFKB1","AFF4","PBX1","CEBPE","FOXD4L1","JAZF1","POU3F1","BOLA2","KLF3","GSX2","RCOR2","RCOR1","KDM2A","ZC3H6","HES5","HHEX","RFX6","HNF4A","NHLH1","PITX1","PRDM7","CIC","NOLC1","HOXB1","DMAP1","IRX4","SPDEF","VSX2","SIX2","CREBZF","SOX8","MXD3","ATF2","PROX1","LHX6","ASCL1","NKX6-3","PLEK","ZNF280C","HNRNPK","NFYC","FOXB2","TCF7","SNAI1","KLF2","CAMTA1","LHX8","NCOA1","UNCX","ENO1","JUNB","FOXA1","BSX","PTF1A","E2F4","FOXD1","TOX4","SMAD3","RFX7","ZBTB4","OVOL1","NRL","DLX4","ZNF292","SIM1","HES4","ZFPM1","CEBPB","RCOR3","TBX2","PKNOX1","OTP","OTX2","MTA2","MIER1","PRDM5","THAP9","PAX1","EZH1","GATA1","TAF7","KDM5B","SLC30A9","TEAD1","SMAD4","HOXC4","GATAD1","HMG20A","HOXB4","HLF","CREB3","MEOX2","STAT5B","NFIX","IRX6","BHLHA15","TBX5","NOC4L","PAX9","TBX1","POU1F1","BOLA1","PRDM12","NKX2-6","KLF13","HOXA1","BTAF1","HMG20B","HIF3A","NFYA","ARID3A","EMX2","MLLT1","KLF6","CREB3L4","NEUROG2","BTG2","AFF3","TAL1","FLI1","ARNT2","ATF6B","ASCL4","PCGF5","NKX6-1","NR2F1","ATF6","SOX21","SCML2","FEV","ARHGAP35","MYCN","FOXB1","SRF","HIF1A","RXRG","TAF6","ELF1","DBP","SUPT4H1","FOXE1","MAFK","ZIC1","GATA3","RAI1","JUND","ZNF277","ZBTB24","HDAC1","PCBP1","NKRF","ZNF511","MSX1")
not_sig_degs_df<-data.frame(not_sig_degs)
names(not_sig_degs_df)[1] <- c("symbol")

not_sig_degs2<- c("POU5F1B","PRRX2","CREB3L3","BCL11A","ELF2","MTA1","PKNOX2","MYBL1","BOLA2B","RCAN1","GATA4","OVOL2","PAX2","PCGF3","PHOX2B","SP6","SUB1","POU4F1","HLX","THAP3","RUNX1","PLAGL2","NKX2-1","E2F2","PRDM15","AHRR","SIX1","GABPA","GRHL2","FOXD4L6","RFX5","TFCP2L1","NKX2-4","FOXD3","ZNF518B","PRDM4","NPAS1","ZGPAT","FOXD4L3","LEF1","KLF12","MEIS1","EGR3","CARHSP1","CIR1","EMX1","OLIG3","ZC3H3","KAT7","HOXA2","TLX1","SOX15","DNMT3B","CBFB","TAF5","GON4L","FOXC1","PAX7","MLLT10","FOXE3","DRGX","HNF4G","TEAD2","FOXG1","HINFP","MLXIPL","TSHZ3","RBPJ","DMTF1","MAFB","RELA","REL")
not_sig_degs_df2<-data.frame(not_sig_degs2)
names(not_sig_degs_df2)[1] <- c("symbol")

#make one big data frame - bind two data frames into one
not_sig_degs_df_tot<- rbind(not_sig_degs_df, not_sig_degs_df2) 
length(not_sig_degs_df_tot$symbol) 

#join with homohydra map to get genes that are expressed but not sig DE expressed
nrow(homohydra.map)
homohydra.map<-inner_join(homohydra.map[,1:6],not_sig_degs_df_tot,by="symbol")
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

homohydra.map<- homohydra.map[!duplicated(homohydra.map$Hydractinia_seqid), ]
length(homohydra.map$Hydractinia_seqid) 

keep<-unique(homohydra.map$Hydractinia_seqid); length(keep) #1 col with 32 enteries 
length(unique(homohydra.map$tot_symbol)) 
length(unique(homohydra.map$symbol)) 

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) 
sort(homohydra.map$symbol)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
options(max.print=999999)
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



dev.off()



#layout key (lmat)
# 1-Heatmap,
# 2-Row dendrogram,
# 3-Column dendrogram,
# 4-Key

