##EdgeR - Hydractinia data - 4 stages/days 4-16-20
  #3 reps of 4 groups which are days

setwd("~/Desktop/R/Hydractinia_transcriptomics_3-19-20")   #<-----set me!

#source("https://bioconductor.org/biocLite.R"); biocLite("tximport"); install.packages("readr"); biocLite("edgeR")    #<-----uncomment & run if any not installed yet

library(tximport); library(readr); library(edgeR); library(limma)

#in your cwd make a text file named libraries_to_stages.txt - this should have 2 columns - the first is the name
#of your sample and the second (separated by a space) is the group number it belongs to (so the stage)

#Make a dir in your cwd named mapping - place all of your output dirs for each sample that was outputted by salmon in this dir
dir <- getwd()
list.files()

#library info file
devstages<-read.table("libraries_to_stages.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4


##pick colors for each library type
stage1col<-"#DE3C00"
stage2col <-"#DE7F00"
stage3col <-"#08548F"
stage4col <-"aquamarine3"


################################
### Multidimensional scaling ###
################################

runEdgeRmds<-function(salmonquantfiles, groups, colors){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
  
  keep <- rowSums(cpm(cts)>10.0) >= length(groups)
  cts<- cts[keep,]
  dim(cts)
  print(colSums(cts))
  
  group <- groups
  y <- DGEList(counts=cts ,group=group)
  y <- calcNormFactors(y)
  y <-estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y, prior.n=6.25) # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(12-4) = 50/8 =6.25
  
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y),col=colors)
  legend("bottomright",legend=paste("Stage ",levels(as.factor(groups))),text.col=colors[seq(1,length(colors),4)]) #change the last number to how many groups you have
  
}
#pdf(file="edgeRplots.pdf",width=8, height=11)
cols<-c(stage1col, stage2col, stage3col, stage4col)
runEdgeRmds(c(dev1files,dev2files, dev3files, dev4files), c(1,1,1,2,2,2,3,3,3,4,4,4),rep(cols,each=3)) #change the last number to how many groups you have




##############################################
### pairwise comparisons between libraries ###
##############################################

runEdgeRpwise<-function(salmonquantfiles, groups,colors){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
  
  keep <- rowSums(cpm(cts)>10.0) >= length(groups)
  cts<- cts[keep,]
  dim(cts)
  print(colSums(cts))
  
  y <- DGEList(counts=cts ,group=groups)
  y <- calcNormFactors(y)
  
  y<-estimateCommonDisp(y)
  # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(36-3) = 50/33 =1.515152
  y <- estimateTagwiseDisp(y, prior.n=6.25)
  
  group<-levels(as.factor(groups))		
  et<-exactTest(y, pair=c(group[1],group[2]))
  tab<-summary(de <- decideTestsDGE(et, p=0.0001, adjust="BH"))
  n<-tab[1]+tab[3]
  detags <- rownames(y)[as.logical(de)]
  
  
  plotSmear(et, de.tags=detags, main="DGE Exact Test")
  abline(h = c(-2, 2), col = "blue")
  abline(h = c(-4, 4), col = "blue")
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y), col=colors)
  textcol<-colors[seq(1,length(colors),3)]
  legend("center",legend=paste("Stage ",group),text.col= textcol,bty="n")
  
  plotBCV(y, main="BCV plot")
  
  meanVarPlot <- plotMeanVar(estimateCommonDisp(y) , show.raw.vars=TRUE,
                             show.tagwise.vars=TRUE,
                             show.binned.common.disp.vars=FALSE,
                             show.ave.raw.vars=FALSE , NBline=TRUE,
                             nbins=100,
                             pch=16,
                             xlab="Mean Expresion (Log10)",
                             ylab="Variance (Log10)",
                             main="Mean-Variance Plot")
  
  #positive FC: higher in group2 than 1 
  return(topTags(et, n=n))
}


#evaluate stage1 vs stage2 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s1vs2<-runEdgeRpwise(c(dev1files,dev2files), c(1,1,1,2,2,2),rep(c(stage1col, stage2col),each=3))
title("Stage 1 vs 2",outer=T)
write.table(s1vs2, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s1vs2_Hs_updated.txt", sep="\t")

#evaluate stage1 vs stage3
par(mfrow=c(2,2),oma=c(1,1,2,0))
s1vs3<-runEdgeRpwise(c(dev1files,dev3files), c(1,1,1,3,3,3), rep(c(stage1col, stage3col),each=3))
title("Stage 1 vs 3",outer=T)
write.table(s1vs3, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s1vs3_Hs_updated.txt", sep="\t")

#evaluate stage1 vs stage4
par(mfrow=c(2,2),oma=c(1,1,2,0))
s1vs4<-runEdgeRpwise(c(dev1files,dev4files), c(1,1,1,4,4,4), rep(c(stage1col, stage4col),each=3))
title("Stage 1 vs 4",outer=T)
write.table(s1vs4, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s1vs4_Hs_updated.txt", sep="\t")

#evaluate stage2 vs stage3
par(mfrow=c(2,2),oma=c(1,1,2,0))
s2vs3<-runEdgeRpwise(c(dev2files,dev3files), c(2,2,2,3,3,3),rep(c(stage2col, stage3col),each=3))
title("Stage 2 vs 3",outer=T)
write.table(s2vs3, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s2vs3_Hs_updated.txt", sep="\t")

#evaluate stage2 vs stage4
par(mfrow=c(2,2),oma=c(1,1,2,0))
s2vs4<-runEdgeRpwise(c(dev2files,dev4files), c(2,2,2,4,4,4),rep(c(stage2col, stage4col),each=3))
title("Stage 2 vs 4",outer=T)
write.table(s2vs4, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s2vs4_Hs_updated.txt", sep="\t")

#evaluate stage3 vs stage4
par(mfrow=c(2,2),oma=c(1,1,2,0))
s3vs4<-runEdgeRpwise(c(dev3files,dev4files), c(3,3,3,4,4,4), rep(c(stage3col, stage4col),each=3))
title("Stage 3 vs 4",outer=T)
write.table(s3vs4, "C:\\Desktop\\R\\Transcriptomics_Actinula_3-17-20\\s3vs4_Hs_updated.txt", sep="\t")

options(max.print = 1000000)

#these objects list the de genes. Note the logFC (log fold change) column. 
#You can then find these seqs in the assemblies
s1vs2
s2vs3
par(mfrow=c(1,1))
barplot(c(nrow(s1vs2), nrow(s1vs3), nrow(s1vs4), nrow(s2vs3), nrow(s2vs4), nrow(s3vs4)),names.arg=c("s1vs2","s1vs3","s1vs4","s2vs3","s2vs4","s3vs4"), main="Number of differentially expressed transcripts")

#dev.off()
