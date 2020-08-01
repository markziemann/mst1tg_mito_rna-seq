# to be run in R4
library("tidyverse")
library("reshape2")
library("DESeq2")
library("gplots")
library("fgsea")
library("mitch")


# read counts
tmp<-read.table("3col.tsv",header=F)
x<-as.matrix(acast(tmp, V2~V1, value.var="V3", fun.aggregate = sum))
x<-as.data.frame(x)
accession<-sapply((strsplit(rownames(x),"\\|")),"[[",2)
symbol<-sapply((strsplit(rownames(x),"\\|")),"[[",6)
x$geneid<-paste(accession,symbol)

xx<-aggregate(. ~ geneid,x,sum)
rownames(xx)<-xx$geneid
xx$geneid=NULL
xx<-round(xx)

colnames(xx) <- gsub(".fq-trimmed.fastq","",colnames(xx))

pdf("Mst1TG_old_plots.pdf")

# BARPLOT OF Mst1 expression (official gene name Stk4)
par(mar = c(10, 8, 8, 4) + 0.2)
DAT <- as.matrix(xx[grepl("Stk4$",rownames(xx)),])
NAME <- rownames(xx)[grep("Stk4$",rownames(xx))]
barplot(DAT,cex.names=0.9,main=NAME,las=2)

plot(cmdscale(dist(t(xx))), bty='l', xlab="Coordinate 1", ylab="Coordinate 2", 
 type = "p", pch=19, cex.axis=1.3, cex.lab=1.3 , col="gray") 
text(cmdscale(dist(t(xx))), labels=colnames(xx),cex=1.3) 


# curate the samplesheet
samplesheet<-as.data.frame(colnames(xx))
colnames(samplesheet) <- "samplename"
samplesheet <- data.frame(samplesheet[grep("Gal",samplesheet$samplename,invert=TRUE),])
colnames(samplesheet) <- "samplename"
samplesheet$tg <- as.numeric(grepl("Tg",samplesheet$samplename))
rownames(samplesheet)<-samplesheet[,1]
samplesheet[,1]=NULL

####################################################
# ctrl vs tg
####################################################
xx <- xx[,which(colnames(xx) %in% rownames(samplesheet) )]

y<-xx[which(rowSums(xx)/ncol(xx)>=(10)),]

dds <- DESeqDataSetFromMatrix(countData = y , colData = samplesheet, design = ~ tg )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(as.data.frame(z),assay(vsd))
dge<-as.data.frame(zz[order(zz$pvalue),])
write.table(dge,file="Mst1TG_old_deseq.tsv",quote=F,sep="\t")

#some plots
sig<-subset(dge,padj<0.05)
SIG=nrow(sig)
DN=nrow(subset(sig,log2FoldChange<0))
UP=nrow(subset(sig,log2FoldChange>0))
HEADER=paste("nTG vs Mst1-TG:", SIG , "DGEs,", UP ,"upregulated,", DN, "downregulated")

plot(log2(dge$baseMean),dge$log2FoldChange,cex=0.6,cex.axis=1.2,cex.lab=1.3, 
 xlab="log2 base mean",
 ,ylab="log2 fold change" ,pch=19,col="#838383")

points(log2(sig$baseMean),sig$log2FoldChange,cex=0.6,pch=19,col="red")

mtext((HEADER),cex=1.2)

top<-head(sig,20)
#text(log2(top$baseMean)+1, top$log2FoldChange, labels = rownames(top),cex=0.7)

#volcano plot

plot(dge$log2FoldChange, -log2(dge$pvalue)+1E-307 ,cex=0.6, cex.lab=1.3,cex.axis=1.2,
 xlim=c(-3,3),xlab="log2 fold change", ylab="-log2 p-value" ,pch=19,col="#838383")

points(sig$log2FoldChange, -log2(sig$pvalue)+1E-307, cex=0.6,pch=19,col="red")	
mtext((HEADER),cex=1.2)

# top N gene heatmap
colfunc <- colorRampPalette(c("blue", "white", "red"))

heatmap.2(  as.matrix(dge[1:50,c(7:ncol(dge))]), col=colfunc(25),scale="row",
 trace="none",margins = c(6,20), cexRow=.6, cexCol=.8,  main="Top 50 genes")

dev.off()

####################################################
# homology mapping
####################################################
# ensembl gene 99 from biomart
gt <- read.csv("../ref/mart_export.txt",sep="\t")
gt <- gt[,c(1,3)] 

####################################################
# gene sets
####################################################
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",
 destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")

####################################################
# mitch
####################################################
gsets<-gmt_import("ReactomePathways.gmt")
dge2 <- dge
rownames(dge2) <- sapply(strsplit(rownames(dge),"\\."),"[[",1)
yy<-mitch_import(dge2, DEtype="deseq2", geneTable=gt)
res <- mitch_calc(yy, gsets)
unlink("Mst1TG_old_mitch.html")
mitch_report(res, outfile="Mst1TG_old_mitch.html")
mitch_plots(res, outfile="Mst1TG_old_mitch.pdf")

####################################################
# mitch with mouse gene sets
###################################################
gsets<-gmt_import("ReactomePathways_mouse.gmt")
gt <- read.csv("../ref/mart_export.txt",sep="\t")
gt <- gt[,c(1,4)]
yy<-mitch_import(dge2, DEtype="deseq2", geneTable=gt)
res <- mitch_calc(yy, gsets)
unlink("Mst1TG_old_mitch_mouse.html")
mitch_report(res, outfile="Mst1TG_old_mitch_mouse.html")
mitch_plots(res, outfile="Mst1TG_old_mitch_mouse.pdf")

####################################################
# heatmaps of top sets
####################################################
topsets <- head(res$enrichment_result,50)$set
mytopsets <- gsets[which(names( gsets ) %in% topsets)]

i=1

pdf("top_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(mytopsets) ) {

  g <- as.vector(unlist(mytopsets[i]))
  symbol <- sapply((strsplit(rownames(dge)," ")),"[[",2)
  z <- dge[which(symbol %in% g),7:ncol(dge)]
  SETNAME=names(mytopsets[i])
  heatmap.2( as.matrix(z), col=colfunc(25),scale="row",
   trace="none",margins = c(6,20), cexRow=.4, cexCol=.8)
  mtext(SETNAME)

}

dev.off()

####################################################
# mitch with gene sets of interest
###################################################
gsets<-gmt_import("pathways_of_interest_mouse.gmt")
gt <- read.csv("../ref/mart_export.txt",sep="\t")
gt <- gt[,c(1,4)]
yy<-mitch_import(dge2, DEtype="deseq2", geneTable=gt)
res <- mitch_calc(yy, gsets)
unlink("Mst1TG_old_mitch_mouse_specificsets.html")
mitch_report(res, outfile="Mst1TG_old_mitch_mouse_specificsets.html")
mitch_plots(res, outfile="Mst1TG_old_mitch_mouse_specificsets.pdf")

####################################################
# heatmaps of pathways of interest
####################################################

gsets<-gmt_import("pathways_of_interest_mouse.gmt")

pdf("pathways_of_interest_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(gsets) ) {

  g <- as.vector(unlist(gsets[i]))
  symbol <- sapply((strsplit(rownames(dge)," ")),"[[",2)
  z <- dge[which(symbol %in% g), 7:ncol(dge)]
  SETNAME=names(gsets[i])

  if ( nrow(z) > 2 ) {
    heatmap.2( as.matrix(z), col=colfunc(25),scale="row",
     trace="none",margins = c(6,20), cexRow=.4, cexCol=.8)
    mtext(SETNAME)
  }

}

dev.off()


####################################################
# heatmap of mitochondrial
####################################################

gsets<-gmt_import("../ref/mouse_mt.gmt")
pdf("mtdna_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
i=1
g <- as.vector(unlist(gsets[i]))
symbol <- sapply((strsplit(rownames(dge)," ")),"[[",2)
z <- dge[which(symbol %in% g), 7:ncol(dge)]
SETNAME=names(gsets[i])

if ( nrow(z) > 2 ) {
  heatmap.2( as.matrix(z), col=colfunc(25),scale="row",
   trace="none",margins = c(6,20), cexRow=.8, cexCol=.8)
  mtext(SETNAME)
}

dev.off()

save.image("run_dgepw.RData")

