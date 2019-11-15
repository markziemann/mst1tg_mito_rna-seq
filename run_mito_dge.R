library("DESeq2")
library("mitch")


# obtain the count matrix
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE106201&format=file&file=GSE106201%5FMst1TgGal3KO%5Fcounts%2Etsv%2Egz",destfile="GSE106201_Mst1TgGal3KO_counts.tsv.gz")
x<-read.table("GSE106201_Mst1TgGal3KO_counts.tsv.gz",header=TRUE,row.names=1)
nTG_cols<-grep ( "wt_" , colnames(x) )
TG_cols<-grep ( "MstTg_" , colnames(x) )
x<-x[,c(nTG_cols,TG_cols)]
x<-x[which(rowMeans(x)>10),]
colnames(x)<-gsub("wt","nTG",colnames(x))
colnames(x)<-gsub("MstTg","TG",colnames(x))

# curate a sample sheet
ss<-data.frame(colnames(x))
rownames(ss)<-ss[,1]
ss$TG<-as.numeric(!grepl("nTG",rownames(ss)))
ss[,1]=NULL


# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = x,
                              colData = ss,
                              design= ~ TG)
dds <- DESeq(dds)
de <- results(dds)

# prepare for mitch analysis
genenames<-sapply( strsplit(rownames(de), "_" ),"[[",2)
gt<-data.frame(rownames(de),genenames)

# gene sets for enrichment analysis
genesets<-gmt_import("https://raw.githubusercontent.com/markziemann/mst1tg_mito_rna-seq/master/reactome.v5.2.symbols_mouse.gmt")
names(genesets)<-gsub("REACTOME_","",names(genesets))
names(genesets)<-gsub("_"," ",names(genesets))

# run mitch analysis
y<-mitch_import(as.data.frame(de),DEtype="deseq2",geneTable=gt)
res<-mitch_calc(y,genesets=genesets)
mitch_plots(res,outfile="plots.pdf")
mitch_report(res,outfile="plots.html")

