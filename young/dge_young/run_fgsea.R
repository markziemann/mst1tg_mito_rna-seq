####################################################
# FGSEA
####################################################

dge$SYMBOL<-sapply((strsplit(rownames(dge)," ")),"[[",2)

res2 <- dge %>%
  dplyr::select(SYMBOL,stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(SYMBOL) %>%
  summarize(stat=mean(stat))
res2

stat<- deframe(res2)

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",
 destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")

gsets<-gmtPathways("ReactomePathways.gmt")

fgseaRes <- fgsea(pathways=gsets, stats=stat, minSize=10, nperm=10000)

fgseaRes<-fgseaRes[order(fgseaRes$pval),]


# peek at upregulated pathways
head(subset(fgseaRes,ES>0),10)

# peek at downregulated pathways
head(subset(fgseaRes,ES<0),10)

pdf("pathway_charts.pdf")

psig<-subset(fgseaRes,padj<=0.05)
plot(fgseaRes$ES,-log10(fgseaRes$pval),pch=19,cex=0.8,xlab="ES",ylab="-log10(p-value)")
points(psig$ES,-log10(psig$pval),pch=19,cex=0.8,xlab="ES",ylab="-log10(p-value)",col="red")
TOTAL=nrow(fgseaRes)
SIG=nrow(psig)
UP=length(which(psig$ES>0))
DN=length(which(psig$ES<0))
HEADER=paste(TOTAL,"gene sets examined,",SIG,"FDR<0.05,",UP,"up-regulated,",DN,"down-regulated")
mtext(HEADER)

# barchart top effect size FDR<0.05
par(mfrow = c(2, 1))
sig_up<-subset(psig,ES>0)
sig_dn<-subset(psig,ES<0)
sig_up<-head(sig_up[order(-sig_up$ES),],20)
sig_dn<-head(sig_dn[order(-sig_dn$ES),],20)

par( mar=c(5,28,4,2)) ; barplot(rev(sig_up$ES),horiz=T,names.arg=rev(sig_up$pathway),las=2, cex.names=0.6 , cex.axis=0.6, xlab="ES",main="up-regulated sets",cex.main=0.6)
par( mar=c(5,28,4,2)) ; barplot(rev(sig_dn$ES),horiz=T,names.arg=rev(sig_dn$pathway),las=2, cex.names=0.6, cex.axis=0.6, xlab="ES",main="down-regulated sets",cex.main=0.6)
dev.off()
# barchart top significant


# Pathway Table
psig$leadingEdge <- vapply(psig$leadingEdge, paste, collapse = ", ", character(1L))
write.table(psig, file = "Pathway_tables.csv", sep = ",")














