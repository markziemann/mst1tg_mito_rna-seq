library("mitch")

# the purpose of this script is to create a volcano plot of yap targets

tf <- gmt_import("../ref/TFofinterest.gmt")

YapTargets <- tf[[grep("48242",names(tf))]]

old <- read.table("Mst1TG_old_deseq.tsv",header=TRUE, row.names=1,sep="\t")

young <- read.table("Mst1TG_young_deseq.tsv",header=TRUE, row.names=1,sep="\t")

old_yap <-  old[which(sapply(strsplit(rownames(old)," "),"[[",2) %in% YapTargets),]

young_yap <- young[which(sapply(strsplit(rownames(young)," "),"[[",2) %in% YapTargets),]

# volcano
make_volcano <- function(de,name) {
    sig <- subset(de,padj<0.05)
    N_SIG=nrow(sig)
    N_UP=nrow(subset(sig,log2FoldChange>0))
    N_DN=nrow(subset(sig,log2FoldChange<0))
    HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn")
    plot(de$log2FoldChange,-log10(de$padj),cex=0.5,pch=19,col="darkgray",
        main=name, xlab="log2 FC", ylab="-log10 pval", xlim=c(-6,6))
    mtext(HEADER)
    grid()
    points(sig$log2FoldChange,-log10(sig$padj),cex=0.5,pch=19,col="red")
}

pdf("charts_2020-08-01.pdf")
make_volcano(old,"15 wk (all)")
make_volcano(old_yap,"15 wk (Yap targets)")

make_volcano(young,"3 wk (all)")
make_volcano(young_yap,"3 wk (Yap targets)")

## Venn diagram
old_up <- rownames(subset(old, padj<0.05 & log2FoldChange>0))
old_dn <- rownames(subset(old, padj<0.05 & log2FoldChange<0))

young_up <- rownames(subset(young, padj<0.05 & log2FoldChange>0))
young_dn <- rownames(subset(young, padj<0.05 & log2FoldChange<0))

old_yap_up <- rownames(subset(old_yap, padj<0.05 & log2FoldChange>0)) 
old_yap_dn <- rownames(subset(old_yap, padj<0.05 & log2FoldChange<0))

young_yap_up <- rownames(subset(young_yap, padj<0.05 & log2FoldChange>0))
young_yap_dn <- rownames(subset(young_yap, padj<0.05 & log2FoldChange<0))




library("eulerr")
v1 <- list("up 3 wk" = young_up , 
           "dn 3 wk" = young_dn ,
           "up 15 wk" = old_up ,
           "dn 15 wk" = old_dn )
plot(euler(v1, shape = "ellipse"), quantities = TRUE, main = "All dysregulated genes")



v1 <- list("up 3 wk" = young_yap_up ,
           "dn 3 wk" = young_yap_dn ,
           "up 15 wk" = old_yap_up ,
           "dn 15 wk" = old_yap_dn )
plot(euler(v1, shape = "ellipse"), quantities = TRUE, main = "Yap target genes")

dev.off()
