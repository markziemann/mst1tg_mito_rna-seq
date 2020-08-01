# to be run in R4
library("tidyverse")
library("reshape2")
library("gplots")
library("mitch")

####################################################
# import profiling data
####################################################
young <- read.csv("Mst1TG_young_deseq.tsv",sep="\t")
rownames(young) <- sapply((strsplit(rownames(young),"\\.")),"[[",1)
old <- read.csv("Mst1TG_old_deseq.tsv",sep="\t")
rownames(old) <- sapply((strsplit(rownames(old),"\\.")),"[[",1)
x <- list("young" = young , "old" = old)

# import these in again for heatmaps
young <- read.csv("Mst1TG_young_deseq.tsv",sep="\t")
old <- read.csv("Mst1TG_old_deseq.tsv",sep="\t")

####################################################
# gene sets
####################################################
gsets <- gmt_import("../ref/ReactomePathways_mouse.gmt")

####################################################
# gene accession to symbol mapping
####################################################
gt <- read.csv("../ref/mart_export.txt",sep="\t")
gt <- gt[,c(1,4)]

####################################################
# mitch priority significance
####################################################
yy <- mitch_import(x, DEtype="deseq2", geneTable=gt)
res <- mitch_calc(yy, gsets, priority="significance")
unlink("Mst1TG_young_and_old_sig_mitch.html")
# top N gene heatmap
mitch_report(res, outfile = "Mst1TG_young_and_old_sig_mitch.html")
mitch_plots(res, outfile = "Mst1TG_young_and_old_sig_mitch.pdf")

####################################################
# heatmaps of top sets by significance
####################################################
# colour scheme
colfunc <- colorRampPalette(c("blue", "white", "red"))

topsets <- head(res$enrichment_result,50)$set
mytopsets <- gsets[which(names( gsets ) %in% topsets)]
yo <- young[,7:ncol(young)]
colnames(yo) <- paste("y",colnames(yo))
ol <- old[,7:ncol(old)]
colnames(ol) <- paste("o",colnames(ol))
mx <- merge(yo,ol,by=0)
rownames(mx) <- mx$Row.names
mx$Row.names = NULL

i = 1
pdf("top_sig_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(mytopsets) ) {
  g <- as.vector(unlist(mytopsets[i]))
  symbol <- sapply((strsplit(rownames(mx)," ")),"[[",2)
  z <- mx[which(symbol %in% g),]
  SETNAME=names(mytopsets[i])
  heatmap.2( as.matrix(z), col=colfunc(25),scale="row",
   trace="none",margins = c(6,20), cexRow=.4, cexCol=.8)
  mtext(SETNAME)
}
dev.off()

####################################################
# mitch priority effect size
####################################################
res <- mitch_calc(yy, gsets, priority = "effect")
unlink("Mst1TG_young_and_old_eff_mitch.html")
mitch_report(res, outfile = "Mst1TG_young_and_old_eff_mitch.html")
mitch_plots(res, outfile = "Mst1TG_young_and_old_eff_mitch.pdf")

####################################################
# heatmaps of top sets by effect size
####################################################
topsets <- head(res$enrichment_result,50)$set
mytopsets <- gsets[which(names( gsets ) %in% topsets)]
i=1
pdf("top_eff_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(mytopsets) ) {
  g <- as.vector(unlist(mytopsets[i]))
  symbol <- sapply((strsplit(rownames(mx)," ")),"[[",2)
  z <- mx[which(symbol %in% g),]
  SETNAME = names(mytopsets[i])
  heatmap.2( as.matrix(z), col = colfunc(25),scale = "row",
   trace="none",margins = c(6,20), cexRow = .4, cexCol = .8)
  mtext(SETNAME)
}
dev.off()


####################################################
# mitch priority discordant effect size
####################################################
res <- mitch_calc(yy, gsets, priority = "SD")
unlink("Mst1TG_young_and_old_dis_mitch.html")
mitch_report(res, outfile = "Mst1TG_young_and_old_dis_mitch.html")
mitch_plots(res, outfile = "Mst1TG_young_and_old_dis_mitch.pdf")

####################################################
# heatmaps of top sets by discordant size
####################################################
topsets <- head(res$enrichment_result,50)$set
mytopsets <- gsets[which(names( gsets ) %in% topsets)]
i=1
pdf("top_dis_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(mytopsets) ) {
  g <- as.vector(unlist(mytopsets[i]))
  symbol <- sapply((strsplit(rownames(mx)," ")),"[[",2)
  z <- mx[which(symbol %in% g),]
  SETNAME = names(mytopsets[i])
  heatmap.2( as.matrix(z), col = colfunc(25),scale = "row",
   trace="none",margins = c(6,20), cexRow = .4, cexCol = .8)
  mtext(SETNAME)
}
dev.off()


####################################################
# gene sets of interest
####################################################
gsets <- gmt_import("../ref/pathways_of_interest_mouse.gmt")

####################################################
# mitch pathways of interest
####################################################
yy <- mitch_import(x, DEtype="deseq2", geneTable=gt)
res <- mitch_calc(yy, gsets, priority="significance")
unlink("Mst1TG_young_and_old_pathwayofinterest_mitch.html")
# top N gene heatmap
mitch_report(res, outfile = "Mst1TG_young_and_old_pathwayofinterest_mitch.html")
mitch_plots(res, outfile = "Mst1TG_young_and_old_pathwayofinterest_mitch.pdf")

####################################################
# heatmaps of pathways of interest
####################################################
# colour scheme
colfunc <- colorRampPalette(c("blue", "white", "red"))

topsets <- head(res$enrichment_result,50)$set
mytopsets <- gsets[which(names( gsets ) %in% topsets)]
yo <- young[,7:ncol(young)]
colnames(yo) <- paste("y",colnames(yo))
ol <- old[,7:ncol(old)]
colnames(ol) <- paste("o",colnames(ol))
mx <- merge(yo,ol,by=0)
rownames(mx) <- mx$Row.names
mx$Row.names = NULL

i = 1
pdf("pathwayofinterest_heatmaps.pdf")
par(mar = c(10, 8, 8, 4) + 0.2)
for ( i in seq_along(mytopsets) ) {
  g <- as.vector(unlist(mytopsets[i]))
  symbol <- sapply((strsplit(rownames(mx)," ")),"[[",2)
  z <- mx[which(symbol %in% g),]
  SETNAME=names(mytopsets[i])
  heatmap.2( as.matrix(z), col=colfunc(25),scale="row",
   trace="none",margins = c(6,20), cexRow=.4, cexCol=.8)
  mtext(SETNAME)
}
dev.off()

