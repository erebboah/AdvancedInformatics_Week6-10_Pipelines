library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(EnhancedVolcano)

load("/data/class/ecoevo283/erebboah/RNAseq/analysis/rld.rda")
load("/data/class/ecoevo283/erebboah/RNAseq/analysis/dds.rda")
res <- results( dds )

fname_MA ="/data/class/ecoevo283/erebboah/RNAseq/figures/plotMA.pdf"
fname_dispEst="/data/class/ecoevo283/erebboah/RNAseq/figures/dispEst.pdf"
fname_hist="/data/class/ecoevo283/erebboah/RNAseq/figures/hist.pdf"
fname_pca="/data/class/ecoevo283/erebboah/RNAseq/figures/pca.pdf"
fname_heatmap="/data/class/ecoevo283/erebboah/RNAseq/figures/heatmap.pdf"
fname_heatmap_topvargenes="/data/class/ecoevo283/erebboah/RNAseq/figures/heatmap_topvargenes.pdf"
fname_volcano_pvsb="/data/class/ecoevo283/erebboah/RNAseq/figures/heatmap_volcano_tissuePvstissueB.pdf"

pdf(file=fname_MA,
    width = 5,
    height = 5) 
plotMA( res, ylim = c(-1, 1), colNonSig = "gray60",
       colSig = "blue",)
dev.off()

pdf(file=fname_dispEst,
    width = 5,
    height = 5) 
plotDispEsts( dds )
dev.off()

pdf(file=fname_hist,
    width = 5,
    height = 5) 
hist( res$pvalue, breaks=20, col="grey" )
dev.off()

sampleDists = dist( t( assay(rld) ) )
# heat map
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(file=fname_heatmap,
    width = 5,
    height = 5) 
heatmap.2( sampleDistMatrix, trace="none", col=colours)
dev.off()

# PCs
# wow you can sure tell tissue apart
pdf(file=fname_pca,
    width = 5,
    height = 5) 
print( plotPCA( rld, intgroup = "TissueCode") )
dev.off()

# heat map with gene clustering
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
pdf(file=fname_heatmap_topvargenes,
    width =5,
    height = 5) 
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

# volcano plot this is an exercise
res <- results(dds,
    contrast = c('TissueCode','P','B'))
res <- lfcShrink(dds,
    contrast = c('TissueCode','P','B'), res=res, type = 'normal')

pdf(file=fname_volcano_pvsb,
    width =10,
    height = 10)     
 EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
dev.off()

