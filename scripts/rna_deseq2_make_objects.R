# create DEseq2 object & run DEseq

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(EnhancedVolcano)

countdata =read.table("/data/class/ecoevo283/erebboah/RNAseq/counts/counts_100samples_flybaseIDs.tsv", stringsAsFactors=F)
sampleInfo = read.table("/data/class/ecoevo283/erebboah/RNAseq/RNAseq384_SampleCoding.txt", header=T, stringsAsFactors=F)
sampleInfo_100 = sampleInfo[paste0("x",sampleInfo$FullSampleName) %in% colnames(countdata),]
sampleInfo_100=sampleInfo_100[match(colnames(countdata), paste0("x",sampleInfo_100$FullSampleName)),]
#table(colnames(countdata)==paste0("x",sampleInfo_100$FullSampleName))
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo_100, design=~TissueCode)
# filter lowly expressed genes http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results( dds )

###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)
save(rld,file="/data/class/ecoevo283/erebboah/RNAseq/analysis/rld.rda")
save(dds,file="/data/class/ecoevo283/erebboah/RNAseq/analysis/dds.rda")
