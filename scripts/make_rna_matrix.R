library(data.table)

setwd("/data/class/ecoevo283/erebboah/RNAseq")

#################### Load metadata #################### 
samples = read.table("/data/class/ecoevo283/erebboah/RNAseq/prefixes_tissueE.txt",header=FALSE)
samples = samples$V1

#################### Read in .txt and make counts matrix across all samples in data folder ####################
thisMat=as.data.frame(read.table(paste0("/data/class/ecoevo283/erebboah/RNAseq/counts/",samples[1],".counts.txt"),sep="",head=T, stringsAsFactors = F))
head(thisMat)
colnames(thisMat) = c("GeneID","Chr","Start","End","Strand","Length","SampleID")
counts=thisMat[,c('GeneID','SampleID')]
head(counts)

for(i in 2:length(samples)){
  thisMat=as.data.frame(read.table(paste0("/data/class/ecoevo283/erebboah/RNAseq/counts/",samples[i],".counts.txt"),sep="",head=T, stringsAsFactors = F))
  colnames(thisMat) = c("GeneID","Chr","Start","End","Strand","Length","SampleID")
  thisMat=thisMat[na.omit(match(counts$GeneID,thisMat$GeneID)),]
  counts=as.data.frame(cbind(counts,thisMat[,"SampleID"]))
}

dim(counts)
rownames(counts) = counts$GeneID
counts$GeneID = NULL
colnames(counts) = samples
head(counts) 

# Save counts matrix with FlyBase IDs
write.table(counts, file="/data/class/ecoevo283/erebboah/RNAseq/counts/counts_tissueE_flybaseIDs.tsv", quote=F)

