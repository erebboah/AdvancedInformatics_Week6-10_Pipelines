mytab = read.table("/data/class/ecoevo283/erebboah/RNAseq/rna_samples.txt",header=TRUE)
write.table(mytab$renamed_SampleName[mytab$TissueCode=="E"],file="/data/class/ecoevo283/erebboah/RNAseq/prefixes_tissueE.txt", col.names = F, row.names = F, quote = F)
