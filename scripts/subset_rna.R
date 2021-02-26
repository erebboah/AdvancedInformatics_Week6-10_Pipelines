mytab = read.table("/data/class/ecoevo283/erebboah/RNAseq/rna_samples.txt",header=TRUE)
random_samples = sample(mytab$SampleNumber,100)
write.table(mytab$renamed_SampleName[random_samples],file="/data/class/ecoevo283/erebboah/RNAseq/prefixes_random.txt", col.names = F, row.names = F, quote = F)
