# Get list of fastq.gz files from DNAseq and make symlinks
# Files are already renamed, no need to change names

system("find /data/class/ecoevo283/public/RAWDATA/DNAseq/*.fq.gz > /data/class/ecoevo283/erebboah/DNAseq/dna_samples.txt")

samples = as.data.frame(read.table("/data/class/ecoevo283/erebboah/DNAseq/dna_samples.txt", stringsAsFactors = F))
samplename = sapply(strsplit(samples[,1], "DNAseq/"), "[[", 2)
samples$name = samplename
colnames(samples) = c("data_location","sample_name")
write.table(samples, file="/data/class/ecoevo283/erebboah/DNAseq/dna_samples.txt",sep="\t", quote=F, row.names = F)

for (i in 1:length(samplename))
{
  system(paste0("ln -s ", samples[i,1]," /data/class/ecoevo283/erebboah/DNAseq/data/",samplename[i]))
}

