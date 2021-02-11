# Get list of fastq.gz files from ATACseq and make symlinks
# Files are almost named nicely, just remove first part

system("find /data/class/ecoevo283/public/RAWDATA/ATACseq/*.fq.gz > /data/class/ecoevo283/erebboah/ATACseq/atac_samples.txt")

samples = as.data.frame(read.table("/data/class/ecoevo283/erebboah/ATACseq/atac_samples.txt", stringsAsFactors = F))
samplename = sapply(strsplit(samples[,1], "4R009_L1_"), "[[", 2)
samples$name = samplename
colnames(samples) = c("data_location","sample_name")
write.table(samples, file="/data/class/ecoevo283/erebboah/ATACseq/atac_samples.txt",sep="\t", quote=F, row.names = F)

for (i in 1:length(samplename))
{
  system(paste0("ln -s ", samples[i,1]," /data/class/ecoevo283/erebboah/ATACseq/data/",samplename[i]))
}

