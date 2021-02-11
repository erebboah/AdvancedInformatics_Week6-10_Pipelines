# Get list of fastq.gz files from RNAseq and make symlinks
# Files are hot off the sequencer nestled in multiple directories

samples = as.data.frame(read.table("/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384_SampleCoding.txt", header = T))
samples$renamed_SampleName = paste0("x",samples$FullSampleName)
write.table(samples, file="/data/class/ecoevo283/erebboah/RNAseq/rna_samples.txt",sep="\t", quote=F, row.names = F)

plex = sapply(strsplit(samples$Multiplexi5index, "_"), "[[", 1)
for (i in 1:length(samples$renamed_SampleName))
{
	path = paste0("/data/class/ecoevo283/public/RAWDATA/RNAseq/RNAseq384plex_flowcell01/Project_",
		plex[i],"/Sample_",samples$SampleNumber[i])
	system(paste0("ln -s ", path, "/*_R1_001.fastq.gz",
		" /data/class/ecoevo283/erebboah/RNAseq/data/",samples$renamed_SampleName[i],"_R1.fq.gz"))
        system(paste0("ln -s ", path, "/*_R2_001.fastq.gz",
                " /data/class/ecoevo283/erebboah/RNAseq/data/",samples$renamed_SampleName[i],"_R2.fq.gz"))
}

