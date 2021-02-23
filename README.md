# Advanced Informatics Week 6 Exercises

The goal of this week was to organize raw data using symlinks.

I used R to read in a metadata file (in the case of the RNA-seq data), or simply used the existing paths to make a metadata file indicating the new names for each file.

I looped through the files and used `system(paste0())` to paste the `ln -s` command to the terminal and make symlinks for each file.

The metadata files with the new names for the samples will be in the directories `DNAseq`, `ATACseq`, or `RNAseq` with the symlinks in sub-directories named `data`.

I ran `fastqc` on one raw data file from the DNA-seq experiment, `ADL06_1_1`. The results are in the `fastqc` folder and the html file can be viewed [here](http://crick.bio.uci.edu/erebboah/ADL06_1_1_fastqc.html). 

# Advanced Informatics Week 7 Exercises

This week, we mapped the 3 datasets to the referene genome, using BWA for DNA-seq and ATAC-seq and HISAT2 for RNA-seq (a subset of samples).

On an interactive node, I indexed the reference file here: `ref="/data/class/ecoevo283/erebboah/dmel-all-chromosome-r6.13.fasta"` for BWA, Picard, and HISAT2. 
```
module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load java/1.8.0
module load hisat2/2.2.1

bwa index $ref
samtools faidx  $ref
java -d64 -Xmx128g -jar /opt/apps/picard-tools/1.87/reateSequenceDictionary.jar R=$ref O=ref/dmel-all-chromosome-r6.13.dict
hisat2-build $ref $ref
```

I made "prefix" files to use for the bash scripts with the provided code (and a subsetted prefix file using an R script in the `scripts` folder):
```
ls /data/class/ecoevo283/erebboah/DNAseq/data/*1.fq.gz | sed 's/_1.fq.gz//' > ../prefixes.txt
ls /data/class/ecoevo283/erebboah/ATACseq/data/*R1.fq.gz | sed 's/_R1.fq.gz//' > ../prefixes.txt
ls /data/class/ecoevo283/erebboah/RNAseq/data/*R1.fq.gz | sed 's/_R1.fq.gz//' > ../prefixes.txt
```
The prefix files are in the `DNAseq`, `ATACseq`, and `RNAseq` folders.

The bash scripts that made use of the prefix files are called `align_dna.sh`, `align_atac.sh`, and `align_rna.sh`.

The scripts folder also contains a brief README about the various R and bash scripts.
