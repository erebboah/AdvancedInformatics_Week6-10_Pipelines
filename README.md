# Advanced Informatics Week 6 Exercises

The goal of this week was to organize raw data using symlinks.

I used R to read in a metadata file (in the case of the RNA-seq data), or simply used the existing paths to make a metadata file indicating the new names for each file.

I looped through the files and used `system(paste0())` to paste the `ln -s` command to the terminal and make symlinks for each file.

The scripts were run as such:
```
Rscript get_symlinks_dna.R
Rscript get_symlinks_atac.R
Rscript get_symlinks_rna.R
```

The scripts output metadata files with the new names for the samples are in the directories `DNAseq`, `ATACseq`, or `RNAseq` with the symlinks in sub-directories named `data`:
```
DNAseq/
    dna_samples.txt
    data/
        ADL06_1_1.fq.gz
        ADL06_1_2.fq.gz
        ADL06_2_1.fq.gz
        ADL06_2_2.fq.gz
        ADL06_3_1.fq.gz
        ADL06_3_2.fq.gz
        ...
ATACseq/
    atac_samples.txt
    data/
        P004_R1.fq.gz
        P004_R2.fq.gz 
        P005_R1.fq.gz
        P005_R2.fq.gz
        P006_R1.fq.gz
        P006_R2.fq.gz
        ...
RNAseq/
    rna_samples.txt
    data/
        x21001B0_R1.fq.gz
        x21001B0_R2.fq.gz
        x21001E0_R1.fq.gz
        x21001E0_R2.fq.gz
        x21001H0_R1.fq.gz
        x21001H0_R2.fq.gz
        ...
```

I also ran `fastqc` on one raw data file from the DNA-seq experiment, `ADL06_1_1`. The results are in the `fastqc` folder and the html file can be viewed [here](http://crick.bio.uci.edu/erebboah/ADL06_1_1_fastqc.html). 

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

The bash scripts that made use of the prefix files are called `align_dna.sh`, `align_atac.sh`, and `align_rna.sh`, and I ran them using:
```
sbatch align_dna.sh
sbatch align_atac.sh
sbatch align_rna.sh
```

The scripts folder also contains a brief README about the various R and bash scripts.

The output is in each subfolder named `mapped`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    mapped/
        ADL06_1.RG.bam
        ADL06_1.RG.bam.bai
        ADL06_1.sort.bam
        ADL06_2.RG.bam
        ADL06_2.RG.bam.bai
        ADL06_2.sort.bam
        ADL06_3.RG.bam
        ADL06_3.RG.bam.bai
        ADL06_3.sort.bam
        ...
ATACseq/
    atac_samples.txt
    prefixes.txt
    data/
    mapped/
        P004.sort.bam
        P005.sort.bam
        P006.sort.bam
        ...
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_tissueE.txt
    data/
    mapped/
        x21001E0.sorted.bam
        x21001E0.sorted.bam.bai
        x21002E0.sorted.bam
        x21002E0.sorted.bam.bai
        x21012E0.sorted.bam
        x21012E0.sorted.bam.bai
        ...
```
