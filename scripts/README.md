### R scripts 
The following scripts make symlinks in the data folders for ATAC-seq, DNA-seq, and RNA-seq:
```
get_symlinks_atac.R
get_symlinks_dna.R
get_symlinks_rna.R
```

`subset_rna.R` is a simple script to subset the prefixes used in mapping the RNA-seq data to only one tissue type (since there are so many samples to map).

### Bash scripts
`align_atac.sh` and `align_dna.sh` align ATAC-seq and DNA-seq using BWA v.0.7.0 and for DNA-seq, add `RG` (read group) tags using Picard 1.87 (and Java 1.8.0) to use in GATK. 
The `fastq` reads become `sam` files when aligned, which are converted to `bam` files and sorted using samtools 1.10 for downstream applications. 
The readgroups `bam` file is also indexed with samtools.

`align_rna.sh` aligns RNA-seq using HISAT2 v.2.2.1, converts `sam` to  `bam`, sorts, and indexes.

The scripts also remove intermediate files such as `sam` and unsorted `bam` to save space.
