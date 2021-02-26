### R scripts 
#### Week 6 (3)
The following scripts create symlinks in the data folders for ATAC-seq, DNA-seq, and RNA-seq, making use of the  `system(paste0())` command:
```
get_symlinks_atac.R
get_symlinks_dna.R
get_symlinks_rna.R
```

#### Week 7 (1)
`subset_rna.R` is a simple script to subset the prefixes used in mapping the RNA-seq data to only 100 random samples (since there are so many samples to map).

#### Week 8
`make_rna_matrix.R` loops through each of the 100 `.counts.txt` files by their ID in `prefixes_random.txt` and joins each sample by row, their FlyBase gene IDs. The matrix is formatted so that each column is the sample number (preceded by "x" as R does not like variables that start with a number), and each row is a FlyBase gene ID.

### Bash scripts
#### Week 7
`align_atac.sh` and `align_dna.sh` align ATAC-seq and DNA-seq using BWA v.0.7.0 and for DNA-seq, add `RG` (read group) tags using Picard 1.87 (and Java 1.8.0) to use in GATK. Importantly, the RG tag must be the same for all replicates of a sample or GATK will fail at the `HaplotypeCaller` command because it thinks the merged `BAM` is made up of multiple unique samples instead of replicates. I used the same prefix file but added the RG with `prefix:0:5` which subsets the sample name to the first 5 characters, removing the `_1`, `_2`, `_3`.

The `fastq` reads become `SAM` files when aligned, which are converted to `BAM` files and sorted using samtools 1.10 for downstream applications. 
The readgroups `BAM` file is also indexed with samtools.

`align_rna.sh` aligns RNA-seq using HISAT2 v.2.2.1, converts `SAM` to  `BAM`, sorts, and indexes.

The scripts also remove intermediate files such as `SAM` and unsorted `BAM` to save space.

#### Week 8
`dna_gatk_step1.sh` uses Picard v.1.87 to merge the RG bam files and GATK v.4.1.9.0 to mark duplicates and call genotypes on each sample.

`dna_gatk_step2.sh` uses GATK v.4.1.9.0 to combine the `GVCF` files into a merged `GVCF`.

`dna_gatk_step3v1.sh` is the first method to call SNPs using `GenotypeGVCFs` from GATK v.4.1.9.0, resulting in a final `VCF` file.

`dna_gatk_step3v2.sh` is the second method to call SNPs using `GenotypeGVCFs` from GATK v.4.1.9.0, resulting in 1,881 final `VCF` files that need to be merged using `cat`.

`count_rna.sh` counts number of reads per gene across all 100 samples (another array job) using subread v2.0.1 in paired-end mode, with a minimum mapping score of 30.

### Python scripts
#### Week 8 
`fasta_generate_regions.py` is an accessory script provided by Dr. Long that splits the genome into 10Mb regions for use in the interval-based method of SNP calling.
