The following R scripts make symlinks in the data folders for ATAC-seq, DNA-seq, and RNA-seq:
```
get_symlinks_atac.R
get_symlinks_dna.R
get_symlinks_rna.R
```

`subset_rna.R` is a simple R script to subset the prefixes used in mapping the RNA-seq data to only one tissue type (since there are so many samples to map).

`align_atac.sh` and `align_dna.sh` are bash scripts to map ATAC-seq and DNA-seq using BWA and for DNA-seq, add `RG` (read group) tags using Picard to use in GATK.
