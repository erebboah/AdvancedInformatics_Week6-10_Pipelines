# Advanced Informatics Week 6 Exercises

The goal of this week was to organize raw data using symlinks.

I used R to read in a metadata file (in the case of the RNA-seq data), or simply used the existing paths to make a metadata file indicating the new names for each file.

I looped through the files and used `system(paste0())` to paste the `ln -s` command to the terminal and make symlinks for each file.

The metadata files with the new names for the samples will be in the directories `DNAseq`, `ATACseq`, or `RNAseq` with the symlinks in sub-directories named `data`.

I ran `fastqc` on one raw data file from the DNA-seq experiment, `ADL06_1_1`. The results are in the `fastqc` folder and the html file can be viewed [here](http://crick.bio.uci.edu/erebboah/ADL06_1_1_fastqc.html). 
