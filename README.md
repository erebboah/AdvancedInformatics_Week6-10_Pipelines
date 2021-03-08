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

This week, we mapped the 3 datasets to the reference genome.

On an interactive node, I indexed the reference file here: `ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"` for BWA, Picard, and HISAT2. 
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
java -d64 -Xmx128g -jar /opt/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=$ref O=ref/dmel-all-chromosome-r6.13.dict
hisat2-build $ref $ref
```

I made "prefix" files to use for the bash scripts with the provided code from Dr. Long's notes. For RNA-seq, I made a subsetted prefix file using `scripts/subset_rna.R`, containing 100 random samples.
```
ls /data/class/ecoevo283/erebboah/DNAseq/data/*1.fq.gz | sed 's/_1.fq.gz//' > ../prefixes.txt
ls /data/class/ecoevo283/erebboah/ATACseq/data/*R1.fq.gz | sed 's/_R1.fq.gz//' > ../prefixes.txt
ls /data/class/ecoevo283/erebboah/RNAseq/data/*R1.fq.gz | sed 's/_R1.fq.gz//' > ../prefixes.txt
Rscript subset_rna.R
```
The prefix files are in the `DNAseq`, `ATACseq`, and `RNAseq` folders.

The bash scripts that made use of the prefix files are called `align_dna.sh`, `align_atac.sh`, and `align_rna.sh`. DNA-seq data is aligned using BWA followed by generation of read groups using Picard. ATAC-seq data data is simply aligned using BWA. RNA-seq data is aligned using HISAT2 and the resulting `SAM` file is converted to `BAM` and sorted. I ran them with the following commands:
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
    data/
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
    prefixes_random.txt
    data/
    mapped/
        x21001B0.sorted.bam
        x21001B0.sorted.bam.bai
        x21001B0.sorted.bam
        x21001H0.sorted.bam.bai
        x21001H0.sorted.bam
        x21001H0.sorted.bam.bai
        ...
```

# Advanced Informatics Week 8 Exercises

The goals for this week were to:
1. work through the DNA-seq analysis pipeline to call SNPs using GATK and generate `VCF` files and 
2. generate counts per gene for RNA-seq data using the subread package.

## DNA-seq GATK pipeline
I made another set of "prefix" files to use for the first step, which is now on a per-sample basis with merged replicates, thus 4 samples: `ADL06`, `ADL09`, `ADL10`, and `ADL14`.
```
ls /data/class/ecoevo283/erebboah/DNAseq/data/*_1_1.fq.gz | sed 's/_1.fq.gz//' > ../prefixes2.txt
```

The prefix file is in the `DNAseq` folder.
I used the code provided in Dr. Long's notes to make 3 bash scripts. I ran the first script to merge sample replicates with Picard `MergeSamFiles`, remove duplicates with GATK `MarkDuplicatesSpark`, and call SNPs with GATK `HaplotypeCaller` to generate one `GVCF` file per sample.
```
sbatch dna_gatk_step1.sh
```

The output is in `DNAseq/gatk`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    prefixes2.txt
    data/
    mapped/
    gatk/
        ADL06.dedup.bam
        ADL06.dedup.bam.bai
        ADL06.dedup.bam.sbi
        ADL06.dedup.metrics.txt
        ADL06.g.vcf.gz
        ADL06.g.vcf.gz.tbi
        ...
```

Next, I ran the second script to combine the `GVCF` output per sample into 1 unified file using GATK `CombineGVCFs`:
```
sbatch dna_gatk_step2.sh
```
The output is `DNAseq/gatk/allsample.g.vcf.gz`. (And index `DNAseq/gatk/allsample.g.vcf.gz.tbi`.)

Finally, the third script performs joint genotyping on the combined `GVCF` file to output a final `VCF` using GATK `GenotypeGVCFs`.
```
sbatch dna_gatk_step3v1.sh
```
The output is `DNAseq/gatk/result.vcf.gz`. (And index `DNAseq/gatk/result.vcf.gz.tbi`.)

We can also call SNPs by intervals to speed up the process. 

Create 10Mb regions using the provided support script `fasta_generate_regions.py`
```
python3 fasta_generate_regions.py ../ref/dmel-all-chromosome-r6.13.fasta 10000000 > ../ref/my_regions_10Mb.txt
```

I ran the following script which uses the `--intervals` option in GATK with the `my_regions_10Mb.txt` file to output a `.vcf.gz` and index for each 10Mb region (1,881 regions total).
```
sbatch dna_gatk_step3v2.sh
```

The output is in `DNAseq/gatk/SNPbyregion`:
```
DNAseq/
    dna_samples.txt
    prefixes.txt
    prefixes2.txt
    data/
    mapped/
    gatk/
       SNPbyregion/
            211000022278031:1-1021.vcf.gz
            211000022278031:1-1021.vcf.gz.tbi
            211000022278032:1-6936.vcf.gz
            211000022278032:1-6936.vcf.gz.tbi
            211000022278033:1-1143.vcf.gz
            211000022278033:1-1143.vcf.gz.tbi
            ...
```


## RNA-seq counts matrix generation
Following the format of the scripts provided for the DNA-seq pipeline, I wrote a bash script to count the number of reads per gene from the 100 sorted `BAM` files using `featureCounts` from the subreads package.
```
sbatch count_rna.sh
```

The output is in `RNAseq/counts`:
```
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_random.txt
    data/
    mapped/
    counts/
        x21001B0.counts.txt
        x21001B0.counts.txt.summary
        x21001H0.counts.txt
        x21001H0.counts.txt.summary
        x21001P0.counts.txt
        x21001P0.counts.txt.summary
        ...
```

I concatenated the `*.counts.txt` files into a matrix with each sample as a column and each gene as a row (100 samples x 17,491 genes) using `scripts/make_rna_matrix.R`. The script loops though each sample using the `prefixes_random.txt` file and joins each `.counts.txt` file by the FlyBase gene ID.
```
Rscript make_rna_matrix.R
```
The output is `RNAseq/counts/counts_100samples_flybaseIDs.tsv`.

For next week, I installed `DESeq2`, `GenomicFeatures`, `Rsamtools`, and `GenomicAlignments` in my R conda enviroment. Some packages installed easily with `install.packages()` and some like `DESeq2` I installed using [conda](https://anaconda.org/bioconda/bioconductor-deseq2), `conda install -c bioconda bioconductor-deseq2`.

# Advanced Informatics Week 9 Exercises

The goals for this week were to:
1. analyze RNA-seq,
2. make some plots of the RNA-seq results, and
3. visualize ATAC-seq using IGV (Integrative Genomics Viewer). 

## RNA-seq analysis
First, I copied the metadata to my directory, `RNAseq/RNAseq384_SampleCoding.txt`.

I made DESeq2 objects using Dr. Long's provided code. 

We are roughly following section 2.4.3, Starting from count tables on [page 10](https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf). 
The `countdata` is `RNAseq/counts/counts_100samples_flybaseIDs.tsv` and the `coldata` (metadata) is `RNAseq/RNAseq384_SampleCoding.txt`, subsetted to match the 100 samples I processed in order of the columns in `countdata`. 
I used `design=~TissueCode`, which tells DESeq2 that the experimental design is simply based on different tissues. There can be complicated designs with multiple experimental variabes. 

Next the script filters out lowly expressed genes where the sum of the counts across all 100 samples is less than 10, before running DESeq2 in a single command: `dds <- DESeq(dds)`.

DESeq2 does the following:
```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 9 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
```

Inspecting the result with `results(dds)` shows log2 fold change, p value, and other useful plotting parameters the last variable in the design formula. 
Specific variables can be called by `results(dds,contrast = c("TissueCode","P","E"))`, for example to compare tissues P and E. 

Then the data was rlog-transformed, which returns a SummarizedExperiment object containing the transformed values in the assay slot.

```
rld = rlog( dds )
mydata = assay(rld)
```

This command took a while; the package recommends using `vst()` instead of `rlog()` for a large number of samples.

To run the script:

```
Rscript rna_deseq2_make_objects.R
```

The output is in `RNAseq/analysis`:
```
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_random.txt
    data/
    mapped/
    counts/
    analysis/
        dds.rda
        rld.rda
```

## RNA-seq analysis results plotting
Next, I made plots with DESeq2, heatmap.2, and EnhancedVolcano, which conveniently takes in a DESeq2 SummarizedExperiment object. The script makes 7 plots. 
```
Rscript rna_deseq2_plots.R
```

The output is in `RNAseq/figures`:
```
RNAseq/
    rna_samples.txt
    prefixes.txt
    prefixes_random.txt
    data/
    mapped/
    counts/
    analysis/
    figures/
        dispEst.png
        heatmap.png
        heatmap_topvargenes.png
        hist.png
        pca.png
        plotMA.png
        volcano_tissuePvstissueB.png
```


The x axis of the [MA plot](https://en.wikipedia.org/wiki/MA_plot) is the average expression over all samples and the y axis the log2 fold change between P and B (last variable comparison by default).
It shows that most significantly differentially expressed (DE) genes have medium expression values. 

![ma plot](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/plotMA.png?raw=true)

The gene dispersion estimate plot shows how the dispersion estimates for the genes (black dots) will be shrunk towards the red trend line. 

![dispersion plot](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/dispEst.png?raw=true)

The histogram of p-values shows that most genes returned by the test for differential expression are significantly differentially expressed.

![histogram](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/hist.png?raw=true)

A heatmap of the [distance matrix](https://en.wikipedia.org/wiki/Distance_matrix) shows a clear distinction between tissues.

![heatmap](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/heatmap.png?raw=true)

Principal component analysis (PCA) of the data labeled by tissue type shows that the tissues are clearly separated by gene expression. 

![pca](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/pca.png?raw=true)

A heatmap of top variable gene expression again shows distinct gene expression in each tissue type. 

![heatmap top var](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/heatmap_topvargenes.png?raw=true)

A volcano plot of genes differentially expressed in tissue P (positive LFC) vs. tissue B (negative LFC) shows that many are significant by both p-value and fold change (default 1e-05 and 1, respectively).
![volcano](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/RNAseq/figures/volcano_tissuePvstissueB.png?raw=true)

## ATAC-seq coverage on UCSC
### Generate bigwig files
I wrote a bash script to convert `bam` to `bedgraph` using `genomeCoverageBed`. Using `-scale` the bedgraphs can be scaled by read depth. The number of reads can be determined by `samtools view -c -F 4` where [4](https://www.samformat.info/sam-format-flag) indicates unmapped reads and [F](http://www.htslib.org/doc/samtools-view.html) indicates to count all reads EXCEPT those with the flag. Samtools is used to pipe the bam file into the UCSC tool. There is no need for a reference when running `genomeCoverageBed` from a `bam` file. 

```
Nreads=`samtools view -c -F 4 ${inpath}mapped/${prefix}.sort.bam`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`

samtools view -b ${inpath}mapped/${prefix}.sort.bam | genomeCoverageBed -ibam - -bg -scale $Scale > ${inpath}coverage/${prefix}.coverage
```

Next I sorted the `bedgraph` by `sortBed` before running `bedGraphToBigWig` to make the final `bigwig` file. 

To make a genome for `bedGraphToBigWig`:
```
cut -f 1,2 /data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta.fai > /data/class/ecoevo283/erebboah/ref/dm6.chrom.sizes
```

Finally I remove the unsorted `bedgraph` files to save space.

To run the script:
```
sbatch bigwig_atac.sh
```

The output is in `ATACseq/coverage`:
```
ATACseq/
    atac_samples.txt
    prefixes.txt
    data/
    mapped/
    coverage/
        P004.bw
        P004.sort.coverage
        P005.bw
        P005.sort.coverage
        P006.bw
        P006.sort.coverage
        ...
```
### Upload to UCSC genome browser
I used our lab's [website](http://crick.bio.uci.edu/erebboah/class/) to host the `bigwig` files.
```
scp *.bw erebboah@crick.bio.uci.edu:/var/www/html/erebboah/class
```

I added 25 ATAC-seq files via the [Upload custom tracks](https://genome.ucsc.edu/cgi-bin/hgCustom) tool on UCSC. The session can be shared with others via this [link](https://genome.ucsc.edu/s/erebboah/atac_coverage_dm6).

A screenshot of a random section of the Drosophila genome (chr2L:1,030,879-1,204,340) with the ATAC-seq `bigwig` tracks can be found at `ATACseq/figures/chr2L-1030879-1204340.png`:

![UCSC](https://github.com/erebboah/AdvancedInformatics_Week6-10_Pipelines/blob/main/ATACseq/figures/chr2L-1030879-1204340.png?raw=true)
