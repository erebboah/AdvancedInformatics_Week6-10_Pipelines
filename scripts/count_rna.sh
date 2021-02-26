#!/bin/bash
#SBATCH --job-name=rna_count      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-100         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs
#SBATCH --output=rna_count-%J.out ## output log file
#SBATCH --error=rna_count-%J.err ## error log file

module load subread/2.0.1
gtf="/data/class/ecoevo283/erebboah/ref/dmel-all-r6.13.gtf"
inpath="/data/class/ecoevo283/erebboah/RNAseq/"
file=$inpath"prefixes_random.txt" # prefix file
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# -p paired end, -T ?? num threads, -t Specify feature type(s) in a GTF annotation, -g Specify attribute type in GTF annotation, -Q minimum mapping quality score, -F format of the provided annotation file, -a gene annot, -o output file name

featureCounts -p -T 5 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o ${inpath}counts/${prefix}.counts.txt ${inpath}mapped/${prefix}.sorted.bam
