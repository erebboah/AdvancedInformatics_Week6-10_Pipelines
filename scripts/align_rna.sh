#!/bin/bash
#SBATCH --job-name=rna_align      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-79         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs
#SBATCH --output=rna_align-%J.out ## output log file
#SBATCH --error=rna_align-%J.err ## error log file

module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load java/1.8.0
module load hisat2/2.2.1


inpath="/data/class/ecoevo283/erebboah/RNAseq/"
file=$inpath"prefixes_tissueE.txt" # prefix file
ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# alignments
hisat2 -p 2 -x $ref -1 ${inpath}data/${prefix}_R1.fq.gz -2 ${inpath}data/${prefix}_R2.fq.gz -S ${inpath}mapped/${prefix}.sam
samtools view -bS ${inpath}mapped/${prefix}.sam > ${inpath}mapped/${prefix}.bam
samtools sort ${inpath}mapped/${prefix}.bam ${inpath}mapped/${prefix}.sorted.bam
samtools index ${inpath}mapped/${prefix}.sorted.bam
rm ${inpath}mapped/${prefix}.sam
rm ${inpath}mapped/${prefix}.bam
