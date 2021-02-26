#!/bin/bash
#SBATCH --job-name=dna_callSNPs_interval      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-1881         ## number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs
#SBATCH --output=dna_callSNPs_interval-%J.out ## output log file
#SBATCH --error=dna_callSNPs_interval-%J.err ## error log file

module load java/1.8.0
module load gatk/4.1.9.0

ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"
cd /data/class/ecoevo283/erebboah/DNAseq/gatk
interval=`head -n $SLURM_ARRAY_TASK_ID ../../ref/my_regions_10Mb.txt | tail -n 1`

/opt/apps/gatk/4.1.9.0/gatk GenotypeGVCFs -R $ref -V allsample.g.vcf.gz --intervals $interval -stand-call-conf 5 -O SNPbyregion/${interval}.vcf.gz
## end second approach to calling SNPs  (you have to put them back together...)
