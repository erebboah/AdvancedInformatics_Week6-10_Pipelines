#!/bin/bash
#SBATCH --job-name=dna_callSNPs      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8    ## number of cores the job needs
#SBATCH --output=dna_callSNPs-%J.out ## output log file
#SBATCH --error=dna_callSNPs-%J.err ## error log file

module load java/1.8.0
module load gatk/4.1.9.0

ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"
cd /data/class/ecoevo283/erebboah/DNAseq/gatk

/opt/apps/gatk/4.1.9.0/gatk GenotypeGVCFs -R $ref -V allsample.g.vcf.gz -stand-call-conf 5 -O result.vcf.gz
## end call SNPs
