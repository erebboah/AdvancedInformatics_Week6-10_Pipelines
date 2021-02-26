#!/bin/bash
#SBATCH --job-name=dna_combineGVCF      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8    ## number of cores the job needs
#SBATCH --output=dna_combineGVCF-%J.out ## output log file
#SBATCH --error=dna_combineGVCF-%J.err ## error log file

module load java/1.8.0
module load gatk/4.1.9.0

ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"
cd /data/class/ecoevo283/erebboah/DNAseq/gatk
# another dirty shell trick, you can't printf -V...
/opt/apps/gatk/4.1.9.0/gatk CombineGVCFs -R $ref $(printf -- '-V %s ' *.g.vcf.gz) -O allsample.g.vcf.gz
## and merge gVCFs 
