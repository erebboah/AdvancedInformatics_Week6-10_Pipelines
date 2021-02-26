#!/bin/bash
#SBATCH --job-name=dna_merge_snps      ## Name of the job.
#SBATCH -A SEYEDAM_LAB       ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --array=1-4         ## number of tasks to launch
#SBATCH --cpus-per-task=2    ## number of cores the job needs
#SBATCH --output=dna_merge_snps-%J.out ## output log file
#SBATCH --error=dna_merge_snps-%J.err ## error log file

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.10.2
module load python/3.8.0
module load java/1.8.0
module load gatk/4.1.9.0
module load picard-tools/1.87
module load bamtools/2.5.1        
module load freebayes/1.3.2       
module load vcftools/0.1.16

inpath="/data/class/ecoevo283/erebboah/DNAseq/"
file=$inpath"prefixes2.txt" # prefix file
ref="/data/class/ecoevo283/erebboah/ref/dmel-all-chromosome-r6.13.fasta"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

# merge lanes within samples
# this is an incredible useful shell trick you wouldn't learn for a few years = $(printf 'I=%s ' ${prefix}*.RG.bam)
# what is it doing?
#echo $(printf 'I=%s ' ${prefix}*.RG.bam)

cd /data/class/ecoevo283/erebboah/DNAseq/mapped

java -jar /opt/apps/picard-tools/1.87/MergeSamFiles.jar $(printf 'I=%s ' ${prefix}*.RG.bam) O=${prefix}.merge.bam

# remove duplicates
/opt/apps/gatk/4.1.9.0/gatk MarkDuplicatesSpark -I ${prefix}.merge.bam -O ../gatk/${prefix}.dedup.bam -M ../gatk/${prefix}.dedup.metrics.txt

# no recalibration unless dealing with human data (low SNP density, many HQ known SNPs)
# call genotype on each sample (call SNPs)
/opt/apps/gatk/4.1.9.0/gatk HaplotypeCaller -R $ref -I ../gatk/${prefix}.dedup.bam --minimum-mapping-quality 30 -ERC GVCF -O ../gatk/${prefix}.g.vcf.gz --sample-name ${prefix}


