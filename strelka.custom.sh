#!/bin/bash
#PBS -l nodes=1:ppn=6,walltime=48:00:00,vmem=30g,mem=220g
#PBS -N calling
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com


set -eux

module load java/8 bwa/0.7.15 samtools/1.9
module load snpEff/4.3
module load python/2.7


## H4H directories
scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/bam_clean/dup_clean

## __DO NOT CHANGE__
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')
bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
#index human_ref_genome
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
bam_clean=$scratch/sickkids1/bam_clean
mills=$scratch/human/gold/Mills.gold.standard.strelka.adapted_noChr.vcf.gz
#index bed file
agilent=$scratch/human//agilent_original/206.clean.id.bed.gz
snpeff_config=$home/scripts/snpEff.config
strelka=$home/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py
## parameters
pipeline=strelka
gMUM=30g
threads=25
snpeff_genome=GRCh37.75



#################################
## calling copy number variants##
#################################


mkdir -p $workdir
cd $workdir

# __change __

sample[1]=$raw/217_P6.final.bam
sample[2]=$raw/217_bu.final.bam



samtools index ${sample[1]}
samtools index ${sample[2]}

# Configure strelka
mkdir -p $workdir/config

${strelka} \
  --normalBam=${sample[2]} \
  --tumorBam=${sample[1]} \
  --referenceFasta=$human_ref_genome \
  --exome \
  --callRegions=$agilent \
  --runDir=$workdir/config


$workdir/config/runWorkflow.py \
	-m local

# annotate

cd $workdir/config/results/variants/
somatic=$workdir/config/results/variants/somatic.snvs.vcf
output=$workdir/config/results/variants/somatic.snvs.annotated
gunzip ${somatic}.gz

java -jar $snpeff_dir/SnpSift.jar annotate \
       -noLog -noDownload \
       -c $snpeff_config \
       -clinvar \
       ${somatic} \
       | java -jar $snpeff_dir/SnpSift.jar annotate \
       		 -noLog -noDownload \
       		 -c $snpeff_config \
       		 -dbsnp - \
       | java -jar $snpeff_dir/snpEff.jar eff \
       		 -noDownload -noLog \
       		 -c $snpeff_config \
       		 -stats ${output}.html \
       		 $snpeff_genome \
       		 > ${output}.vcf


echo -e "\n\nAnalyses are done, without complaints."
