#!/bin/bash
#PBS -l nodes=1:ppn=6,walltime=48:00:00,vmem=30g,mem=220g
#PBS -N calling
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com

##Task CNV calling using CNVkit

set -eux

module load java/8 bwa/0.7.15 samtools/1.9
module load snpEff/4.3
module load python3/3.4.3
module load R/3.5.0
module load CNVkit/0.9.3

## H4H directories
scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/bam_clean/dup_clean


## __DO NOT CHANGE__
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')
bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
bam_clean=$scratch/sickkids1/bam_clean
mills=$scratch/human/gold/Mills.gold.standard.adapted_noChr.vcf
agilent=$scratch/human/agilent_original/74_clean.bed
snpeff_config=$home/scripts/snpEff.config


mkdir -p $workdir
cd $workdir

	cnvkit.py batch $raw/206_P6.final.bam --normal $raw/206_Bu.final.bam \
		  --targets $agilent \
		  --fasta $human_ref_genome \
		  --short-names \
		  --output-reference $workdir/my_reference.cnn \
		  -d $workdir/output.2.cnv.txt


echo -e "\n\nAnalyses are done, without complaints."
