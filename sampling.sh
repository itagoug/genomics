#!/bin/bash
#PBS -l nodes=1:ppn=32,walltime=48:00:00,vmem=30g,mem=220g
#PBS -N sampling reads
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com

## Task: sampling fastq to run it during debbuging 
set -eux

module load seqtk

## H4H directories
scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/bam_clean/set1
data=$scratch/sickkids1/DATA/custom
label=$raw/217_P6


# __CHANGE__
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')
bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
bam_clean=$scratch/sickkids1/bam_clean
mills=$scratch/human/gold/Mills.gold.standard.adapted_noChr.vcf
agilent=$scratch/human/agilent_original/206.clean.bed
snpeff_config=$home/scripts/snpEff.config
LocatIt=$home/AGeNT

mkdir -p $workdir
cd $workdir


seqtk sample -s1234567 $raw/217_P6_R1.fastq 500000 > $workdir/217_test_R1.fastq
seqtk sample -s1234567 $raw/217_P6_R2.fastq 500000 > $workdir/217_test_R2.fastq

echo -e "\n\nAnalyses are done, without complaints."
