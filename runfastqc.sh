#!/bin/bash
#PBS -l nodes=1:ppn=2,walltime=24:00:00,vmem=30g,mem=220g
#PBS -N quality
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com

#This batch to run quality control of FASQT files##

set -eux

module load java/7
module load fastqc/0.11.5

## H4H directories
## Misc tools parameters
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')

scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/DATA
workdir=$scratch/sickkids1/analysis/$pbs
home=/cluster/home/itagoug
user_databases=$scratch/databases
human.ref.genome=$scratch/human/hg19/GCF_000001405.25_GRCh37.p13_genomic.fna




reads[1]=$raw/9_A02_R1.fastq.gz
reads[2]=$raw/9_A02_R2.fastq.gz

mkdir -p $workdir
cd $workdir

for i in {1..2}; do
	zcat ${reads[$i]} | fastqc ${reads[$i]} --outdir=$workdir
done
