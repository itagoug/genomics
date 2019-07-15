#!/bin/bash
#PBS -l nodes=1:ppn=16,walltime=48:00:00,vmem=30g,mem=220g
#PBS -N qlimap
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com


set -eux

module load qualimap/2.2
module load java/8

## H4H directories
### __CHANGE__
scratch=/cluster/projects/tiedemannlab/itagoug
bam_clean=$scratch/sickkids1/bam_clean/set1
raw=$scratch/sickkids1/bam_clean/set1
label=12_D02
sample=$raw/$label*bam


#for i in {9 11 }; do
#label=$i



#label=$raw/9_A02


## __DO NOT CHANGE__
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')
bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
#GRCh37-lite.fa link http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
bam_clean=$scratch/sickkids1/bam_clean
mills=$scratch/human/gold/Mills.gold.standard.adapted_noChr.vcf
#agilent=$scratch/human/agilent_original/padded.clean.bed
agilent_corrected=$scratch/human/agilent_original/padded.clean_qualimap.bed
snpeff_config=$home/scripts/snpEff.config

## parameters
pipeline=bcftools
gMUM=30g
threads=25
snpeff_genome=GRCh37.75



#################################
## calling copy number variants##
#################################

#find any data with final.bam and excute by using samtools and bcftools
#sample_label: for the list of sample: basename, cut -f2 -d'. mean cut the second fragment after . of the sample name and run the second sample
#exemple for sample.9_A02.final.bam cut 9_A02 and run with the second sample 10_B02
mkdir -p $workdir
cd $workdir

#awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' $agilent > $agilent_corrected

qualimap bamqc -bam $sample -c \
	--java-mem-size=$gMUM \
	-gd HUMAN - hg19 \
	-gff $agilent_corrected \
	-hm 3 \
	-nr 1000 \
	-nt $threads \
	-nw 400 \
	-oc $workdir/$label.genome_coverage.txt \
	-os $workdir/$label.stats.txt \
	-outdir $workdir \
	-outfile $workdir/$label.report.pdf \
	-outformat HTML \
	-p non-strand-specific


done
