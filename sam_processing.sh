#!/bin/bash
#PBS -l nodes=1:ppn=25,walltime=24:00:00,vmem=30g,mem=220g
#PBS -N bamProc
#PBS -q all
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M inestagoug@gmail.com


set -eux

module load java/8 bwa/0.7.15 samtools/1.9
module load bamtools/2.4.2 snpEff/4.3
module load tabix/0.2.6

## __CHANGE__
scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/analysis/905365
bam_clean=$scratch/sickkids1/bam_clean/dup_clean
label=$raw/217_P6

## PARAMETERS
pipeline=bcftools
gMUM=30g
threads=25
snpeff_genome=GRCh37.75


## __DO NOT CHANGE__
#raw=$scratch/sickkids1/testing/raw.samples
pbs=$(echo $PBS_JOBID | cut -f1 -d '.')
bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
_lb=$(basename $label)
_output=$workdir/$_lb

#human_ref_genome=$scratch/human/hg19/GCF_000001405.25_GRCh37.p13_genomic.fna change it to GRCh37-lite
#GRCh37-lite.fa link http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
mills=$scratch/human/gold/Mills.gold.standard.adapted_noChr.vcf
agilent=$scratch/human/agilent_original/206.clean.bed
snpeff_config=$home/scripts/snpEff.config



mkdir -p $workdir
cd $workdir

# if the ref.fai does not existe then create the sequence dictionary
if [ ! -e $human_ref_genome.fai ]; then
		## index files
		bwa index ${human_ref_genome}

		module load gatk/4.0.5.1
		java -Xmx$gMUM -jar $gatk_dir/gatk-package-4.0.5.1-local.jar \
			CreateSequenceDictionary \
			-R ${human_ref_genome}

		samtools faidx ${human_ref_genome}

fi





#############
## PART II ##
#############
## improvement
## reduce the number of miscalls
module load gatk/3.8
java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${human_ref_genome} \
	-I $raw/217_P6.sorted.bam \
	-known $mills \
	-o $_output.intervals

java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${human_ref_genome} \
	-targetIntervals $_output.intervals \
	-I $raw/217_P6.sorted.bam \
	-known $mills \
	-o $_output.realigned.bam


	## remove duplicate reads
	module load picard/2.10.9
	java -Xmx$gMUM -jar $picard_dir/picard.jar MarkDuplicates \
		VALIDATION_STRINGENCY=LENIENT \
		I=$_output.realigned.bam \
		O=$_output.rmdup.bam \
		M=$_output.rmdup.metrics


		##############
		## PART III ##
		##############

		_listOfBams=$(find $workdir -name "*rmdup.bam")


		for _bam in $_listOfBams; do
		 if [ ! -e $_bam.bai ]; then

					samtools index -@$threads $_bam

		 fi
		done

	## realign indels
	java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ${human_ref_genome} \
		-I $_output.rmdup.bam \
		-known $mills \
		-o $_output.rmdup.intervals

	java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${human_ref_genome} \
		-I $_output.rmdup.bam \
		-targetIntervals $_output.rmdup.intervals \
		-known $mills \
		-o $_output.final.bam

	#index bam
	#make a list with anyfile end witg final.bam and index it
		_listOfFinal=$(find $workdir -name "*.final.bam")

		for _bam in $_listOfFinal; do
		 if [ ! -e $_bam.bai ]; then

					samtools index -@$threads $_bam

		 fi
		done


	#move the final file to bam_clean folder
	#we will excute the others tools such as bcftools or varscan for the analysis

	mv ${_output}.final* $bam_clean


	echo -e "\n\nAnalyses are done, without complaints."
