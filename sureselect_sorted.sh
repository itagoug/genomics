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
raw=$scratch/sickkids1/DATA/custom
bam_clean=$scratch/sickkids1/bam_clean/dup_clean
i2file=$scratch/sickkids1/DATA/custom
label=$raw/217_P7

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
LocatIt=$scratch/debugging/AGeNT



## open bam archive
#for the gz file do this
#for var in R1 R3; do

#gzip -d ${label}_$var.fastq.gz

#done


## RUN analysis
############
## PART I ##
############
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


##runing multiple samples


for i in $(_lb *.fast.gz | rev | cut -c 13- | rev | uniq)
do
tophat -o /path/to/_output/${i}
done





## align reads
#samtools workflow
#use the genome reference and the reads to generate a sam file
bwa mem -R '@RG\tID:A02\tSM:itagoug\tLB:library1' \
	${human_ref_genome} \
	${label}_R1.fastq \
	${label}_R3.fastq \
	> $_output.sam

	##remove duplicate using molecular barecode
	java -Xmx$gMUM -jar $LocatIt/LocatIt_v4.0.1.jar \
		-X $workdir -t $workdir -N 3000000 \
		-PM:xm,Q:xq,q:nQ,r:nR,I:ni \
		-m 1 -U -IS -OB -i -r -c 2500 \
		-l $agilent -o $_output.fixmate.bam $_output.sam ${label}_R2.fastq.gz


	## sorting Bam files
	samtools sort \
	  -O bam \
		-o $_output.sorted.bam \
	  -T $workdir/tmp \
		 $_output.fixmate.bam

	#make a files with any sample end with sorted.bam,named _listOfBams.
	#in this list look for any sample not end with .bai and add an index
	 _listOfBams=$(find $workdir -name "*sorted.bam")

	for _bam in $_listOfBams; do
		if [ ! -e $_bam.bai ]; then

				 samtools index -@$threads $_bam

		fi
	done



	#############
	## PART II ##
	#############
	## improvement
	## reduce the number of miscalls

	module load gatk/3.8
	java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ${human_ref_genome} \
		-I $_output.sorted.bam \
		-known $mills \
		-o $_output.intervals

	java -Xmx$gMUM -jar $gatk_dir/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${human_ref_genome} \
		-targetIntervals $_output.intervals \
		-I $_output.sorted.bam \
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
