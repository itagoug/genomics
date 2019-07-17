#!/bin/bash
#SBATCH -t 03:00:00
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=61440M
#SBATCH -J varscan
#SBATCH -p himem


set -eux
module load java/8 bwa/0.7.15 samtools/1.9
module load varscan/2.4.2 snpEff/4.3
module load tabix/0.2.6

## H4H directories
scratch=/cluster/projects/tiedemannlab/itagoug
raw=$scratch/sickkids1/bam_clean/dup_clean
#label=$raw/9_A02


## __DO NOT CHANGE__
#pbs=$(echo $PBS_JOBID | cut -f1 -d '.')

pbs=$SLURM_JOB_ID

bamfiles=$scratch/sickkids1/bamfiles
home=/cluster/home/itagoug
user_databases=$scratch/databases
workdir=$scratch/sickkids1/analysis/$pbs
#GRCh37-lite.fa link http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/
human_ref_genome=$scratch/human/lite/GRCh37-lite.fa
bam_clean=$scratch/sickkids1/bam_clean
mills=$scratch/human/gold/Mills.gold.standard.adapted_noChr.vcf
agilent=$scratch/human/agilent_original/74_clean.bed
snpeff_config=$home/scripts/snpEff.config

## parameters
pipeline=varscan2
gMUM=30g
threads=25
snpeff_genome=GRCh37.75



#################################
## calling copy number variants##
#################################

mkdir -p $workdir
cd $workdir


# __change __
_individual=167_P6
_sa=167_P6
_sb=167_NKT

# __do not change __
sample[1]=$raw/$_sa*bam
sample[2]=$raw/$_sb*bam


module load varscan/2.4.2



			output[1]=output.1.$pipeline.samtools.mpileup.agilent.ind_$_individual.$pbs

			module load samtools/1.3.1

			samtools mpileup -B --ignore-RG -l $agilent -f $human_ref_genome ${sample[1]} ${sample[2]} \
			| java -jar $varscan_dir/VarScan.jar mpileup2snp - \
			--min-var-freq 0.0000001 \
			--p-value 0.5 \
			--output-vcf 1 \
					> ${output[1]}.comparison.vcf

output[2]=output.2.$pipeline.snpeff_annot.$snpeff_genome.clinvar.dbsnp.ind_$_individual.$pbs
java -jar $snpeff_dir/SnpSift.jar annotate \
-noLog -noDownload \
-c $snpeff_config \
-clinvar \
${output[1]}.comparison.vcf \
| java -jar $snpeff_dir/SnpSift.jar annotate \
		 -noLog -noDownload \
		 -c $snpeff_config \
		 -dbsnp - \
| java -jar $snpeff_dir/snpEff.jar eff \
		 -noDownload -noLog \
		 -c $snpeff_config \
		 -stats ${output[2]}.comparison.html \
		 $snpeff_genome \
		 > ${output[2]}.comparison.vcf

sleep 100


grep -v "^#" ${output[2]}.comparison.vcf | grep -v "rs" | grep -iv "^[a-z]"|  awk 'match($3, /\./) {print $0}' | sed -e 's/WT.*ANN=.|//g' -e 's/|ENS.*ADR\t/\t/g' | cut -f 1-2,4-5,8-10 | sed -e 's/;/\t/g' -e 's/|/\t/g' -e's/:/\t/g' -e 's/\tADP=/\t/g' -e 's/%//g' | cut -f 1-8,12-16,26-30 | awk 'NF>17' | grep -v "=" | egrep -iv "^\.|muc|usp|rp11" > ${output[2]}.var.txt




_summary=$home/variants

cd $workdir
cp *.var.txt $_summary/varscan
cp *.stats.txt $_summary/stats_sickkids
cp *.genes.txt $_summary/genes_sickkids

for i in 1 2;
do
	gzip ${output[$i]}.comparison.vcf
done

echo -e "\n\nAnalyses are done, without complaints."
