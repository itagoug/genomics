## Genome analysis

Scripts found in this repo will start analyses on __genomic data__ (fasta, fastq, bam, bed). Most tools used are well referenced in the field and used in bioinformatics in general.

Genome analysis software will:
- align and convert FASQT to bam.
- remove duplicate reads for SureSelect custom capture.
- check coverage.
- Variant calling (SNVs, indels).
- Germline filtering analysis.
-Gene and variant annotation

## Testing

Runnings tests were performed on PBS and Slurm on H4H cluster at UHN, Toronto Canada.

#strelka2:
Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs.
# https://github.com/Illumina/strelka
#Varscan2
