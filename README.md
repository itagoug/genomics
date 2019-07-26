## Genome analysis

Scripts found in this repo will start analyses on __genomic data__ (fasta, fastq, bam, bed). Most tools used are well referenced in the field and used in bioinformatics in general.

Genome analysis software will:
- align and convert FASQT to bam (samtools).
- remove duplicate reads (agilent sureselect).
- check coverage (Qualimap).
- Variant calling (Varscan2, bcftools and strelka2).

## Testing

Runnings tests were performed on PBS and Slurm.
