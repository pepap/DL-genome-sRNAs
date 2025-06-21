#!/bin/bash

metamodule add fastx-0.0.14

FAMRGZ=(
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLfoot1_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLfoot2_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLfoot3_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLhead1_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLhead2_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLhead3_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLjuv1_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLjuv2_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLjuv3_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLovo1_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLovo2_L1-all.re-collapsed.fa.gz
 /storage/brno12-cerit/home/pepap/brno1/TobiasBer/02.smallRNAseq-DL-newKit-20240504/FASTQ/CUTADAPT/DLovo3_L1-all.re-collapsed.fa.gz
)

zcat ${FAMRGZ[@]} | fastx_collapser -v -i - -o all-DL-smallRNAseq-merged.re-collapsed.fa

gzip -v all-DL-smallRNAseq-merged.re-collapsed.fa

