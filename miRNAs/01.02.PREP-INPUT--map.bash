#!/bin/bash

metamodule load star-2.7.7a
metamodule add samtools-1.9

exeCMD="STAR"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.fasta"
IND="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/STAR--2.7.7a_index"
STOR="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/01.merge_input_reads"

FASTQS=(
 all-DL
)
OUTBAM=(
 all-DL
)
LIB_PE=(
 SINGLE
)
SUFF="-smallRNAseq-merged.re-collapsed.fa.gz"

i=0

${exeCMD} --runMode alignReads \
          --runThreadN  4 \
          --genomeDir ${IND} \
          --genomeLoad LoadAndRemove \
          --readFilesIn ${FASTQS[$i]}${SUFF} \
          --readFilesPrefix ${STOR}/ \
          --readFilesCommand zcat \
          --limitBAMsortRAM  20000000000 \
          --outFileNamePrefix ${OUTBAM[$i]}.se. \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 99999 \
          --outFilterMismatchNoverLmax 0.1 \
          --outFilterMatchNminOverLread 0.66 \
          --alignSJoverhangMin   999 \
          --alignSJDBoverhangMin 999 

samtools index ${OUTBAM[$i]}.se.Aligned.sortedByCoord.out.bam

metamodule add fastQC-0.11.5
fastqc ${OUTBAM[$i]}.se.Unmapped.out.mate1
gzip   ${OUTBAM[$i]}.se.Unmapped.out.mate1

