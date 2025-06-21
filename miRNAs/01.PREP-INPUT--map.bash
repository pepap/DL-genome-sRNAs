#!/bin/bash

metamodule add fastx-0.0.14

FAMRGZ=(
 INPUT_READ_FILE1.fa.gz
 INPUT_READ_FILE2.fq.gz
 ...
)

zcat ${FAMRGZ[@]} | fastx_collapser -v -i - -o all-DL-smallRNAseq-merged.re-collapsed.fa

gzip -v all-DL-smallRNAseq-merged.re-collapsed.fa

metamodule load star-2.7.7a
metamodule add samtools-1.9

exeCMD="STAR"
REF="/PATH/TO/INPUT/GENOMIC/FASTA"
IND="/PATH/TO/STAR/INDEX"
STOR="all-DL-smallRNAseq-merged.re-collapsed.fa.gz"
OUTBAM="all-DL"

${exeCMD} --runMode alignReads \
          --runThreadN  4 \
          --genomeDir ${IND} \
          --genomeLoad LoadAndRemove \
          --readFilesIn ${ICOLFAGZ} \
          --readFilesCommand zcat \
          --limitBAMsortRAM  20000000000 \
          --outFileNamePrefix ${OUTBAM}.se. \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 99999 \
          --outFilterMismatchNoverLmax 0.1 \
          --outFilterMatchNminOverLread 0.66 \
          --alignSJoverhangMin   999 \
          --alignSJDBoverhangMin 999 

samtools index ${OUTBAM}.se.Aligned.sortedByCoord.out.bam
