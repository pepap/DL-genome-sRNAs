#!/bin/bash

metamodule add samtools-1.9

#Test if bam file is coordinate-sorted using samtools package
# @pepap : input is a single BAM file
# @pepap : if OK, returns "SO:coordinate"
samtools view -H $1 | \
grep '^@HD' | \
awk -F'\t' '{
 for(i=1;i<=NF;i++) if($i ~ /^SO:/) print $i
            }'

