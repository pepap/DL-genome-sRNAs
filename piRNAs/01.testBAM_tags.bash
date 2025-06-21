#!/bin/bash

metamodule add samtools-1.9

#Test if bam file has NM and NH tags using samtools package
# @pepap : input is a single BAM file
# @pepap : if OK, returns "NM and NH tags found"
samtools view $1 | \
head -n 1000 | \
awk '{
      if($0 ~ /NM:i:/ && $0 ~ /NH:i:/) print "NM and NH tags found"
 else if($0 ~ /NM:i:/)                 print "!! Only NM tag found !!"
 else if($0 ~ /NH:i:/)                 print "!! Only NH tag found !!"
 else                                  print "!! Neither NM nor NH tag found !!"
 exit
     }'

