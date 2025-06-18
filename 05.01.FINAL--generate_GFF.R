#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.smrna.fa"
MAP="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.map.txt"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.EDITED.fasta"
SRS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/SINGLEREF"
EXE="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/bin/mireap.pl"

XCHANGES="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/04.xPARSE/xCHANGES"

RNDN=${RANDOM}
PAT1=(`awk '{ print $1 }' ${XCHANGES}`)
PAT2=(`awk '{ print $3 }' ${XCHANGES}`)
PAT3=(`awk '{ print $4 }' ${XCHANGES}`)

cp ${MAP} tmp.all.map.${RNDN}.txt
for ((i=0;i<`cat ${XCHANGES} | wc -l`;i++))
do
 awk -v "p2=${PAT2[$i]}" -v "p3=${PAT3[$i]}" '{ if( ($3==p2) && ($4==p3) ) { print "#"$0 } else { print $0 } }' tmp.all.map.${RNDN}.txt > all.map.corrected.txt
 mv all.map.corrected.txt tmp.all.map.${RNDN}.txt
done

cat ${XCHANGES}                    > _tmpCOMP
grep "^#" tmp.all.map.${RNDN}.txt >> _tmpCOMP
cat _tmpCOMP | sort -k2,4
rm -f _tmpCOMP
mv tmp.all.map.${RNDN}.txt all.map.corrected.txt

${EXE} -i ${FAS} -m all.map.corrected.txt -r ${REF} -o ODIR.SUCC -t DerLae

