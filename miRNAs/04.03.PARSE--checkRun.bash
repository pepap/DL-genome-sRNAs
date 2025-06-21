#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.smrna.fa"
MAP="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.map.txt"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.fasta"
SRS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/SINGLEREF"
EXE="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/bin/mireap.pl"
  
iMAP=$1
iSRS=`echo ${iMAP} | awk -F"." '{ printf( "%s.%s", $3, $4) }'`

rm -rf ODIR.checkRun
${EXE} -i ${FAS} -m ${iMAP} -r ${SRS}/${iSRS}.fa -o ODIR.checkRun -t ${iSRS}.iSTR

