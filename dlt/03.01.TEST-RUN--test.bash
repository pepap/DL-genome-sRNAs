#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.smrna.fa"
MAP="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.map.txt"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.fasta"
SRS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/SINGLEREF"
EXE="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/bin/mireap.pl"

#iSRS="scaffold.16"
#iSRS="scaffold.17"
#iSRS="scaffold.19"
#iSRS="scaffold.27"
#iSRS="scaffold.7"
#iSRS="scaffold.75"
iSRS="scaffold.8"
echo ${iSRS}
awk '{ count[$2]++ } END{ for(i in count) { print i, count[i] } }' ${MAP} | sort -n -k1,1 | grep -w "${iSRS}" > tmp.srefs.list

#T# grep -w ${iSRS} ${MAP} | awk '{ if( $NF=="+" ) { print $0 } }' > mFWD/tmp.map.${iSRS}.fwd.txt
 grep -w ${iSRS} ${MAP} | awk '{ if( $NF=="-" ) { print $0 } }' > mREV/tmp.map.${iSRS}.rev.txt

#>> @pepap : for these selected reads the prediction fails with unknown technical error
 if [[ $iSRS == "scaffold.11" ]]
 then
  sed -i -e "31,36d;41,49d;68,69d" mREV/tmp.map.${iSRS}.rev.txt
 fi
 if [[ $iSRS == "scaffold.13" ]]
 then
  sed -i -e "22d;26d;28d;534,535d;540,541d" mREV/tmp.map.${iSRS}.rev.txt
 fi
 if [[ $iSRS == "scaffold.1x" ]]
 then
  sed -i -e "22d;26d;28d;534,535d;540,541d" mREV/tmp.map.${iSRS}.rev.txt
 fi
#<<

#T# HTOT=`cat mREV/tmp.map.${iSRS}.rev.txt | wc -l`
 HTOT=`cat mFWD/tmp.map.${iSRS}.fwd.txt | wc -l`
 for ((h=1;h<=${HTOT};h++))
 do echo $h
#T# head -n $h mREV/tmp.map.${iSRS}.rev.txt | tail -n 1 > TMP.MAP.TXT
 head -n $h mFWD/tmp.map.${iSRS}.fwd.txt | tail -n 1 > TMP.MAP.TXT

 echo " > FWD"
 rm -f ODIR.TEST/mireap-${iSRS}.fwd.*
 ${EXE} -i ${FAS} -m TMP.MAP.TXT -r ${SRS}/${iSRS}.fa -o ODIR.TEST -t ${iSRS}.fwd
#T# echo " > REV"
#T# rm -f ODIR.TEST/mireap-${iSRS}.rev.*
#T# ${EXE} -i ${FAS} -m TMP.MAP.TXT -r ${SRS}/${iSRS}.fa -o ODIR.TEST -t ${iSRS}.rev

 done

# rm -f tmp.map.${iSRS}.fwd.txt tmp.map.${iSRS}.rev.txt

