#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.smrna.fa"
MAP="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.map.txt"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.fasta"
SRS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/SINGLEREF"
EXE="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/bin/mireap.pl"

iSRS=$1
iSTR=$2

echo " => ${iSRS} : ${iSTR}"

if [[ ! -f tmp.map.${iSRS}.iSTR.txt ]]
then
 grep -w ${iSRS} ${MAP} | awk -v "iSTR=${iSTR}" '{ if( $NF==iSTR ) { print $0 } }' > tmp.map.${iSRS}.iSTR.txt
fi

TOTN=`cat tmp.map.${iSRS}.iSTR.txt | wc -l`

cp tmp.map.${iSRS}.iSTR.txt parse.map.iSTR.txt

LINE=0
ITER=1
while [[ ${LINE} == 0 ]]
do echo " iter : $ITER | number of lines = ${TOTN}"

NH=`echo $((${TOTN}/2))`
NT=`echo $((${TOTN}-${NH}))`

echo $NH $NT
cp parse.map.iSTR.txt parse.map.iSTR.txt.bcp

head -n ${NH} parse.map.iSTR.txt > parse.map.iSTR.txt.TMP
rm -f ODIR.iSTR/mireap-${iSRS}.iSTR.*
${EXE} -i ${FAS} -m parse.map.iSTR.txt.TMP -r ${SRS}/${iSRS}.fa -o ODIR.iSTR -t ${iSRS}.iSTR
if [[ $? != 0 ]]
then
mv parse.map.iSTR.txt.TMP parse.map.iSTR.txt
else
tail -n ${NT} parse.map.iSTR.txt > parse.map.iSTR.txt.TMP
rm -f ODIR.iSTR/mireap-${iSRS}.iSTR.*
${EXE} -i ${FAS} -m parse.map.iSTR.txt.TMP -r ${SRS}/${iSRS}.fa -o ODIR.iSTR -t ${iSRS}.iSTR
fi

if [[ $? != 0 ]]
then
mv parse.map.iSTR.txt.TMP parse.map.iSTR.txt
fi

TOTN=`cat parse.map.iSTR.txt | wc -l`

if [[ ${TOTN} == 1 ]]
then
LINE=1
fi

let ITER=ITER+1

done

FHEAD=`awk '{ print $1 }' parse.map.iSTR.txt`
echo "parse.map.iSTR.txt :"
cat parse.map.iSTR.txt
grep "^>${FHEAD} " -A 1 ${FAS}

##for ((N=1;N<=${TOTN};N++))
#for ((N=1;N<=2;N++))
#do echo $N
#
# head -n $N tmp.map.${iSRS}.iSTR.txt | tail -n 1 > tmp.inp.map.txt
# rm -f ODIR.iSTR/mireap-${iSRS}.iSTR.*
# ${EXE} -i ${FAS} -m tmp.inp.map.txt -r ${SRS}/${iSRS}.fa -o ODIR.iSTR -t ${iSRS}.iSTR
#
#done

