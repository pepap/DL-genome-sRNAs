#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.smrna.fa"
MAP="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/BAM.smallRNAseq/00.stats/04.mireap/02.PREP/all.map.txt"
REF="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/Deroceras_laeve_CZ.1.0.fasta"
SRS="/storage/brno12-cerit/home/pepap/brno1/TobiasBer/10.DLgenome-CZ-1.0/SINGLEREF"
EXE="/storage/brno12-cerit/home/pepap/brno1/pLIBRARY/mireap_0.2/bin/mireap.pl"

#>> @pepap : running all at once (not used)
#!#${EXE} -i ${FAS} -m ${MAP} -r ${REF} -o ODIR.1 -t DLtest

awk '{ count[$2]++ } END{ for(i in count) { print i, count[i] } }' ${MAP} | sort -n -k1,1 > tmp.srefs.list
#awk '{ count[$2]++ } END{ for(i in count) { print i, count[i] } }' ${MAP} | sort -n -k1,1 | grep -w "HiC.scaffold.11" > tmp.srefs.list

mkdir mFWD mREV

#for iSRS in `cat tmp.srefs.list | head -n 1 | awk '{ print $1 }'`
for iSRS in `cat tmp.srefs.list | awk '{ print $1 }'`
do echo ${iSRS}

 grep -w ${iSRS} ${MAP} | awk '{ if( $NF=="+" ) { print $0 } }' > mFWD/tmp.map.${iSRS}.fwd.txt
 grep -w ${iSRS} ${MAP} | awk '{ if( $NF=="-" ) { print $0 } }' > mREV/tmp.map.${iSRS}.rev.txt

##>> @pepap : for these selected reads the prediction fails with unknown technical error
#wc -l mFWD/tmp.map.${iSRS}.fwd.txt
#wc -l mREV/tmp.map.${iSRS}.rev.txt
# if [[ $iSRS == "scaffold.11"  ]]
# then
#  sed -i -e "31,36d;41,49d;68,69d"                                                                                           mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.13" ]]
# then
#  sed -i -e "22d;26d;28d;534,535d;540,541d"                                                                                  mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.16" ]]
# then
#  sed -i -e "259,260d;336d"                                                                                                  mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.17" ]]
# then
#  sed -i -e "182,206d;208,232d;805d;808d;819,821d"                                                                           mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.19" ]]
# then
#  sed -i -e "27,28d"                                                                                                         mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.27" ]]
# then
#  sed -i -e "46,57d;59,62d;246,248d;250d"                                                                                    mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.7" ]]
# then
#  sed -i -e "531,536d"                                                                                                       mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.75" ]]
# then
#  sed -i -e "1,2d"                                                                                                           mREV/tmp.map.${iSRS}.rev.txt
# fi
# if [[ $iSRS == "scaffold.8" ]]
# then
#  sed -i -e "308,313d;318,326d"                                                                                              mREV/tmp.map.${iSRS}.rev.txt
#  sed -i -e "58d;60d;64d;72,73d;75,77d;80,81d;88,90d;93,95d;98,101d;181d;332,335d;338,339d;348,355d;617d;622d;627d;636,637d" mFWD/tmp.map.${iSRS}.fwd.txt
# fi
#wc -l mFWD/tmp.map.${iSRS}.fwd.txt
#wc -l mREV/tmp.map.${iSRS}.rev.txt
##<<

 echo " > FWD"
 rm -f ODIR.TEST/mireap-${iSRS}.fwd.*
 ${EXE} -i ${FAS} -m mFWD/tmp.map.${iSRS}.fwd.txt -r ${SRS}/${iSRS}.fa -o ODIR.TEST -t ${iSRS}.fwd
 echo " > REV"
 rm -f ODIR.TEST/mireap-${iSRS}.rev.*
 ${EXE} -i ${FAS} -m mREV/tmp.map.${iSRS}.rev.txt -r ${SRS}/${iSRS}.fa -o ODIR.TEST -t ${iSRS}.rev

# rm -f tmp.map.${iSRS}.fwd.txt tmp.map.${iSRS}.rev.txt

done

