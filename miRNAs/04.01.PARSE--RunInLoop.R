#!/bin/bash

#01#CHR="scaffold.1"
#02#CHR="scaffold.12"
#03#CHR="scaffold.13"
#04#CHR="scaffold.26"
CHR="scaffold.99"
STR="-"
#!#STR="+"

READ=1
TIER=1
while [[ ${READ} == 1 ]]
do echo ""; echo " ** tier : $TIER **"; echo ""

 bash parserun_mireap ${CHR} ${STR}
 cat parse.map.iSTR.txt >> xCHANGES

 PAT1=`awk '{ print $1 }' parse.map.iSTR.txt`
 PAT2=`awk '{ print $3 }' parse.map.iSTR.txt`
 PAT3=`awk '{ print $4 }' parse.map.iSTR.txt`
#!# awk -v "p1=${PAT1}" -v "p2=${PAT2}" -v "p3=${PAT3}" '{ if( ($1==p1) && ($3==p2) && ($4==p3) ) { print "#"$0 } else { print $0 } }' tmp.map.${CHR}.iSTR.txt > XXX
 awk -v "p2=${PAT2}" -v "p3=${PAT3}" '{ if( ($3==p2) && ($4==p3) ) { print "#"$0 } else { print $0 } }' tmp.map.${CHR}.iSTR.txt > XXX
 mv XXX tmp.map.${CHR}.iSTR.txt

 bash checkRun tmp.map.${CHR}.iSTR.txt

 if [[ $? != 0 ]]
 then
  rm -rf ODIR.* parse.map.iSTR.txt*
 else
  READ=0
 fi

 let TIER=TIER+1

done

