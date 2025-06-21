#!/bin/bash

metamodule add gsl
metamodule add viennaRNA-2.4.17

export PERL5LIB="/PATH/TO/MIREAP/LIBRARY/mireap_0.2/lib:${PERL5LIB}"

FAS="/PATH/TO/all.smrna.fa"
MAP="/PATH/TO/all.map.txt"
REF="/PATH/TO/INPUT/GENOMIC/FASTA"
EXE="/PATH/TO/MIREAP/mireap_0.2/bin/mireap.pl"
