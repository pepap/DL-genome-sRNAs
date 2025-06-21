#!/bin/bash

R --no-save < 02.02.FORMAT-INPUT--bam2prep.R

sed -i "s/_/./g" all.map.txt

R --no-save < 02.03.FORMAT-INPUT--splitREF.R
