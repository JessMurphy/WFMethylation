#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/

for i in {1..22..1}
do

cut -d " " -f 1 < Chr${i}_significant.txt > Chr${i}_QTL.txt

done
