#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/bed_files/

for i in {1..21..1}
do

perl -p -i -e 's/ /\t/g' Chr${i}.bed

(head -n 1 Chr${i}.bed && tail -n +2 Chr${i}.bed | sort -n -k 2) > Chr${i}_sorted.bed

bgzip Chr${i}_sorted.bed && tabix -p bed Chr${i}_sorted.bed.gz

done