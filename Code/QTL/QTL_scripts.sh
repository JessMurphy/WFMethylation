#!/bin/bash

# generate QTL scripts for each chromosome
for i in {1..22..1}
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/chr_scripts/QTL.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/chr_scripts/Chr${i}_QTL.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/chr_scripts/Chr${i}_QTL.sh

done