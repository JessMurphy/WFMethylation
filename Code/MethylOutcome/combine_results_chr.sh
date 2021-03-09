#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/

# for the chromosomes broken in half
for i in 3 4 5 6 8 10 11 12 16 17 19
do

# remove the header of the second files
sed -i '1d' Chr${i}_resultsB.txt

# combine the two chromosome files into one
cat Chr${i}_resultsA.txt Chr${i}_resultsB.txt > Chr${i}_results.txt

done

# for the chromosomes broken into thirds
for i in 1 2 7
do

# remove the header of the second file
sed -i '1d' Chr${i}_resultsBB.txt

# remove the header of the third file
sed -i '1d' Chr${i}_resultsCC.txt

# combine the three chromosome files into one
cat Chr${i}_resultsAA.txt Chr${i}_resultsBB.txt Chr${i}_resultsCC.txt > Chr${i}_results.txt

done
