#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/

# for the chromosomes broken in half
for i in 3 4 5 6 8 10 11 12 16 17 19
do

# remove the header of the second files
sed -i '1d' Chr${i}_resultsB_LGAZ.txt
sed -i '1d' Chr${i}_resultsB_WGAZ.txt
sed -i '1d' Chr${i}_resultsB_HCGAZ.txt
sed -i '1d' Chr${i}_resultsB_WLGAZ.txt

# combine the two chromosome files into one
cat Chr${i}_resultsA_LGAZ.txt Chr${i}_resultsB_LGAZ.txt > Chr${i}_results_LGAZ.txt
cat Chr${i}_resultsA_WGAZ.txt Chr${i}_resultsB_WGAZ.txt > Chr${i}_results_WGAZ.txt
cat Chr${i}_resultsA_HCGAZ.txt Chr${i}_resultsB_HCGAZ.txt > Chr${i}_results_HCGAZ.txt
cat Chr${i}_resultsA_WLGAZ.txt Chr${i}_resultsB_WLGAZ.txt > Chr${i}_results_WLGAZ.txt

done

# for the chromosomes broken into thirds
for i in 1 2 7
do

# remove the header of the second file
sed -i '1d' Chr${i}_resultsBB_LGAZ.txt
sed -i '1d' Chr${i}_resultsBB_WGAZ.txt
sed -i '1d' Chr${i}_resultsBB_HCGAZ.txt
sed -i '1d' Chr${i}_resultsBB_WLGAZ.txt

# remove the header of the third file
sed -i '1d' Chr${i}_resultsCC_LGAZ.txt
sed -i '1d' Chr${i}_resultsCC_WGAZ.txt
sed -i '1d' Chr${i}_resultsCC_HCGAZ.txt
sed -i '1d' Chr${i}_resultsCC_WLGAZ.txt

# combine the three chromosome files into one
cat Chr${i}_resultsAA_LGAZ.txt Chr${i}_resultsBB_LGAZ.txt Chr${i}_resultsCC_LGAZ.txt > Chr${i}_results_LGAZ.txt
cat Chr${i}_resultsAA_WGAZ.txt Chr${i}_resultsBB_WGAZ.txt Chr${i}_resultsCC_WGAZ.txt > Chr${i}_results_WGAZ.txt
cat Chr${i}_resultsAA_HCGAZ.txt Chr${i}_resultsBB_HCGAZ.txt Chr${i}_resultsCC_HCGAZ.txt > Chr${i}_results_HCGAZ.txt
cat Chr${i}_resultsAA_WLGAZ.txt Chr${i}_resultsBB_WLGAZ.txt Chr${i}_resultsCC_WLGAZ.txt > Chr${i}_results_WLGAZ.txt

done
