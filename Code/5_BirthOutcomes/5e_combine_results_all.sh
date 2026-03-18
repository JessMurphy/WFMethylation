#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/

for i in {2..22..1}
do

# remove the header of the second files
sed -i '1d' Chr${i}_filtered_LGAZ.txt
sed -i '1d' Chr${i}_filtered_WGAZ.txt
sed -i '1d' Chr${i}_filtered_HCGAZ.txt
sed -i '1d' Chr${i}_filtered_WLGAZ.txt

done

cat Chr1_filtered_LGAZ.txt Chr2_filtered_LGAZ.txt Chr3_filtered_LGAZ.txt Chr4_filtered_LGAZ.txt Chr5_filtered_LGAZ.txt Chr6_filtered_LGAZ.txt Chr7_filtered_LGAZ.txt Chr8_filtered_LGAZ.txt Chr9_filtered_LGAZ.txt Chr10_filtered_LGAZ.txt Chr11_filtered_LGAZ.txt Chr12_filtered_LGAZ.txt Chr13_filtered_LGAZ.txt Chr14_filtered_LGAZ.txt Chr15_filtered_LGAZ.txt Chr16_filtered_LGAZ.txt Chr17_filtered_LGAZ.txt Chr18_filtered_LGAZ.txt Chr19_filtered_LGAZ.txt Chr20_filtered_LGAZ.txt Chr21_filtered_LGAZ.txt Chr22_filtered_LGAZ.txt > LGAZ_results.txt

cat Chr1_filtered_WGAZ.txt Chr2_filtered_WGAZ.txt Chr3_filtered_WGAZ.txt Chr4_filtered_WGAZ.txt Chr5_filtered_WGAZ.txt Chr6_filtered_WGAZ.txt Chr7_filtered_WGAZ.txt Chr8_filtered_WGAZ.txt Chr9_filtered_WGAZ.txt Chr10_filtered_WGAZ.txt Chr11_filtered_WGAZ.txt Chr12_filtered_WGAZ.txt Chr13_filtered_WGAZ.txt Chr14_filtered_WGAZ.txt Chr15_filtered_WGAZ.txt Chr16_filtered_WGAZ.txt Chr17_filtered_WGAZ.txt Chr18_filtered_WGAZ.txt Chr19_filtered_WGAZ.txt Chr20_filtered_WGAZ.txt Chr21_filtered_WGAZ.txt Chr22_filtered_WGAZ.txt > WGAZ_results.txt

cat Chr1_filtered_HCGAZ.txt Chr2_filtered_HCGAZ.txt Chr3_filtered_HCGAZ.txt Chr4_filtered_HCGAZ.txt Chr5_filtered_HCGAZ.txt Chr6_filtered_HCGAZ.txt Chr7_filtered_HCGAZ.txt Chr8_filtered_HCGAZ.txt Chr9_filtered_HCGAZ.txt Chr10_filtered_HCGAZ.txt Chr11_filtered_HCGAZ.txt Chr12_filtered_HCGAZ.txt Chr13_filtered_HCGAZ.txt Chr14_filtered_HCGAZ.txt Chr15_filtered_HCGAZ.txt Chr16_filtered_HCGAZ.txt Chr17_filtered_HCGAZ.txt Chr18_filtered_HCGAZ.txt Chr19_filtered_HCGAZ.txt Chr20_filtered_HCGAZ.txt Chr21_filtered_HCGAZ.txt Chr22_filtered_HCGAZ.txt > HCGAZ_results.txt

cat Chr1_filtered_WLGAZ.txt Chr2_filtered_WLGAZ.txt Chr3_filtered_WLGAZ.txt Chr4_filtered_WLGAZ.txt Chr5_filtered_WLGAZ.txt Chr6_filtered_WLGAZ.txt Chr7_filtered_WLGAZ.txt Chr8_filtered_WLGAZ.txt Chr9_filtered_WLGAZ.txt Chr10_filtered_WLGAZ.txt Chr11_filtered_WLGAZ.txt Chr12_filtered_WLGAZ.txt Chr13_filtered_WLGAZ.txt Chr14_filtered_WLGAZ.txt Chr15_filtered_WLGAZ.txt Chr16_filtered_WLGAZ.txt Chr17_filtered_WLGAZ.txt Chr18_filtered_WLGAZ.txt Chr19_filtered_WLGAZ.txt Chr20_filtered_WLGAZ.txt Chr21_filtered_WLGAZ.txt Chr22_filtered_WLGAZ.txt > WLGAZ_results.txt