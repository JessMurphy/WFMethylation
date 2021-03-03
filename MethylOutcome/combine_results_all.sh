
cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/

# remove the header of every file except the first one
for i in {2..22..1}
do

sed -i '1d' Chr${i}_filtered.txt

done

# combine chromosome results into one file
cat Chr1_filtered.txt Chr2_filtered.txt Chr3_filtered.txt Chr4_filtered.txt Chr5_filtered.txt Chr6_filtered.txt Chr7_filtered.txt Chr8_filtered.txt Chr9_filtered.txt Chr10_filtered.txt Chr11_filtered.txt Chr12_filtered.txt Chr13_filtered.txt Chr14_filtered.txt Chr15_filtered.txt Chr16_filtered.txt Chr17_filtered.txt Chr18_filtered.txt Chr19_filtered.txt Chr20_filtered.txt Chr21_filtered.txt Chr22_filtered.txt > EWAS_results.txt