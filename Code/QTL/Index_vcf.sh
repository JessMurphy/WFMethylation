#!/bin/bash

# copy vcf files
# cp /nfs/storage/math/gross-s2/projects/guatemala/MEGA_backbone_Full_Project_08012018_Deliverable/Genetics/Ian_QC/imputation/filtered/final/jessica_qtl/*.vcf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/vcf_files/

# set working directory
cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/vcf_files/

# sort and index the files
for i in {1..22..1}
do

bcftools sort -o Chr${i}_sorted.vcf WFguat_chr${i}.vcf

bgzip Chr${i}_sorted.vcf && tabix -p vcf Chr${i}_sorted.vcf.gz

done
