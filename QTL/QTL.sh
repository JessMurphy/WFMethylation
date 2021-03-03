#!/bin/bash

# set working directory
cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/

# run QTLTools for mychrom
/nfs/apps/math/links/QTLtools cis --vcf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/vcf_files/Chrmychrom_sorted.vcf.gz --bed /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/bed_files/Chrmychrom_sorted.bed.gz --cov /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/QTL_covs.txt --permute 1000 --out /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chrmychrom_permutations.txt --seed 12345