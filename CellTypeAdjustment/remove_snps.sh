#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/

# find overlap between snps and CpGs
fgrep -xf chrpos_clean.txt CpGs_noX.txt > removed_snps.txt

# remove overlapped snps from methylation matrix
fgrep -wvf removed_snps.txt methylation_filtered_noX.txt > methylation_noX_noSNPs.txt



