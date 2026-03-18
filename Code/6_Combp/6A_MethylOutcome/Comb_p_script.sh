######################

### First create the file 
#### chrom start end
##### make the end 1 more than the start
##### Turn off scientific notation: options(scipen=999)


#### first column must be: (with the #)
#### Just added manually
#chrom	start	end	rawp

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep

nano ToSort_in...   

#### sort with bedtools. MUST sort with bedtools, even if its already sorted

#### the bedtools program runs from the home directory
cd

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ToSort_interaction4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/Sorted_interaction4combp.bed 

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ToSort_arm4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/Sorted_arm4combp.bed 

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ToSort_bmi4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/Sorted_bmi4combp.bed 

#### comb-p requires the chromosome to be 'chr1' 
##### Fix in R (same script)
options(scipen=999) ## turn off scientific notation


#####
# Run Comb-p
#####

#### navigate the the correct folder
cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS

comb-p pipeline -c 4 --seed 0.05 --dist 300 -p results/combp_interaction_results --region-filter-p 0.05 --anno hg19 FilePrep/Sorted_interaction4combp.bed

#### note that the command 'bedtools' does not pull up bedtools here, and that is why the annotation doesn't work!
#### The results will write and then there will be a bunch of errors
#### The errors are related to bedtools/annotating the results




comb-p pipeline -c 4 --seed 0.05 --dist 300 -p results/combp_arm_results --region-filter-p 0.05 --anno hg19 FilePrep/Sorted_arm4combp.bed


comb-p pipeline -c 4 --seed 0.05 --dist 300 -p results/combp_bmi_results --region-filter-p 0.05 --anno hg19 FilePrep/Sorted_bmi4combp.bed




### try to make a Manhattan plot
#### works on my computer

comb-p manhattan InteractionPvalues4combp_sorted.bed -c 4 \
    --image interaction.manhattan.png



###############
## ME's
###############

### Add the #

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep

nano ...

#### the bedtools program runs from the home directory
cd

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_ToSort_interaction4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_Sorted_interaction4combp.bed 

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_ToSort_arm4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_Sorted_arm4combp.bed 

bedtools2/bin/bedtools sort -i /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_ToSort_bmi4combp.bed > /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep/ME_Sorted_bmi4combp.bed 

#### comb-p requires the chromosome to be 'chr1' 
##### Fix in R (same script)
options(scipen=999) ## turn off scientific notation


#####
# Run Comb-p
#####

#### navigate the the correct folder
cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS

comb-p pipeline -c 4 --seed 0.05 --dist 300 -p ME_results/ME_combp_interaction_results --region-filter-p 0.05 --anno hg19 FilePrep/ME_Sorted_interaction4combp.bed

#### note that the command 'bedtools' does not pull up bedtools here, and that is why the annotation doesn't work!
#### The results will write and then there will be a bunch of errors
#### The errors are related to bedtools/annotating the results




comb-p pipeline -c 4 --seed 0.05 --dist 300 -p ME_results/ME_combp_arm_results --region-filter-p 0.05 --anno hg19 FilePrep/ME_Sorted_arm4combp.bed


comb-p pipeline -c 4 --seed 0.05 --dist 300 -p ME_results/ME_combp_bmi_results --region-filter-p 0.05 --anno hg19 FilePrep/ME_Sorted_bmi4combp.bed


