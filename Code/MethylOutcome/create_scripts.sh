

for i in 9 13 14 15 18 20 21 22
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysis.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysis.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysis.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysis.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysis.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysis.R


done

# break the following chromosomes in half

for i in 3 4 5 6 8 10 11 12 16 17 19
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisA.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisA.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisA.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisA.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisA.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisA.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisB.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisB.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisB.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisB.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisB.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisB.R


done

# break the following chromosomes into thirds

for i in 1 2 7
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisAA.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/birth_outcomes/analysis/Chr${i}_analysisAA.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisAA.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisAA.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/birth_outcomes/analysis/Chr${i}_analysisAA.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisAA.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisBB.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/birth_outcomes/analysis/Chr${i}_analysisBB.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisBB.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisBB.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/birth_outcomes/analysis/Chr${i}_analysisBB.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisBB.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisCC.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/birth_outcomes/analysis/Chr${i}_analysisCC.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisCC.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr_analysisCC.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/birth_outcomes/Chr${i}_analysisCC.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/analysis/Chr${i}_analysisCC.R


done