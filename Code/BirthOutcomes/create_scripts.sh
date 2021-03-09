

for i in 9 13 14 15 18 20 21 22
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_pheno.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_pheno.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_pheno.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_pheno.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_pheno.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_pheno.R


done

# break the following chromosomes in half

for i in 3 4 5 6 8 10 11 12 16 17 19
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoA.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoA.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoA.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoA.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoA.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoA.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoB.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoB.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoB.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoB.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoB.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoB.R


done

# break the following chromosomes into thirds

for i in 1 2 7
do

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoAA.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoAA.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoAA.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoAA.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoAA.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoAA.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoBB.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoBB.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoBB.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoBB.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoBB.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoBB.R

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoCC.sh /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoCC.sh

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoCC.sh

cp /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr_phenoCC.R /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoCC.R

sed -i "s/mychrom/$i/g" /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/Chr${i}_phenoCC.R


done