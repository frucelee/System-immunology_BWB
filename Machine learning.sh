#!/bin/bash
#
#SBATCH --job-name=machine learning
#SBATCH --output=MC_mpi.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

for j in $(seq 1 10);
do for i in $(seq 1 144);
do
plink --bfile /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/genotypes --keep resampling/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID --chr-set 29 --allow-extra-chr
plink --bfile ${j%}_GCTA_ID --indiv-sort f resampling/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID_1 --chr-set 29 --allow-extra-chr
gcta64 --mlma-loco --bfile ${j%}_GCTA_ID_1 --pheno resampling/${j%}_pheno_.txt --out mc_${i%}_loco_${j%} --thread-num 10 --qcovar resampling/${j%}_GCTA_qcov.txt --mpheno $i --covar resampling/${j%}_GCTA_cov.txt
awk '$9<=5E-05' mc_${i%}_loco_${j%}.loco.mlma >mc/sig/mc_${i%}_loco_${j%}_sig.subset.txt
#rm -r *bed *bim * fam * log *nosex *mlma
done
done

cd mc/sig
for j in $(seq 1 10);
do for i in $(seq 1 144);
do
awk '{print $2}' mc_${i%}_loco_${j%}_sig.subset.txt > SNPs.txt
vcftools --vcf geno_mind_test_vcf.vcf --snps SNPs.txt --recode --recode-INFO-all --out test1
plink --vcf test1.recode.vcf --keep-allele-order --double-id --make-bed --out test1_1 --chr-set 29 --allow-extra-chr --set-missing-var-ids @#
plink --bfile test1_1 --recodeA --out test1_1_geno --chr-set 29 --allow-extra-chr
head test1_1_geno.raw| cut -d" " -f 1-9
cat test1_1_geno.raw | cut -d" " -f2,7- | sed 's/_[A-Z]//g' >mc/mc_${i%}_loco_S${j%}_genotype.txt
mv mc_${i%}_loco_${j%}_sig.subset.txt mc/sig/sig
#rm -r *bed *bim * fam * log *nosex *mlma *vcf
done
done

for i in $(seq 1 144);
do
mv mc/mc_${i%}_loco_1_genotype.txt mc/sig/S1
mv mc/mc_${i%}_loco_2_genotype.txt mc/sig/S2
mv mc/mc_${i%}_loco_3_genotype.txt mc/sig/S3
mv mc/mc_${i%}_loco_4_genotype.txt mc/sig/S4
mv mc/mc_${i%}_loco_5_genotype.txt mc/sig/S5
done

##prediction using elastic net
Rscript machine_learning.R 
sed -i 's/S1/S2/g' ${id} machine_learning.R
Rscript machine_learning.R
sed -i 's/S2/S3/g' ${id} machine_learning.R
Rscript machine_learning.R
sed -i 's/S3/S4/g' ${id} machine_learning.R
Rscript machine_learning.R
sed -i 's/S4/S5/g' ${id} machine_learning.R
Rscript machine_learning.R