#!/bin/bash
#
#SBATCH --job-name=GCTA_h2
#SBATCH --output=GCTA_h2_mpi.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

##Real estimate for GWAS heritability
for i in $(seq 1 144);
do
gcta64 --grm /scratch/users/s/h/shifang/GCTA/test --pheno /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/T1_trait_new.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/GCTA_qcovar.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/GCTA_covar.txt
done

##Resampling 10 times
for j in $(seq 1 10);
do 
for i in $(seq 1 144);
do
plink --bfile /scratch/users/s/h/shifang/GCTA/genotypes --keep /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID --chr-set 29 --allow-extra-chr

plink --bfile ${j%}_GCTA_ID --indiv-sort f /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID_1 --chr-set 29 --allow-extra-chr

gcta64 --bfile ${j%}_GCTA_ID_1 --autosome-num 29 --make-grm --out test

gcta64 --grm test --pheno /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_pheno_.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_qcov.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_cov.txt

find . -type f -name "*.hsq" >111.txt
find . -type f -name "*.hsq" -exec awk 'NR==5' {} \; > output.txt

paste 111.txt output.txt | awk '{print $1,$2,$3,$4}' >/scratch/users/s/h/shifang/GCTA/h2/resampling/test/output/mergedfile_${j%}.txt
rm *.bed *.bim *.fam *.log *.bin *.id *.hsq *.txt
done;
done

##permutation test with 1000 times 
for j in $(seq 1 144);
do for i in $(seq 1 1000); 
do
gcta64 --grm /scratch/users/s/h/shifang/GCTA/test --pheno /scratch/users/s/h/shifang/GCTA/h2/perm/raw/${j%}.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/GCTA_qcovar.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/GCTA_covar.txt 
done
rm *.log
find . -type f -name "*.hsq" -exec awk 'NR==5' {} \; > /scratch/users/s/h/shifang/GCTA/h2/perm/perm/perm_${j%}_output.txt
rm *.hsq
done
