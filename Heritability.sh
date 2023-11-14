#!/bin/bash

##Real estimate for GWAS heritability
for i in $(seq 14 25) $(seq 27 29) $(seq 33 37) $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206) $(seq 11 13);
do
gcta64 --grm /scratch/users/s/h/shifang/GCTA/test --pheno /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/T1_trait_new.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/GCTA_qcovar.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/GCTA_covar.txt
done

##Resampling 10 times
for j in $(seq 1 10);
do for i in 55 59 72 76 93 94 95 106 110 123 127 128 129 131 140 144 145 146 157 166 167 172 178 179 180 195 208 212 214 219 221 223 229 1 2 7 8 9 57 74 108 125 142 159 176 193 210 227 248 $(seq 14 25) $(seq 27 29) $(seq 33 37) $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206);
do
plink --bfile /scratch/users/s/h/shifang/GCTA/genotypes --keep /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID --chr-set 29 --allow-extra-chr

plink --bfile ${j%}_GCTA_ID --indiv-sort f /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID_1 --chr-set 29 --allow-extra-chr

gcta64 --bfile ${j%}_GCTA_ID_1 --autosome-num 29 --make-grm --out test

gcta64 --grm test --pheno /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_pheno_.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_qcov.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/h2/resampling/raw/${j%}_GCTA_cov.txt
done
find . -type f -name "*.hsq" >111.txt
find . -type f -name "*.hsq" -exec awk 'NR==5' {} \; > output.txt

paste 111.txt output.txt | awk '{print $1,$2,$3,$4}' >/scratch/users/s/h/shifang/GCTA/h2/resampling/test/output/mergedfile_${j%}.txt
rm *.bed *.bim *.fam *.log *.bin *.id *.hsq *.txt
done;
done

##permutation test with 1000 times 
for j in 11 12 13 55 59 72 76 93 94 95 106 110 123 127 128 129 131 140 144 145 146 157 166 167 172 178 179 180 195 208 212 214 219 221 223 229 1 2 7 8 9 57 74 108 125 142 159 176 193 210 227 248 $(seq 14 25) $(seq 27 29) $(seq 33 37) $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206) ;
do for i in $(seq 1 1000); 
do
gcta64 --grm /scratch/users/s/h/shifang/GCTA/test --pheno /scratch/users/s/h/shifang/GCTA/h2/perm/raw/${j%}.txt --reml --qcovar /scratch/users/s/h/shifang/GCTA/GCTA_qcovar.txt --out h_${i%} --thread-num 10 --mpheno $i --covar /scratch/users/s/h/shifang/GCTA/GCTA_covar.txt 
done
rm *.log
find . -type f -name "*.hsq" -exec awk 'NR==5' {} \; > /scratch/users/s/h/shifang/GCTA/h2/perm/perm/perm_${j%}_output.txt
rm *.hsq
done
