#!/bin/bash

for i in $(seq 1 144); do
    echo MMMMMMMMMMMMMMMMMMMMMMMMMM  $i  MMMMMMMMMMMMMMMMMMMMMMMMMMM
    date ;
    gcta64 --mlma-loco --bfile genotypes --pheno /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/T1_trait_new.txt --out WG_LOCO_T1_${i} --qcovar /scratch/users/s/h/shifang/GWAS/GCTA_cov.txt --mpheno ${i} --thread-num 11 --autosome-num 29 ;
    date ;
done