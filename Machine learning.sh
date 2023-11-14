#!/bin/bash

for j in $(seq 1 10);
do for i in 55 59 72 76 93 94 95 106 110 123 127 128 129 131 140 144 145 146 157 166 167 172 178 179 180 195 208 212 214 219 221 223 229 57 74 108 125 142 159 176 193 210 227 248 $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206);
do
plink --bfile /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/geno_mind_test --keep /scratch/users/s/h/shifang/machine_learning/resampling/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID --chr-set 29 --allow-extra-chr
plink --bfile ${j%}_GCTA_ID --indiv-sort f /scratch/users/s/h/shifang/machine_learning/resampling/${j%}_GCTA_ID.txt --make-bed --out ${j%}_GCTA_ID_1 --chr-set 29 --allow-extra-chr
gcta64 --mlma-loco --bfile ${j%}_GCTA_ID_1 --pheno /scratch/users/s/h/shifang/machine_learning/resampling/${j%}_pheno_.txt --out mc_${i%}_loco_${j%} --thread-num 12 --qcovar /scratch/users/s/h/shifang/machine_learning/resampling/${j%}_GCTA_qcov.txt --mpheno $i --covar /scratch/users/s/h/shifang/machine_learning/resampling/${j%}_GCTA_cov.txt
awk '$9<=5E-05' mc_${i%}_loco_${j%}.loco.mlma >/scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig/mc_${i%}_loco_${j%}_sig.subset.txt
#rm -r *bed *bim * fam * log *nosex *mlma
done
done

cd /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig
for j in $(seq 1 10);
do for i in 55 59 72 76 93 94 95 106 110 123 127 128 129 131 140 144 145 146 157 166 167 172 178 179 180 195 208 212 214 219 221 223 229 57 74 108 125 142 159 176 193 210 227 248 $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206);
do
awk '{print $2}' mc_${i%}_loco_${j%}_sig.subset.txt > SNPs.txt
vcftools --vcf /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/geno_mind_test_vcf.vcf --snps SNPs.txt --recode --recode-INFO-all --out test1
plink --vcf test1.recode.vcf --keep-allele-order --double-id --make-bed --out test1_1 --chr-set 29 --allow-extra-chr --set-missing-var-ids @#
plink --bfile test1_1 --recodeA --out test1_1_geno --chr-set 29 --allow-extra-chr
head test1_1_geno.raw| cut -d" " -f 1-9
cat test1_1_geno.raw | cut -d" " -f2,7- | sed 's/_[A-Z]//g' >/scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig/geno/mc_${i%}_loco_${j%}_genotype.txt
mv mc_${i%}_loco_${j%}_sig.subset.txt /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig/sig
#rm -r *bed *bim * fam * log *nosex *mlma *vcf
done
done

for i in 55 59 72 76 93 94 95 106 110 123 127 128 129 131 140 144 145 146 157 166 167 172 178 179 180 195 208 212 214 219 221 223 229 57 74 108 125 142 159 176 193 210 227 248 $(seq 43 53) $(seq 60 70) $(seq 97 104) $(seq 133 138) $(seq 148 155) $(seq 111 121) $(seq 182 189) $(seq 196 206);
do
mv /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig/mc_${i%}_loco_1_genotype.txt /scratch/users/s/h/shifang/imputation/chr23/DATA/Merged/All_merged/GZ_TBIX/mc/sig/T1
done

cat resampling_times | while read id
do
Rscript --vanilla machine_learning.R ${id}
done