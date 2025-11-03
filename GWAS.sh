#!/bin/bash
#
#SBATCH --job-name=pQTL mapping
#SBATCH --output=pQTL_mpi.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000

##Using a leave-one-chrom-out mixed linear model with GRM as the random variable and age and season as fixed effects.
for i in $(seq 1 144); do
    echo MMMMMMMMMMMMMMMMMMMMMMMMMM  $i  MMMMMMMMMMMMMMMMMMMMMMMMMMM
    date ;
    gcta64 --mlma-loco --bfile /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/genotypes --pheno /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/T1_trait_new.txt --out WG_LOCO_T1_${i} --qcovar /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/GCTA_cov.txt --mpheno ${i} --thread-num 10 --autosome-num 29 ;
    date ;
done

##extract significant and suggestive significant loci
cat config | while read id
do
awk '$9<=1.30890052356021E-07' ${id}.loco.mlma >/sig/${id}.subset.txt
awk '$9<=2.61780104712042E-06' ${id}.loco.mlma >suggestive/${id}.subset.txt
done
find sig -type f -empty -delete
find suggestive -type f -empty -delete

##condition analysis
cd sig
mkdir cojo
ls *txt|awk '{print substr($0,1,length($0)-11)}' | sort | uniq > config
cat config | while read id
do
awk '{print $1}' ${id}.subset.txt>chr.txt
sort chr.txt | uniq > output_file
cat output_file | while read id1
do
awk '{print $2,$4,$5,$6,$7,$8,$9}' /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/${id}.loco.mlma >condition.mlma
awk -F '\t' -v OFS='\t' '{ $(NF+1) = 246; print }' condition.mlma >outfile
gcta64  --bfile /scratch/users/s/h/shifang/GWAS/immunophenotypes/T1/genotypes  --chr ${id1} --autosome-num 29 --thread-num 10 --maf 0.05 --cojo-file outfile --cojo-slct --out cojo/${id}_${id1} --cojo-p 1.30890052356021E-07
done
rm -r chr.txt output_file condition.mlma
done

