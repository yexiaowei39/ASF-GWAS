group=1

for chr in $(seq 1 18)
do

  # 1 FST
  $vcftools --gzvcf $vcf --weir-fst-pop $case --weir-fst-pop $control --out ${group}_fst_chr$chr --fst-window-size 100000 --fst-window-step 10000

  # 2 assoc / allele frequence chi-square test
  plink_data_out=chr${chr}_case_pheno 
  $plink --bfile $plink_data_out --assoc --adjust --out $plink_out
  
done
