group=1

for chr in $(seq 1 18)
do

  # 1 FST
  vcffile_all=chr${chr}_merge.vcf.gz   # vcf data
  $vcftools --gzvcf $vcffile_all --weir-fst-pop $case.txt --weir-fst-pop $control.txt --out ${group}_fst_chr$chr --fst-window-size 100000 --fst-window-step 10000

  # 2 assoc / allele frequence chi-square test
  plink_data_out=chr${chr}_case_pheno  # plink data  # with phenotype
  $plink --bfile $plink_data_out --double-id --allow-no-sex --assoc --adjust --out $plink_out
  
done
