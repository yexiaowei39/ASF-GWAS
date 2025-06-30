group=1
mkdir ./fst/case/$group
mkdir ./fst/control/$group


for chr in $(seq 1 18)
do

  # 1 FST
  vcffile_all=chr${chr}_merge.vcf.gz   # vcf data
  
  $vcftools --gzvcf $vcffile_all --weir-fst-pop $case.txt --weir-fst-pop $control1.txt --out ./fst/case/$group/${group}_fst_chr$chr --fst-window-size 100000 --fst-window-step 10000
  $vcftools --gzvcf $vcffile_all --weir-fst-pop $case.txt --weir-fst-pop $control2.txt --out ./fst/control/$group/${group}_fst_chr$chr --fst-window-size 100000 --fst-window-step 10000


  # 2 assoc / allele frequence chi-square test
  plink_control_1=./assoc/control/1/1_out/chr${chr}_out
  
  $plink --bfile $plink_data_out --double-id --allow-no-sex --assoc --adjust --out $plink_control_1
  
done
