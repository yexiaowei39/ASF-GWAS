group=1

for chr in $(seq 1 18)
do

  plink_case_out=chr${chr}_case_pheno           ## individual with phenotype

  ## GWAS

    # 1 make grm
    $gcta --bfile $plink_case_out --make-grm --make-grm-gz --thread-num 10 --out chr${chr}_case_grm_1

    # 2 merge grm   ## leave one chr out
    $gcta --mgrm-gz grm_case_chr${chr}.txt --make-grm-gz --out case_kin_chr$chr
    Rscript gemma_input.R case_kin_chr$chr.grm.gz case_kin_chr$chr.txt

    # 3 regression
    $gemma -bfile $plink_case_out -k case_kin_chr$chr.txt -lmm 1 -o case_chr$chr -outdir out -c $gwas_cov
    
done
