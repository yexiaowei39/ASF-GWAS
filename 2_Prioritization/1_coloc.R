library(tidyverse)
library(coloc)
library(data.table)

          gene.coloc<-list(
            pvalues=eqtl_coloc$p,
            beta=eqtl_coloc$beta,
            varbeta=eqtl_coloc$se^2,
            snp=eqtl_coloc$variant_id,
            type="quant",
            N=eqtl_coloc$ma_samples,
            MAF=eqtl_coloc$af)
          
          gwas.coloc<-list(
            pvalues=gwas_coloc$p,
            beta=gwas_coloc$beta,
            varbeta=gwas_coloc$se^2,
            snp=snpid,
            type="cc",
            N=n,
            s=s,
            MAF=gwas_coloc$af)
          
          
          coloc_x = coloc.abf(gene.coloc, gwas.coloc)
