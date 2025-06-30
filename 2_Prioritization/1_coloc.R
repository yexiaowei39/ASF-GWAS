library(tidyverse)
library(coloc)
library(data.table)
tissues<-c("Adipose","Artery","Blastocyst","Blastomere","Blood","Brain","Cartilage",
           "Colon","Duodenum","Embryo","Fetal_thymus","Frontal_cortex","Heart",
           "Hypothalamus","Ileum","Jejunum","Kidney","Large_intestine","Liver",
           "Lung","Lymph_node","Macrophage","Milk","Morula","Muscle","Oocyte",
           "Ovary","Pituitary","Placenta","Small_intestine","Spleen",
           "Synovial_membrane","Testis","Uterus")

for (group in 1:4) {
  
  gene_gwas<-read.table(paste0(group,'_gene_gwas_filter.txt'),header=T)   # all loci
  genes<-unique(gene_gwas$gene_id)
  gwas<-read.table(paste0(group,'_significant_loci_1MB.txt'),header=T)    # gwas upstream and downstream 1MB loci
  
  s <- ifelse(group == 1, 0.77,
              ifelse(group %in% c(2, 3), 0.3,
                     ifelse(group == 4, 0.47, NA)))  
  
  
  for(tissue in tissues){
    
    eqtl_dir<-paste0("/eQTL/",tissue)
    eqtls<-list.files(eqtl_dir,full.names = T)   
    eqtls2<-str_glue("{eqtl_dir}/{genes}.txt")
    eqtls<-eqtls[eqtls%in%eqtls2]
    
    
    if(length(eqtls)>0){
      res<-map(eqtls, ~{
        eqtl<-fread(.x,h=T)
        pos_eqtl<-data.frame(              
          SNP=str_match(eqtl$variant_id,"(\\d+_\\d+)_")[,2],
          ref=str_match(eqtl$variant_id,"\\d+_\\d+_(\\w+)_\\w+")[,2],
          alt=str_match(eqtl$variant_id,"\\d+_\\d+_\\w+_(\\w+)")[,2])
        
        gwas2<-gwas[gwas$SNP%in%pos_eqtl$SNP,]              
        
        
        if(nrow(gwas2)>0){
          
          gwas2<-merge(gwas2,pos_eqtl,by="SNP")              ## merge gwas & eqtl loci
          gwas2<-gwas2[complete.cases(gwas2),]
          not_cons<-gwas2$allele1!=gwas2$ref
          
          gwas_coloc<-gwas2
          gwas_coloc[not_cons,]$beta<--gwas_coloc[not_cons,]$beta
          snpid<-paste0(gwas_coloc$SNP,"_",gwas_coloc$ref,"_",gwas_coloc$alt)
          eqtl_coloc<-eqtl[match(snpid,eqtl$variant_id),]
          
          
          # coloc
          gene.coloc<-list(
            pvalues=eqtl_coloc$pval_nominal,
            beta=eqtl_coloc$slope,
            varbeta=eqtl_coloc$slope_se^2,
            snp=eqtl_coloc$variant_id,
            type="quant",
            N=eqtl_coloc$ma_samples,
            MAF=eqtl_coloc$af)
          
          gwas.coloc<-list(
            pvalues=gwas_coloc$p_wald,
            beta=gwas_coloc$beta,
            varbeta=gwas_coloc$se^2,
            snp=snpid,
            type="cc",
            N=474,
            s=s,
            MAF=gwas_coloc$af)
          
          
          coloc_x = coloc.abf(gene.coloc, gwas.coloc)
          
          res<-data.frame(Tissue=tissue,Gene=eqtl_coloc$phenotype_id[1],Group=group,
                          nsnps=coloc_x$summary[1],
                          PP0=coloc_x$summary[2],
                          PP1=coloc_x$summary[3],
                          PP2=coloc_x$summary[4],
                          PP3=coloc_x$summary[5],
                          PP4=coloc_x$summary[6])
          
          return(res)
          
        }else{return(NULL)}
        
      })
      
      res <- bind_rows(res)
      
      if(nrow(res)>0){
        write.table(res,paste0("./1_output/",group,"_",tissue,".txt"),row.names=F,quote=F)
      }
      
    }
    
    cat("Tissue: ", tissue, " ", nrow(res), "\n")
  }
}

