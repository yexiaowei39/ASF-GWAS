# library(GALLO)
# library(dplyr)
# qtf.inp <- import_gff_gtf(db_file="/disk212/yexw/index/Sus_scrofa.Sscrofa11.1.105.gtf",file_type="gtf")
# 
# for (group in 2:4) {
#   # 1.所有位点GWAS摘要统计量 
#   gwas_all<-read.table(paste0("/disk212/yexw/wgs/ASF/3_2/2_2_TWAS/0_input/",group,'_input.txt.gz'),header=T)
#   
#   # 2.过滤后的GWAS位点
#   gwas.data<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_gwas_filter.txt'),header=T)
#   
#   marker_gwas<-data.frame(
#     SNP=paste0(gwas.data$chr,"_",gwas.data$ps),
#     CHR=gwas.data$chr,
#     BP=gwas.data$ps,
#     P=gwas.data$p_wald,
#     MAF=gwas.data$af
#   )
#   gene_gwas<-find_genes_qtls_around_markers(db_file=qtf.inp,marker_file=marker_gwas,
#                                             method="gene",marker="snp",interval=1000000,nThreads=20)
#   write.table(gene_gwas,paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_gene_gwas_filter.txt'),quote=F,row.names=F)
#   
#   # 3.匹配显著位点上下游1MB的所有位点
#   gwas_significant <- marker_gwas %>%
#     mutate(region_start = BP - 1e6,  # 上下游1MB
#            region_end = BP + 1e6)
#   
#   gwas <- gwas_all %>%
#     rowwise() %>% # 行遍历
#     filter(any(gwas_significant$CHR == chr & 
#                  ps >= gwas_significant$region_start & 
#                  ps <= gwas_significant$region_end))
#   gwas$SNP<-paste0(gwas$chr,"_",gwas$ps)
#   
#   write.table(gwas,paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_significant_loci_1MB.txt'),quote=F,row.names=F)
# }


# eQTL #######################
library(tidyverse)
library(coloc)
library(data.table)
tissues<-c("Adipose","Artery","Blastocyst","Blastomere","Blood","Brain","Cartilage",
           "Colon","Duodenum","Embryo","Fetal_thymus","Frontal_cortex","Heart",
           "Hypothalamus","Ileum","Jejunum","Kidney","Large_intestine","Liver",
           "Lung","Lymph_node","Macrophage","Milk","Morula","Muscle","Oocyte",
           "Ovary","Pituitary","Placenta","Small_intestine","Spleen",
           "Synovial_membrane","Testis","Uterus")
# group<-1
# tissue<-c("Macrophage")

for (group in 1:4) {
  
  gene_gwas<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_gene_gwas_filter.txt'),header=T)
  genes<-unique(gene_gwas$gene_id)
  gwas<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_significant_loci_1MB.txt'),header=T)
  
  s <- ifelse(group == 1, 0.77,
              ifelse(group %in% c(2, 3), 0.3,
                     ifelse(group == 4, 0.47, NA)))  
  
  
  for(tissue in tissues){
    
    eqtl_dir<-paste0("/disk224/zz/TB_metagenome/mGWAS_Metaphlan_v4/gemma/coloc2/eQTL/",tissue)
    # eqtl_dir<-c("/disk212/yexw/wgs/ASF/3_2/2_3_colol/1")
    eqtls<-list.files(eqtl_dir,full.names = T)   
    eqtls2<-str_glue("{eqtl_dir}/{genes}.txt")
    eqtls<-eqtls[eqtls%in%eqtls2]
    
    
    if(length(eqtls)>0){
      res<-map(eqtls, ~{
        eqtl<-fread(.x,h=T)
        # eqtl<-read.table(eqtls,h=T)        ## 读入eqtl文件
        pos_eqtl<-data.frame(              ## eqtl 位点信息
          SNP=str_match(eqtl$variant_id,"(\\d+_\\d+)_")[,2],
          ref=str_match(eqtl$variant_id,"\\d+_\\d+_(\\w+)_\\w+")[,2],
          alt=str_match(eqtl$variant_id,"\\d+_\\d+_\\w+_(\\w+)")[,2])
        
        gwas2<-gwas[gwas$SNP%in%pos_eqtl$SNP,]               ## 合并gwas和eqtl位点
        
        
        if(nrow(gwas2)>0){
          
          gwas2<-merge(gwas2,pos_eqtl,by="SNP")              ## 合并gwas和eqtl位点
          gwas2<-gwas2[complete.cases(gwas2),]
          not_cons<-gwas2$allele1!=gwas2$ref
          
          # 对于颠倒的allele，取beta的相反数
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
        write.table(res,paste0("/disk212/yexw/wgs/ASF/3_2/2_3_colol/1_output/",group,"_",tissue,".txt"),row.names=F,quote=F)
      }
      
    }
    
    cat("Tissue: ", tissue, " ", nrow(res), "\n")
  }
}



# 可视化 #########################################################
library(GALLO)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(data.table)
library(ggbio)
library(GenomicRanges)
# qtf.inp <- import_gff_gtf(db_file="/disk212/yexw/index/Sus_scrofa.Sscrofa11.1.105.gtf",file_type="gtf")  
plink<-"/disk191/chenzt/software/plink"


# test
# group<-1
# tissue <- "Cartilage"
# gene <- "ENSSSCG00000023771"

# 定义画图函数
plotColoc<-function(group,tissue,gene){

    eqtl_dir<-paste0("/disk211/weiran/others/data/eQTL/",tissue)
    eqtls<-list.files(eqtl_dir,full.names = T)   
    eqtls2<-str_glue("{eqtl_dir}/{gene}.txt")
    eqtls<-eqtls[eqtls%in%eqtls2]
      
      
        outdir<-paste0("/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/",group,"_",tissue,"_",gene_name)
        dir.create(outdir,showWarnings = F)
    
        eqtl<-read.table(paste0("/disk211/weiran/others/data/eQTL/",tissue,"/",gene,".txt"),header=T)
        pos_eqtl<-data.frame(
          SNP=str_match(eqtl$variant_id,"(\\d+_\\d+)_")[,2],
          ref=str_match(eqtl$variant_id,"\\d+_\\d+_(\\w+)_\\w+")[,2],
          alt=str_match(eqtl$variant_id,"\\d+_\\d+_\\w+_(\\w+)")[,2],
          pval_nominal=eqtl$pval_nominal
        )
        gwas2<-gwas[gwas$SNP%in%pos_eqtl$SNP,]
        gwas2<-merge(gwas2,pos_eqtl,by="SNP")    #合并gwas和eqtl位点
      
        gwas_coloc<-gwas2
    
        snpid<-paste0(gwas_coloc$SNP,"_",gwas_coloc$ref,"_",gwas_coloc$alt)
        eqtl_coloc<-eqtl[match(snpid,eqtl$variant_id),]
        
        
        # 作图
        ## LD: GWAS significant loci 左右1Mb的
        # tag<-gene_gwas[gene_gwas$gene_id==gene,]$SNP
        tag<-gwas2[which.min(gwas2$p_wald)[1],]$SNP
        chr<-as.numeric(str_match(tag,"(\\d+)_\\d+")[,2])
        pos<-as.numeric(str_match(tag,"\\d+_(\\d+)")[,2])
        tag<-paste0("chr",chr,":",pos)
        
        system(paste0(plink,' --bfile /disk212/yexw/wgs/ASF/3_2/0_data/case/',group,'/chr',chr,'_case_pheno --ld-snp ',tag,' --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out ',outdir,'/gwas'))
        ld<-read.table(str_glue("{outdir}/gwas.ld"),h=T)
        
        gws<-filter(gwas_coloc,rs%in%ld$SNP_B)
        gws$logp<--log10(gws$p_wald)
        gws$r2<-ld$R2[match(gws$rs,ld$SNP_B)]
        gws<-dplyr::select(gws,chr,ps,logp,r2)
        gws$target_snp<-tag
        gws$col[gws$r2<=0.2]<-0.2
        gws$col[gws$r2>0.2 & gws$r2<=0.4]<-0.4
        gws$col[gws$r2>0.4 & gws$r2<=0.6]<-0.6
        gws$col[gws$r2>0.6 & gws$r2<=0.8]<-0.8
        gws$col[gws$r2>0.8 & gws$r2<=1]<-1
      
        pos_eqtl$chr<-as.numeric(str_match(pos_eqtl$SNP,"(\\d+)_\\d+")[,2])
        pos_eqtl$BP<-as.numeric(str_match(pos_eqtl$SNP,"\\d+_(\\d+)")[,2])
        pos_eqtl$logp<--log10(pos_eqtl$pval_nominal)
        pos_eqtl$SNP_B<-paste0("chr",pos_eqtl$chr,":",pos_eqtl$BP)
        pos_eqtl$r2<-ld$R2[match(pos_eqtl$SNP_B,ld$SNP_B)]
        
        eqtl2<-dplyr::select(pos_eqtl,chr,BP,logp,r2)
        eqtl2<-eqtl2[complete.cases(eqtl2),]
        eqtl2$target_snp<-tag
        eqtl2$col[eqtl2$r2<=0.2]<-0.2
        eqtl2$col[eqtl2$r2>0.2 & eqtl2$r2<=0.4]<-0.4
        eqtl2$col[eqtl2$r2>0.4 & eqtl2$r2<=0.6]<-0.6
        eqtl2$col[eqtl2$r2>0.6 & eqtl2$r2<=0.8]<-0.8
        eqtl2$col[eqtl2$r2>0.8 & eqtl2$r2<=1]<-1
    
        overlap<-merge(gws,eqtl2,by.x="ps",by.y="BP")
        
        subset=overlap[overlap$ps==pos,]
        # res=cor.test(overlap$logp.y, overlap$logp.x)
        # p=round(res$p.value,3)
        # r=round(res$estimate,2) # 0.45

      
        gwas_plot<-
            ggplot(overlap,aes(ps/1000000,logp.x))+
            geom_point(aes(color=factor(col.x)),alpha=0.8)+
            geom_point(data=subset,aes(ps/1000000,logp.x),shape=23,size=4,fill="purple")+
            xlab("")+
            ylab(expression(GWAS-log(italic(P))))+
            # ggtitle(group)+
            theme(panel.background=element_rect(fill='white', colour='black'),
                  legend.position = "none",
                  axis.text.x=element_blank(),
                  plot.margin=unit(c(0.4,0.2,0.2,0.4),"cm"))+
            scale_color_manual(values = c("0.2" = "blue4", "0.4" = "skyblue", "0.6" = "darkgreen", 
                                          "0.8" = "orange", "1" = "red"))+
            geom_text_repel(data=subset,aes(label=target_snp.x),
                            size=3,nudge_x=0.3,
                            box.padding=unit(0.35, "lines"),
                            point.padding=unit(0.3, "lines"))
        
          eqtl_plot<-
            ggplot(overlap,aes(ps/1000000,logp.y))+
            geom_point(aes(color=factor(col.y)),alpha=0.8)+
            scale_color_manual(values = c("0.2" = "blue4", "0.4" = "skyblue", "0.6" = "darkgreen", 
                                          "0.8" = "orange", "1" = "red"))+
            geom_point(data=subset,aes(ps/1000000,logp.y),shape=23,size=4,fill="purple")+
            xlab("")+
            ylab(expression(eQTL-log(italic(P))))+
            ggtitle(paste0(tissue,": ",gene_name))+
            theme(panel.background=element_rect(fill='white', colour='black'),
                  legend.position = "none",
                  axis.text.x=element_blank(),
                  plot.margin=unit(c(0.2,0.2,0.2,0.4),"cm"))+
            geom_text_repel(data=subset,aes(label=target_snp.x),
                            nudge_x=0.3,size=3,
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.3, "lines"))
        
          
          overlap_plot<-
            ggplot(overlap,aes(logp.y,logp.x))+
            geom_point(aes(color=factor(col.x)),alpha=0.8)+ #x eqtl; y gwas
            scale_color_manual(values = c("0.2" = "blue4", "0.4" = "skyblue", "0.6" = "darkgreen", 
                                          "0.8" = "orange", "1" = "red"))+
            geom_point(data=subset,aes(logp.y,logp.x),shape=23,size=4,fill="purple")+
            xlab(expression(eQTL-log(italic(P))))+
            ylab(expression(GWAS-log(italic(P))))+
            labs(color = bquote(~r^2))+
            theme(panel.background=element_rect(fill='white', colour='black'),
                  plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm"))+
            # guides(color=guide_legend(override.aes=list(shape=15,col=cols,size=4),title=bquote(~r^2)))+
            geom_smooth(method=lm,se=FALSE)+
            geom_text_repel(data=subset,aes(label=target_snp.x),
                            size=3,nudge_y=0.5,
                            box.padding=unit(0.35, "lines"),
                            point.padding=unit(0.3, "lines"))
          
          
          
          
          # 准备基因位置数据
          allgenes<-gene_gwas[gene_gwas$CHR==chr,]
          allgenes<-allgenes[!duplicated(allgenes$gene_id),]
          # 筛选符合条件的基因
          selected_genes <- allgenes %>%
            filter(start_pos <= max(overlap$ps) &  # 确保基因起始位置小于等于 overlap 最大位点
                     end_pos >= min(overlap$ps))     # 确保基因结束位置大于等于 overlap 最小位点
          selected_genes$gene_name<-ifelse(!is.na(selected_genes$gene_name),selected_genes$gene_name,selected_genes$gene_id)
          
          genes <- data.frame(seqnames = paste0("chr",chr),
                              start = selected_genes$start_pos,
                              end = selected_genes$end_pos,
                              gene_name = selected_genes$gene_name)
          
          # 将数据转换为 GRanges 对象
          genes_gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
          
          # 为基因分组：按每行 5 个基因排列
          desired_rows <- 5 
          genes$group <- cut(seq_along(genes$gene_name), 
                             breaks = desired_rows, 
                             labels = FALSE)
          genes$color<-ifelse(genes$gene_name==gene_name,"red","black")
          
          gene_pos<-ggplot()+
            geom_rect(data=genes,aes(xmin=start/1000000,xmax=end/1000000,ymin=group-0.3,ymax=group+0.3,color=color), 
                      fill="lightblue") +
            geom_text_repel(data=genes,aes(x=(start/1000000+end/1000000)/2,y=group,label=gene_name,color=color),
                            size=3,box.padding=unit(0.35,"lines"),point.padding=unit(0.3,"lines"))+
            scale_color_manual(values=c("red"="red","black"="black"))+
            xlab(paste("chr",chr,"position (Mb)"))+
            ylab("gene position")+
            theme(panel.background=element_rect(fill='white', colour='black'),
                  axis.text.y = element_blank(), 
                  # axis.title.y = element_blank(),
                  legend.position = "none",
                  plot.margin=unit(c(0.2,0.2,0.4,0.4),"cm"))
          
          
          
      
    
          # 拼图
          p<-ggarrange(ggarrange(gwas_plot,eqtl_plot,gene_pos,nrow=3,labels=c("(a)","(b)","(c)"),
                                 font.label=list(color="black",face="bold",size=14),
                                 align = "v"),
                       ggarrange(overlap_plot,labels="(d)",
                                 font.label=list(color="black",face="bold",size=14)),
                       ncol=2,widths=c(1,1))+
              theme(plot.margin=margin(10,10,10,10),units = "in")
          p<-annotate_figure(p,top=text_grob(paste0("Comparison ",group),color="black",face="bold",size=16,
                                          hjust=0.5,                   # 水平对齐方式：居中
                                          vjust=1))                     # 垂直对齐方式：靠近顶部+ 
        
          file_name_plot<-paste0(outdir,"/coloc.pdf")
          ggsave(file_name_plot,plot = p,device="pdf",width=10,height=8,units = "in") 
          
          
          write.table(overlap,file=paste0(outdir,"/gwas_eqtl_ovlp.txt"),row.names=F,quote=F)
}


coloc_sig_all<-NULL
# coloc_sig_all<-read.table("/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc.txt",header = T)
# outputs<-list.files("/disk212/yexw/wgs/ASF/3_2/2_3_colol/1_output/",full.names = T)   

tissues<-c("Adipose","Artery","Blastocyst","Blastomere","Blood","Brain","Cartilage",
           "Colon","Duodenum","Embryo","Fetal_thymus","Frontal_cortex","Heart",
           "Hypothalamus","Ileum","Jejunum","Kidney","Large_intestine","Liver",
           "Lung","Lymph_node","Macrophage","Milk","Morula","Muscle","Oocyte",
           "Ovary","Pituitary","Placenta","Small_intestine","Spleen",
           "Synovial_membrane","Testis","Uterus")
# group<-3
# tissue<-"Brain"
# gene<-"ENSSSCG00000007687"

# for (output in outputs) {

for (group in c(1:4)) {
  
# 组内所有显著位点附近 1MB位点
        # gwas<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_significant_loci_1MB.txt'),header=T)
        # gene_gwas<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_gene_gwas_filter.txt'),header=T)
        # gwas.data<-read.table(paste0('/disk212/yexw/wgs/ASF/3_2/2_3_colol/0_input/',group,'_gwas_filter.txt'),header=T)
        
        
for (tissue in tissues) {
  
    # coloc<-read.table(output,header=T)
    coloc<-read.table(paste0("/disk212/yexw/wgs/ASF/3_2/2_3_colol/1_output/",group,"_",tissue,".txt"),header=T)
    coloc_sig<-coloc[coloc$PP4>=0.75 | coloc$PP3>=0.75,]
    
  if(nrow(coloc_sig)>0){
    # for(i in 1:nrow(coloc_sig)){
    #   
    #     # group<-coloc_sig$Group[i]
    #     tissue<-coloc_sig$Tissue[i]
    #     gene<-coloc_sig$Gene[i]
    #     
    #     gene_name<-gene_gwas$gene_name[match(gene,gene_gwas$gene_id)][1]
    #     gene_name<-ifelse(!is.na(gene_name),gene_name,gene)
    #     
    #     if(!file.exists(paste0("/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/",group,"_",tissue,"_",gene_name))){
    #       plotColoc(group,tissue,gene)
    #     }
    #     
    #     sig_loci<-gene_gwas[gene_gwas$gene_id==gene,]
    #     sig_loci<-sig_loci[which.min(sig_loci$P),]$SNP
    #     
    #     
    #     coloc_sig$CHR<-as.numeric(str_match(sig_loci,"(\\d+)_\\d+")[,2])
    #     coloc_sig$BP<-as.numeric(str_match(sig_loci,"\\d+_(\\d+)")[,2])
    # }
    
    coloc_sig_all<-rbind(coloc_sig_all,coloc_sig)
  }else{print(paste0(group,"_",tissue," no sig"))}
}
  
}

write.table(coloc_sig_all,"/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc_PP4.txt",row.names=F,quote=F)
write.table(coloc_sig_all,"/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc_PP3&4.txt",row.names=F,quote=F)


###### 统计数量 ##########################
library(dplyr)
coloc_sig_all<-read.table("/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc.txt",header = T)
gene_id_all<-read.table("/disk222/yexw/2_ASF/2_RNA/gene_id.txt",header=T)
gene_id_all<-gene_id_all[,c("gene_id","gene_name")]
gene_id_all<-gene_id_all[!duplicated(gene_id_all$gene_id),]
gene_id_all$gene_name<-ifelse(is.na(gene_id_all$gene_name),gene_id_all$gene_id,gene_id_all$gene_name)

coloc_sig_gene <- coloc_sig_all %>%
  left_join(gene_id_all, by = c("Gene" = "gene_id"))
write.table(coloc_sig_gene,"/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc_PP3&4.txt",row.names=F,quote=F)

coloc_sig_gene<-read.table("/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/all_significant_coloc_PP3&4.txt",header=T)
coloc_sig_gene<-coloc_sig_gene[complete.cases(coloc_sig_gene$gene_name),]
coloc_counts <- coloc_sig_gene %>%
  group_by(gene_name) %>%
  summarise(gene_count = n())
write.table(coloc_counts,"/disk212/yexw/wgs/ASF/3_2/2_3_colol/2_plot/coloc_PP3&4_count.txt",row.names=F,quote=F,sep="\t")



