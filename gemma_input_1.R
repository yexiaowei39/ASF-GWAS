library(data.table)
options <- commandArgs(trailingOnly = TRUE)

grm<-fread(options[1],h=F)
size <- 474
grm_mat<-matrix(0,size,size)

grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
grm_mat<-t(grm_mat)
grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4

write.table(grm_mat,file=options[2],row.names=F,col.names = F,quote=F)
