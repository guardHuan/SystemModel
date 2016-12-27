tumor<-read.csv('../d_tumor.csv')
normal<-read.csv('../d_normal.csv')
rownames(tumor)<-tumor[,1]
rownames(normal)<-normal[,1]
tumor<-tumor[,-1]
normal<-normal[,-1]

library('synRNASeqNet')
counts<-as.matrix(t(tumor))
tumor_gene_cor<-parMIEstimate(counts,method='KNN',nchips=48)
counts<-as.matrix(t(normal))
normal_gene_cor<-parMIEstimate(counts,method='KNN',nchips=48)

normal_gene_cor[lower.tri(normal_gene_cor)]<-NA
diag(normal_gene_cor)<-NA
tumor_gene_cor[lower.tri(tumor_gene_cor)]<-NA
diag(tumor_gene_cor)<-NA

library(reshape)

list_tumor_cor<-melt.array(tumor_gene_cor)
list_normal_cor<-melt.array(normal_gene_cor)
list_tumor_cor2<-list_tumor_cor[complete.cases(list_tumor_cor$value)&complete.cases(list_normal_cor$value),]
list_normal_cor2<-list_normal_cor[complete.cases(list_tumor_cor$value)&complete.cases(list_normal_cor$value),]

#pearson_cors<-cbind.data.frame(genepair=paste(list_tumor_cor2$X1,list_tumor_cor2$X2,sep='-'),normal=list_normal_cor2$value,tumor=list_tumor_cor2$value)
pearson_cors<-cbind.data.frame(tumor_name=list_tumor_cor2$X1,normal_name=list_tumor_cor2$X2,normal=list_normal_cor2$value,tumor=list_tumor_cor2$value)
write.csv(pearson_cors,'multinfo_all_cors.csv')

#p_order_normal_cors<-pearson_cors$normal[order(pearson_cors$normal,decreasing = TRUE)]
#p_order_normal_cors<-p_order_normal_cors[seq(1,length(p_order_normal_cors),5)]

#p_order_tumor_cors<-pearson_cors$tumor[order(pearson_cors$tumor,decreasing = TRUE)]
#p_order_tumor_cors<-p_order_tumor_cors[seq(1,length(p_order_tumor_cors),5)]

#p_order_cors<-cbind.data.frame(tumor=p_order_tumor_cors,normal=p_order_normal_cors)
#write.csv(p_order_cors,'multinfo_order_cors.csv')

