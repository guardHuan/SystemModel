library(parallel)

tumor<-read.csv('../d_tumor.csv')
normal<-read.csv('../d_normal.csv')
rownames(tumor)<-tumor[,1]
rownames(normal)<-normal[,1]
tumor<-tumor[,-1]
normal<-normal[,-1]


n<-length(normal)
pair <- t(combn(1:n, 2))
nt <- detectCores()
cl <- makeCluster(type = "SOCK", rep("localhost", nt - 1)) 
ichunk <- split(as.data.frame(pair), 2:nt)


f1 <- function(id, x) {
  nr <- nrow(id)
  res <- matrix(0, nr, 1)
  k <- 1
  for (i in 1:nr) {
    res[k, ] <-energy::dcor(x[, id[i, 1]], x[, id[i, 2]])
    k <- k + 1
  }
  cbind(id, res)
}

ke_normal_cors <- clusterApply(cl, ichunk, f1, normal)
ke_tumor_cors <- clusterApply(cl, ichunk, f1, tumor)

stopCluster(cl)


list_tumor_cor<-data.frame()
for(i in 1:length(ke_tumor_cors)){
  
  list_tumor_cor<-rbind(list_tumor_cor,ke_tumor_cors[[i]])
}

list_normal_cor<-data.frame()
for(i in 1:length(ke_normal_cors)){
  
  list_normal_cor<-rbind(list_normal_cor,ke_normal_cors[[i]])
}

list_tumor_cor2<-list_tumor_cor[complete.cases(list_tumor_cor$res)&complete.cases(list_normal_cor$res),]
list_normal_cor2<-list_normal_cor[complete.cases(list_tumor_cor$res)&complete.cases(list_normal_cor$res),]

tu_name<-colnames(tumor[,list_tumor_cor2$V2])
no_name<-colnames(tumor[,list_normal_cor2$V1])
#pearson_cors<-cbind.data.frame(genepair=paste(no_name,tu_name,sep='-'),normal=list_normal_cor2$res,tumor=list_tumor_cor2$res)
pearson_cors<-cbind.data.frame(tumor_name=tu_name,normal_name=no_name,normal=list_normal_cor2$res,tumor=list_tumor_cor2$res)
write.csv(pearson_cors,'dcor_all_cors_new.csv')

#p_order_normal_cors<-pearson_cors$normal[order(pearson_cors$normal,decreasing = TRUE)]
#p_order_normal_cors<-p_order_normal_cors[seq(1,length(p_order_normal_cors),5)]

#p_order_tumor_cors<-pearson_cors$tumor[order(pearson_cors$tumor,decreasing = TRUE)]
#p_order_tumor_cors<-p_order_tumor_cors[seq(1,length(p_order_tumor_cors),5)]

#p_order_cors<-cbind.data.frame(tumor=p_order_tumor_cors,normal=p_order_normal_cors)
#write.csv(p_order_cors,'dcor_order_cors.csv')


