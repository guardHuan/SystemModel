tumor<-read.csv('../d_tumor.csv')
normal<-read.csv('../d_normal.csv')
rownames(tumor)<-tumor[,1]
rownames(normal)<-normal[,1]
tumor<-tumor[,-1]
normal<-normal[,-1]

########################################################Covr function
Covr <- function(x, method="HS", h.method="dpill", plot.true=FALSE, ...){
 
   require(princurve)
   require(KernSmooth)

# Computing the Principal Curve
    PC <- switch(method,
             hs=, HS = principal.curve(x,plot.true=plot.true,...),
             pcop=, PCOP = pcop(x,plot.true=plot.true,...)$pcop.f2
          )

# Calculating plug-in bandwidths   
##   h1 <- dpill(PC$lambda,PC$s[,1]) # dpill does not work for x and y
##   h2 <- dpill(PC$lambda,PC$s[,2]) # arrays equally ordered


# Estimating derivatives of the PC
#   "locpoly(...)" does not admit a xgrid array of points
#   where evaluating the regression function.
#   I use my "locpolreg(...)" function instead.
   h1 <- diff(range(PC$lambda))/10
   h2 <- 2*max(diff(sort(PC$lambda)))
   hlambda <- max(h1,h2) 
   lpr1 <- locpolreg_no_plot(PC$lambda,PC$s[,1],p=2,r=1,h=hlambda)
   lpr2 <- locpolreg_no_plot(PC$lambda,PC$s[,2],p=2,r=1,h=hlambda)
   # (Note: lpr*%x and lpr*$y are sorted according to lpr*$x)


# Computing angles alpha(s)
   alpha <- atan2(lpr1$mxgr, lpr2$mxgr)
   aux <- sort(PC$tag, ind=T)
   # recovering original ordering
#points(PC$lambda, lpr2$mxgr[aux$ix],col=2)

   alpha <- alpha[aux$ix]
 
# Estimating conditional orthogonal variance.
#  e2 <- apply((x - PC$s)^2,1,sum)
   nn <- dim(x)[1]
   e2 <- array(apply((x - PC$s)^2,1,sum),dim=c(nn,1))[,1]
   if (h.method != "dpill"){
      he2_1 <- diff(range(PC$lambda))/5
   }else{
      he2 <- dpill(PC$lambda, e2,trim=.01)
      cte.h <- 2.214 # h1 and h2 are computed for Gaussian kernel, and
                     # locpolreg uses Epanechnikov kernel
      he2_1 <- he2*cte.h
      # prevention against NaN produced by dpill
      if (is.na(he2_1)) he2_1 <- diff(range(PC$lambda))/5
   }   
   he2_2 <- 2*max(diff(sort(PC$lambda)))
   he2 <- max(he2_1,he2_2) 
   CondOrthVar <- locpolreg_no_plot(PC$lambda, e2,h=he2)$mxgr
   # recovering original ordering and preventing for negative estimations
   CondOrthVar <- pmax(CondOrthVar[aux$ix],0) 

# Var. over the PC
   VS <- var(PC$lambda)

# Computing local varinaces, covariance and correlation
   LV1 <- VS * cos(alpha)^2 + CondOrthVar*sin(alpha)^2
   LV2 <- VS * sin(alpha)^2 + CondOrthVar*cos(alpha)^2
   LCovGC <- (VS - CondOrthVar)*sin(alpha)*cos(alpha)
   LCorGC <- LCovGC/sqrt(LV1*LV2)

# Computing Generalized Covariance and Correlation along GC
   CovGC <- sqrt(mean(LCovGC^2))
   CorGC <- sqrt(mean(LCorGC^2))

# Value:
   return(list(PC=PC,CovGC=CovGC,CorGC=CorGC,LCovGC=LCovGC,LCorGC=LCorGC))
}

locpolreg_no_plot <- function(x,y,h=(max(x)-min(x))/5,p=1,r=0,xg=x,...){
   if (sum(diff(xg)<0)>0) xg <- sort(xg)
      
   n <- length(x);
   m <- length(xg);
   mxgr <- seq(1,m)*0;
   S <- matrix(0,nrow=m,ncol=n)

   factr <- max(1,prod(1:r))

   for (i in seq(1,m)){
      Ih <- abs(x-xg[i])<h;
      n <- sum(Ih);     
      xh <- x[Ih]-xg[i];
      Dp <- matrix(1,nrow=n,ncol=p+1);
      if (p>0){for (j in 1:p) Dp[,j+1] <- xh^j}
      Wx <- epan(xh/h)/h;
      Wm <- Wx%*%ones(1,p+1);
      Dpp <- Wm*Dp;
      Si <- solve(t(Dp)%*%Dpp)%*%t(Dpp);
      beta <- Si%*%y[Ih];
      mxgr[i] <- factr*beta[r+1];
      S[i,Ih] <- Si[r+1,]
   }
  
return(list(mxgr=mxgr,S=S))
}

epan <- function(x){.75*(x+1)*(1-x)}
ones <- function(n,m){matrix(1,nrow=n,ncol=m)}
##########################################################

cols<-ncol(tumor)
normal_gene_cor<-matrix(nrow = cols,ncol = cols)
colnames(normal_gene_cor)<-colnames(tumor)
rownames(normal_gene_cor)<-colnames(tumor) 
tumor_gene_cor<-matrix(nrow = cols,ncol = cols)
colnames(tumor_gene_cor)<-colnames(tumor)
rownames(tumor_gene_cor)<-colnames(tumor)
for (i in 1:ncol(normal_gene_cor)){
  for (j in 1:i){
    result = tryCatch({normal_gene_cor[i,j] <- Covr(as.matrix(data.frame(normal[,i],normal[,j])),method="HS")$CorGC},
                      error=function(e){
                        #cat("ERROR :",conditionMessage(e),"\n")
                        #cat(i,j,'\n','disease_error.txt')
                        normal_gene_cor[i,j] <- NA})
    normal_gene_cor[j,i] <- normal_gene_cor[i,j]
  }
}

for (i in 1:ncol(tumor_gene_cor)){
  for (j in 1:i){
    result = tryCatch({tumor_gene_cor[i,j] <- Covr(as.matrix(data.frame(tumor[,i],tumor[,j])),method="HS")$CorGC},
                      error=function(e){
                        #cat("ERROR :",conditionMessage(e),"\n")
                        #cat(i,j,'\n','disease_error.txt')
                        tumor_gene_cor[i,j] <- NA})
    tumor_gene_cor[j,i] <- tumor_gene_cor[i,j]
  }
}


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
write.csv(pearson_cors,'corgc_all_cors.csv')

#p_order_normal_cors<-pearson_cors$normal[order(pearson_cors$normal,decreasing = TRUE)]
#p_order_normal_cors<-p_order_normal_cors[seq(1,length(p_order_normal_cors),5)]

#p_order_tumor_cors<-pearson_cors$tumor[order(pearson_cors$tumor,decreasing = TRUE)]
#p_order_tumor_cors<-p_order_tumor_cors[seq(1,length(p_order_tumor_cors),5)]

#p_order_cors<-cbind.data.frame(tumor=p_order_tumor_cors,normal=p_order_normal_cors)
#write.csv(p_order_cors,'sp_order_cors.csv')
