# Script to pre-process GWAS data

rm(list = ls())
# setwd('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/')
if (!require(pacman)) {install.packages('pacman')}

pacman::p_load(ARTP2,readr,tidyr,dplyr,R.matlab,parallel,glmnet,fastR)

cor.col <- function(j,x) {
  c <- cor(x[,j],x[,-j])
  return(which(c==1))
}

df <- read.bed('Data/gwas/maize/imputemaize.bed', 'Data/gwas/maize/imputemaize.bim', 'Data/gwas/maize/imputemaize.fam',encode012 = TRUE)

nzero <- apply(df,2,function(x){sum(x==0)})
subject <- rownames(df)
variant <- names(df)
p <- ncol(df)

dat <- read_tsv('Data/gwas/maize/Romay_etal_2013_GenomeBiol_phenotypes-130503.txt')

subject2 <- dat$Complete_name

subject1 <- subject

tmp <- regexpr('\\.',subject1)
tmp <- sapply(tmp,function(x){return(x[[1]])})

subject1 <- substr(subject1,tmp+1,nchar(subject1))
subject1 <- gsub('\\.',":",subject1)

tmp <- regexpr(':',subject2)
tmp <- sapply(tmp,function(x){return(x[[1]])})

subject2 <- substr(subject2,tmp+1,nchar(subject2))

common <- intersect(subject1,subject2)

dat$name <- subject2
dat <- dat[,names(dat) %in% c('name','GDD_DTS')]

df$name <- subject1
df <- df %>% inner_join(dat,by='name')
df <- df[!is.na(df$GDD_DTS),]
y <- df$GDD_DTS
X <- as.matrix(df[,1:p])
Xorig = X
X <- scale(X)
y <- scale(y)
gc()

# drops <- rep(F,p)
N <- nrow(X)

save(X, file="design_matrix_Xnew.RData")
save(y, file="response_ynew.RData")

# clean out correlation 1's
# tol = 1
# for (j in seq(from=1000,to=p,by=1000)) {
#   print(j)
#   # X0tX <- t(X[,((j-1000)+1):j])%*%X[,((j-1000)+1):p] 
#   X0tX <- crossprod(X[,((j-1000)+1):j],X[,((j-1000)+1):p])
#   ut <- upper.tri(X0tX,diag=F)
#   X0tX <- X0tX*(1*ut)
#   odmax <- apply(X0tX,2,max)
#   dropj <- rep(F,p)
#   dropj[((j-1000)+1):p] <- (abs(odmax/(N-1))>=tol)
#   drops <- drops | dropj
# }

#X <- Xorig[,!drops]

# writeMat(con='Data/gwas/maize/maize.mat',X=X,y=y)
# save(X,y,file='Data/gwas/maize/maize.RData')
# t1 <- proc.time()
# fit.lass <- cv.glmnet(X,y)
# t2 <- proc.time()
# delta.t <- t2-t1
# bets <- fit.lass$glmnet.fit$beta[,fit.lass$lambda==fit.lass$lambda.min]
# 
# save(file='Outputs/lasso_maize.RData',bets,fit.lass,delta.t)

# dat.horse <- readMat('Outputs/post_reg_horse_maize_approx_4_2267_98385.mat')
# which.horse <- dat.horse$keep.id[1:100]
# 
# which.common <- intersect(which.pos,which.horse)
# bets[which.common]
# 
# par(mfrow=c(1,2))
# hist(bets[which.common])
# hist(bets[bets>0])


# l.id <- as.list(seq(p))
# t1 <- proc.time()
# tmp <- cor.col(1,X)
# t2 <- proc.time()
# t2 - t1
# cor.all <- lapply(l.id,cor.col,x=X)




