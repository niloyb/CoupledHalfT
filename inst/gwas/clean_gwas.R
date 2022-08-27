# Script to pre-process GWAS data

rm(list = ls())
# setwd('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/')
if (!require(pacman)) {install.packages('pacman')}

pacman::p_load(ARTP2,readr,tidyr,dplyr,R.matlab,parallel,glmnet,fastR)

cor.col <- function(j,x) {
  c <- cor(x[,j],x[,-j])
  return(which(c==1))
}

# df <- read.bed('Data/gwas/maize/imputemaize.bed', 'Data/gwas/maize/imputemaize.bim', 'Data/gwas/maize/imputemaize.fam',encode012 = TRUE)
df <- read.bed('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/Data/gwas/maize/imputemaize.bed', 
               '/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/Data/gwas/maize/imputemaize.bim', 
               '/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/Data/gwas/maize/imputemaize.fam',
               encode012 = TRUE)

nzero <- apply(df,2,function(x){sum(x==0)})
subject <- rownames(df)
variant <- names(df)
p <- ncol(df)

dat <- read_tsv('/Users/niloybiswas/Dropbox/horseshoe_coupling/Code/Data/gwas/maize/Romay_etal_2013_GenomeBiol_phenotypes-130503.txt')

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
# Normalizing design matrix
X <- scale(X)
# Using log(abs()) of response if all continuous responses same sign
if((mean(y>0)==1)|(mean(y>0)==0)){y <- log(abs(y))}
# Normalizing response
y <- scale(y)
gc()

# save(X, file="/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/processed_maize_data/design_matrix_Xnew.RData")
# save(y, file="/Users/niloybiswas/Google Drive/My Drive/Niloy_Files/Harvard/PhD/Research/high_dim_datasets/processed_maize_data/response_ynew.RData")