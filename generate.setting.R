###################################################################################################
#setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")

rm(list=ls())
###################################################################################################
## Generating upper diagnol correlation vector by tru.cov.prop and true.trt.prop and
## three correlation structure of withing important (cor1), btw important and unimportant (cor2)
## and within unimportant (cor3)
cor.vec.generate=function(n,p,cor1,cor2,cor3){
    #n: number of variables
    # p: proportion of important variables in n.
    n1=n*p
    n2=n*(1-p)
    corr.vec<-list()
    for(i in 1:(n1-1)){
        corr.vec[[i]]<-c(rep(cor1,n1-i),rep(cor2,n2))
    }
    corr.vec[[n1]]<-rep(cor2,n2)
    corr.vec[[n1+1]]<- rep(cor3,n2*(n2-1)/2)
    return(unlist(corr.vec))
}
###################################################################################################

library("devtools")
library("githubinstall")
library("mvtnorm")
library("corpcor")
library("psych")
library("Matrix")
library("BinNor")
library("ICC")
library("miscTools")
library("car")
library("plyr")
library("MASS")
library("gtools")
library("nnet")
library("geepack")
library("bindata")
library("splitstackshape")
library("boot")


##################################################################################################
##################################################################################################
txt <-" 1- First Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.75, True.cov.prop = 0.25\n
Correlation structure for trt :\n 
Within important medium-high btw medium-high \n 
for correlation of Covariates medium-high\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("1")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/1")
writeLines(txt,"setting.txt")
############################################################
## Covariate Setting
set.seed(1)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.9
cor3=0.9
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


#corr.vec<-vector()
corr.vec<-rep(0.8,d*(d-1)/2)
#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.75
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(1)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                          no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){

data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)


###########################################################################################################
###########################################################################################################
txt <-" 2- second Setting\n
 P=20 d=8 T=12\n 
 True.trt.prop = 0.25, True.cov.prop = 0.25\n
 Within important 0 btw 0 within unimportant 0\n 
for correlation of Covariates 0.5\n 
 n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("2")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/2")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(2)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0
cor2=0
cor3=0
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(2)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)


###########################################################################################################
###########################################################################################################
txt <-" 3- Third Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0 btw 0 within unimportant 0\n 
for correlation of Covariates 0.8\n 
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("3")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/3")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(3)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0
cor2=0
cor3=0
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(3)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)


###########################################################################################################
###########################################################################################################
txt <-" 4- 4th Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.4 btw 0.7 within unimportant 0.7\n
Cov Corr = 0.5\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("4")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/4")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(4)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.4
cor2=0.7
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(4)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)

###########################################################################################################
###########################################################################################################
txt <-" 5- 5th Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.4 btw 0.7 within unimportant 0.7\n
Cov Corr = 0.8\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("5")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/5")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(5)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.4
cor2=0.7
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(5)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
###########################################################################################################
###########################################################################################################
txt <-" 6- 6th Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.0.7 btw 0.4 within unimportant 0.4\n
Cov Corr = 0.5\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("6")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/6")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(6)

no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.4
cor3=0.4
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(6)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 7- 7th Setting\n
 P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.7 btw 0.4 within unimportant 0.4\n
Cov Corr = 0.8\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("7")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/7")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(7)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.4
cor3=0.4
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(7)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 8- 8th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.7 btw 0.4 within unimportant 0.7\n
Cov Corr = 0.5\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("8")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/8")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(8)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.4
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(8)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 9- 9th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.7 btw 0.4 within unimportant 0.7\n
Cov Corr = 0.8\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("9")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/9")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(9)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.4
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(9)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 10- 10th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.4 btw 0.7 within unimportant 0.4\n
Cov Corr = 0.5\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("10")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/10")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(10)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.4
cor2=0.7
cor3=0.4
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(10)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 11- 11th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.4 btw 0.7 within unimportant 0.4\n
Cov Corr = 0.8\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("11")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/11")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(11)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.4
cor2=0.7
cor3=0.4
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(11)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 12- 12th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.7 btw 0.7 within unimportant 0.7\n
Cov Corr = 0.5\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("12")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/12")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(12)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.7
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.5,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(12)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
#####################################################################################################
#####################################################################################################
txt <-" 13- 13th Setting\n
P=20 d=8 T=12\n 
True.trt.prop = 0.25, True.cov.prop = 0.25\n
Within important 0.7 btw 0.7 within unimportant 0.7\n
Cov Corr = 0.8\n
n=500"
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")
dir.create("13")
setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases/13")
writeLines(txt, "setting.txt")
############################################################
## Covariate Setting
set.seed(13)
no.bin=7; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
true.cov.prop = 0.25
cor1=0.7
cor2=0.7
cor3=0.7
#no.true.cov = floor(quantile(1:(d/4),2*true.cov.prop))
#true.cov.ot.index    <-floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.5)+no.true.cov) #Outcome and treatments
#true.cov.o.index     <-floor(quantile(1:d,.75)+1) : floor(quantile(1:d,.75)+no.true.cov)                       #only outcome
#n1=d/2


corr.vec<-rep(0.8,d*(d-1)/2)

#corr.vec<-cor.vec.generate(n=d,p=true.cov.prop,cor1,cor2,cor3)
#cor1<- rep(0.5,n1*(n1-1)/2)   # within Unimporatnt (Cov.1 and cov.2)
#cor2<- rep(0.8,n1**2)         # Between Imporatant and unimportant
#cor3<- rep(0.5,n1*(n1-1)/2)   # Within important (Cov.3 and cov.4)
#corr.vec=c(cor1,cor2,cor3)
#cmat = lower.tri.to.corr.mat(corr.vec,d)
############################################################
## Trt Setting

no.trt=12
true.trt.prop=0.25
bincorr.vec<-cor.vec.generate(n=no.trt,p=true.trt.prop,cor1,cor2,cor3)
#n2=no.trt/2
#cort1<- rep(0.5,n2*(n2-1)/2)   # within imporatnt trt 
#cort2<- rep(0.8,n2**2)         # Between Imporatant and unimportant
#cort3<- rep(0.5,n2*(n2-1)/2)   # within unimportant 
#bincorr.vec <-c(cort1,cort2,cort3)
#bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)
#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)

############################################################
## Rest of the setting
no.stage=2            # number of time slices
no.rows=500           # no of Observation

## not suprising that not every seed works with this code. So don't get dissapointed if you are getting
## Error in dimnames(x) <- dn : Try with different seeds until you get results.
set.seed(13)
data1<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, 
                              no.trt, no.stage, bincorr.vec, true.trt.prop, true.cov.prop)


dat<-cbind(data1[["cov"]][["cov"]],data1[["trt"]][["stg.0"]],data1[["trt"]][["stg.1"]],data1[["y"]])
bootstrap<-matrix(sample(1:500,50000,replace = T),nrow = 500,ncol=100)
bootstrap<-apply(bootstrap,MARGIN = 2,FUN = sort)

data<-list()
for (i in 1:99){
    
    data[[i]]<-dat[bootstrap[,i],]
}

data[[100]]<-dat

cov.cor<-cbind(data1[["cov"]][["upper corr vec"]],corr.vec)
trt.cor<-cbind(data1[["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
true.out<-data1[["True.outcome"]][["stg.1"]]
#data1[["True.outcome"]][["stg.1"]]

for (i in 1:100){
    write.csv(data[[i]],file=paste0("data",i,".csv"),row.names = F)
}

write.csv(true.out,file = "true.outcome.csv",row.names=F)
write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
