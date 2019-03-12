###################################################################################################
setwd("c:/Users/axa8855/OneDrive/Research/Proposal/Generating Simulated Data/Generated_data/Proposal")

#input of the function
rm(list=ls())


set.seed(5)
no.rows=500
no.bin=19; no.nor=1
d=no.bin+no.nor
mean.vec.nor=0; var.nor=1
prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1)
corr.vec=c(round(runif(d*(d-1)/4,min=0.01,max=0.3),digits = 2),round(runif(d*(d-1)/4,min=0.5,max=0.7),digits = 2))
cmat = lower.tri.to.corr.mat(corr.vec,d)
no.trt=6
no.stage=3
# min and max of the correlation can be input of the data generating function
#bincorr.vec = round(runif(no.trt*(no.trt-1)/2,min=0.05,max=0.8),digits = 2)
d1=floor(median(1:(no.trt*(no.trt-1)/2)))
d2=(no.trt*(no.trt-1)/2)-d1
bincorr.vec = c(rep(0.4,d1),rep(0.7,d2))
bincorr = lower.tri.to.corr.mat(bincorr.vec,no.trt)

#sigma.star=compute.sigma.star(no.bin=19, no.nor=1, prop.vec.bin=round(runif(19,min=0.3,max=0.8),digits = 1),
#                              corr.mat=cmat)
#mydata=jointly.generate.binary.normal(no.rows,no.bin,no.nor,prop.vec.bin,
#                                     mean.vec.nor,var.nor, sigma_star=sigma.star$sigma_star,
#                                     continue.with.warning=TRUE)


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




a<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, corr.vec, no.trt,no.stage,bincorr.vec)
