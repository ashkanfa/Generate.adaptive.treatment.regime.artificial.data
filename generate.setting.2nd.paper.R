###################################################################################################
#setwd("C:/Users/RandKID/OneDrive/Research/First Paper/Generate.simulated.cases/Simulate.Cases")


#rm(list=ls())
###################################################################################################
## Generating upper diagnol correlation vector by tru.cov.prop and true.trt.prop and
## three correlation structure of withing important (cor1), btw important and unimportant (cor2)
## and within unimportant (cor3)
#cor.vec.generate=function(n,p,cor1,cor2,cor3){
    #n: number of variables
    # p: proportion of important variables in n.
#    n1=n*p
#    n2=n*(1-p)
#    corr.vec<-list()
#    for(i in 1:(n1-1)){
#        corr.vec[[i]]<-c(rep(cor1,n1-i),rep(cor2,n2))
#    }
#    corr.vec[[n1]]<-rep(cor2,n2)
#    corr.vec[[n1+1]]<- rep(cor3,n2*(n2-1)/2)
#    return(unlist(corr.vec))
#}
#define necessary functions
low=function(no.low){
    round(runif(no.low,min=0,max=0.2),digits = 1)
}

med=function(no.med){
    round(runif(no.med,min=0.3,max=0.6),digits = 1)
}

high=function(no.high){
    round(runif(no.high,min=0.7,max=0.9),digits = 1)
}

cor.vec.generate=function(p, true.p, cor1="med",cor2="med",cor3="med"){
    #n: number of variables
    # p: fraction/proportion of important variables in n.
    
    
    p1=p*true.p     #no of true predictors
    p2=p*(1-true.p) #no of spurious predictors
    
    n1=p1*(p1-1)/2
    n2=p1*p2
    n3=p2*(p2-1)/2
    
    cor1.vec<-NULL
    cor2.vec<-NULL
    cor3.vec<-NULL
    cor.vec<-list()
    
    cor1.vec= switch(as.character(cor1),low=low(n1)
                     ,med=med(n1)
                     ,high=high(n1))
    
    cor2.vec= switch(as.character(cor2),low=low(n2)
                     ,med=med(n2)
                     ,high=high(n2))
    
    cor3.vec= switch(as.character(cor3),low=low(n3)
                     ,med=med(n3)
                     ,high=high(n3))
    
    for(i in 1:(p1-1)){
        
        cor.vec[[i]]<-c(cor1.vec[1:(p1-i)],cor2.vec[1:p2])
        cor1.vec<-cor1.vec[-(1:(p1-i))]
        cor2.vec<-cor2.vec[-(1:p2)]
    }
    cor.vec[[p1]]<-cor2.vec #if the code sorks correlctly the remaining elements in cor2.vec should be the right number of elements
    
    cor.vec[[(p1+1)]]<- cor3.vec
    
    if(length(unlist(cor.vec))!=p*(p-1)/2){
        stop("The generated correlation matrix has a misspecified dimension!")
    }
    
    return(unlist(cor.vec))
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

setwd("C:/Users/RandKID/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases")
setwd("C:/Users/RandKID/OneDrive/Research/second Paper")
############################################################
## Covariate Setting
#set.seed(1)
#no.bin=7; no.nor=1
#d=no.bin+no.nor
#mean.vec.nor=0; var.nor=1
#prop.vec.bin=round(runif(no.bin,min=0.2,max=0.8),digits = 1)
#true.cov.prop = 0.25
#cor1=0.7
#cor2=0.9
#cor3=0.9

#we assume that the correaltion structure of trts and covariates are the same, which is not realistic and someone can explore the what if scenarios 
# in which these two are different

case.1<-expand.grid(covar=16,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                   min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                   cor3="med",no.int.co=1,no.int.cot=0,no.int.t=2 )

case.2<-expand.grid(covar=16,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=1,no.int.t=2 )

case.3<-expand.grid(covar=16,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=2,no.int.cot=0,no.int.t=3 )

case.4<-expand.grid(covar=16,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=2,no.int.t=3 )
####
case.5<-expand.grid(covar=20,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=1,no.int.cot=0,no.int.t=2 )

case.6<-expand.grid(covar=20,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=1,no.int.t=2 )

case.7<-expand.grid(covar=20,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=2,no.int.cot=0,no.int.t=3 )

case.8<-expand.grid(covar=20,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=2,no.int.t=3 )
####
case.9<-expand.grid(covar=16,trt=12, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=1,no.int.cot=0,no.int.t=2 )

case.10<-expand.grid(covar=16,trt=12, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=1,no.int.t=2 )

case.11<-expand.grid(covar=16,trt=12, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=2,no.int.cot=0,no.int.t=3 )

case.12<-expand.grid(covar=16,trt=12, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                    min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                    cor3="med",no.int.co=0,no.int.cot=2,no.int.t=3 )

#no interaction 
case.13<-expand.grid(covar=16,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                     min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                     cor3="med",no.int.co=0,no.int.cot=0,no.int.t=0 )

case.14<-expand.grid(covar=20,trt=8, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
                     min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
                     cor3="med",no.int.co=0,no.int.cot=0,no.int.t=0 )

case.15<-expand.grid(covar=16,trt=12, n = 100, true.cov.prop=c(0.25,0.5), true.trt.prop=c(0.5,0.75),
            min.coef=0.6,max.coef=1,cor1=c("low","med","high"),cor2=c("low","med","high"), 
            cor3="med",no.int.co=0,no.int.cot=0,no.int.t=0 )

case<-rbind(case.1,case.2,case.3,case.4,case.5,case.6,case.7,case.8,case.9,case.10,case.11,case.12,case.13,case.14,case.15)


sim.data<-NULL
dat.1<-NULL
dat.2<-NULL
bootstrap<-NULL
big.data.stg1<-NULL
big.data.stg2<-NULL
for (i in seq(dim(case)[1])){
    #tryCatch({
#for (i in seq(10)){    
    #creating directory for the setting
    #setwd("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases")
    #dir.create(paste0(i))
    #setwd(paste0("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases/",i))
    #write.csv(case[i,],file="setting.csv",row.names = F)
    
    
    
    set.seed(i)
    no.rows=case$n[i]
    no.bin =(case$covar[i]) - 1
    no.nor = 1
    prop.vec.bin = round(runif(no.bin,min=0.2,max=0.8),digits = 1)
    mean.vec.nor=0
    var.nor=1
    #make this constant, this is not the focus of the study
    
    # interms of corr ver generation for covariates since second half of data is 
    # the one among which predicts the outcome true.p should always be 0.5
    corr.vec = cor.vec.generate(p=case$covar[i], true.p=0.5, cor1=case$cor3[i],cor2=case$cor2[i],cor3=case$cor1[i])
    #corr.vec = round(runif(case$covar[i]*(case$covar[i]-1)/2,min=0.3,max=0.6),digits = 1)
    no.trt=case$trt[i]
    no.stage=3
    bincorr.vec   =cor.vec.generate(p=case$trt[i], true.p=case$true.trt.prop[i], cor1=case$cor1[i],cor2=case$cor2[i],cor3=case$cor3[i])
    true.trt.prop =case$true.trt.prop[i]
    true.cov.prop =case$true.cov.prop[i]
    min.coef      =case$min.coef[i]
    max.coef      =case$max.coef[i]
    no.int.co     =case$no.int.co[i]
    no.int.cot    =case$no.int.cot[i]
    no.int.t      =case$no.int.t[i] 
     
    
    sim.data[[i]]<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, 
                                      corr.vec, no.trt,no.stage,bincorr.vec, true.trt.prop,true.cov.prop,
                                      min.coef, max.coef, no.int.co, no.int.cot, no.int.t)
    
    
        #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }



for(i in seq(dim(case)[1])){
    #creating directory for the setting
    setwd("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases")
    dir.create(paste0(i))
    setwd(paste0("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases/",i))
    write.csv(case[i,],file="setting.csv",row.names = F)
    
    dat.1[[i]]    <-cbind(sim.data[[i]][["cov"]][["cov"]],sim.data[[i]][["trt"]][["stg.0"]],sim.data[[i]][["trt"]][["stg.1"]],sim.data[[i]][["y"]][,1:2])
    dat.2[[i]]    <-cbind(sim.data[[i]][["cov"]][["cov"]],sim.data[[i]][["trt"]][["stg.0"]],sim.data[[i]][["trt"]][["stg.1"]], sim.data[[i]][["trt"]][["stg.2"]],sim.data[[i]][["y"]])
    
    bootstrap[[i]]<-matrix(sample(1:no.rows,no.rows*100,replace = T),nrow = no.rows,ncol=100)
    bootstrap[[i]]<-apply(bootstrap[[i]],MARGIN = 2,FUN = sort)
    
    data.stg.1<-NULL
    data.stg.2<-NULL
    
    for (j in seq(99)){
        
        data.stg.1[[j]]<-dat.1[[i]][bootstrap[[i]][,j],]
        data.stg.2[[j]]<-dat.2[[i]][bootstrap[[i]][,j],]
    }
    
    data.stg.1[[100]]<-dat.1[[i]]
    data.stg.2[[100]]<-dat.2[[i]]
    
    #cov.cor  <- cbind(sim.data[[i]][["cov"]][["upper corr vec"]],corr.vec)
    #trt.cor  <- cbind(sim.data[[i]][["trt"]][["uppervec.cor.stg.1"]],bincorr.vec)
    true.out.stg.1 <- sim.data[[i]][["True.outcome"]][["stg.1"]]
    true.out.stg.2 <- sim.data[[i]][["True.outcome"]][["stg.2"]]
    #data1[["True.outcome"]][["stg.1"]]
    
    for (i in 1:100){
        write.csv(data.stg.1[[i]],file=paste0("data.stg1.",i,".csv"),row.names = F)
        write.csv(data.stg.2[[i]],file=paste0("data.stg2.",i,".csv"),row.names = F)
    }
    
    write.csv(true.out.stg.1,file = "true.outcome.stg1.csv",row.names=F)
    write.csv(true.out.stg.2,file = "true.outcome.stg2.csv",row.names=F)
    #write.csv(cov.cor,file = "cov.cor.csv",row.names=F)
    #write.csv(trt.cor,file = "trt.cor.csv",row.names=F)
    
    #data<-NULL
    
    #data<-list(data.stg.1, data.stg.2)
    
    
    
    
    
    big.data.stg1[[i]]<-data.stg.1
    big.data.stg1[[i]]<-data.stg.2
}




for(i in seq(dim(case)[1])){
    #creating directory for the setting
    #setwd("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases")
    #dir.create(paste0(i))
    #setwd(paste0("C:/Users/randkid/OneDrive/Research/second Paper/Generate.simulated.cases/Simulate.Cases/",i))
    #write.csv(case[i,],file="setting.csv",row.names = F)
    
    #dat.1[[i]]    <-cbind(sim.data[[i]][["cov"]][["cov"]],sim.data[[i]][["trt"]][["stg.0"]],sim.data[[i]][["trt"]][["stg.1"]],sim.data[[i]][["y"]][,1:2])
    #dat.2[[i]]    <-cbind(sim.data[[i]][["cov"]][["cov"]],sim.data[[i]][["trt"]][["stg.0"]],sim.data[[i]][["trt"]][["stg.1"]], sim.data[[i]][["trt"]][["stg.2"]],sim.data[[i]][["y"]])
    
    bootstrap[[i]]<-matrix(sample(1:no.rows,no.rows*100,replace = T),nrow = no.rows,ncol=100)
    bootstrap[[i]]<-apply(bootstrap[[i]],MARGIN = 2,FUN = sort)
    
    data.stg.1<-NULL
    data.stg.2<-NULL
    
    for (j in seq(99)){
        
        data.stg.1[[j]]<-dat.1[[i]][bootstrap[[i]][,j],]
        data.stg.2[[j]]<-dat.2[[i]][bootstrap[[i]][,j],]
    }
    
    data.stg.1[[100]]<-dat.1[[i]]
    data.stg.2[[100]]<-dat.2[[i]]

    
    
    big.data.stg1[[i]]<-data.stg.1
    big.data.stg1[[i]]<-data.stg.2
}

