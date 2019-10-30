
library(MASS)
library(gtools)
library(nnet)
library(geepack)
library(bindata)
library(splitstackshape)
library(BinNor)
library("devtools")
library("githubinstall")

generate.adaptive.data = function(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, 
                                  corr.vec, no.trt,no.stage,bincorr.vec, true.trt.prop,true.cov.prop,
                                  min.coef=0.4,max.coef=0.7){
  #true.trt.prop: is the proportion of the true trts in true outcome model to the total numberof trts
  #true.cov.prop: is the proportion of the true predictors to the total number of covariates
  
  
  
#############################################################################################################################################################

sum.row<-function(x){ 
    v<-matrix(nrow = dim(x)[1],ncol=1)
    for (i in 1:dim(x)[1])
        v[i,]=sum(x[i,])
    return(v)
}

rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}
#############################################################################################################################################################
if ((no.bin+no.nor)%%4!=0){stop("Sum of the Covariates should be multiplies of 4!\n")}
if ((no.bin+no.nor)==4 & true.cov.prop<0.5){stop("when Sum of the Covariates is 4, the true.cov.prop should be 0.5!\n")}
if (true.cov.prop > 0.5)  {stop("The maximum allowable true poportion of covariate predictors is 0.5!\n")}
#############################################################################################################################################################

slot_0 = mod.jointly.generate.binary.normal( no.rows, no.bin, no.nor, prop.vec.bin, mean.vec.nor, var.nor, corr.vec)
#no.rows=Number of subjects
#no.bin=Number of binary variables (trt vars and non-trt covariates e.g. sex )
#no.nor=Number of normally distributed variables
#prop.vec.bin=Vector of marginal proportions for binary variables
#mean.vec.nor=Vector of means for normal variables
#var.nor=Vector of variances for normal variables
#corr.vec=Specified correlations among all variables

d=no.bin+no.nor
cov.data<-slot_0[[1]]
colnames(cov.data) <- paste0("x", 1:d)


# covariates are divided into 4 categories.
# 1- Spurious (Irrelevant covariates),           cov.1
# 2- Only predictors of trts                     cov.2
# 3- Predictors of Both outcome and trt          cov.3
# 4- Only predictors of outcome                  cov.4

cov.1    <-cov.data [ , 1:floor(quantile(1:d,.25))]                           #Spurious
cov.2    <-cov.data [ , floor(quantile(1:d,.25)+1) : floor(median(1:d))]      #Trts
cov.3    <-cov.data [ , floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.75))] #Outcome and treatments
cov.4    <-cov.data [ , floor(quantile(1:d,.75)+1) : d]                       #only outcome


#defining marginal probability of 6 trts that are correlated and have
#1-half them contribute to outcome that means y0-> t -> y1
#2-and another half not                       y0-> t

#outcome should depend on cov.3 and cov.4
# trts should only depend on cov.2 and cov.3 (x6-x10 x11-x15)


# no.trt can be an input of a function
# no.trt is the total number of trts, later on half them will linearly contribute to response
# and second half will not.


#############################################################################################################################################################
## Defining Response

y          <-matrix(nrow = no.rows, ncol = no.stage)
colnames(y)<-paste0("y.",0:(no.stage-1))


## coef.y is  the coefficient of y in linear association with trt of the next stage.
## please note that this is different from coefficient of y in previous stage in linear association with y of the next stage
coef.y.trt  <- list()
for (i in 1:no.trt){
    coef.y.trt[[paste("trt.",i,sep="")]] = matrix(nrow = 1 ,ncol = no.stage-1)
    for (j in 1:(no.stage-1)){
        coef.y.trt[[paste("trt.",i,sep="")]][,j] = round(runif(1,min = -1,max=1),digits = 1)
    }

}



#pre-evaluation y
y[,1]   <- rnorm(n=no.rows,0,1)
#############################################################################################################################################################
## Define Trt list and the coefficients

trt<-list()
for (i in 1:no.stage){
    trt[[paste("stg.",i-1,sep = "")]]           <-matrix(nrow = no.rows, ncol=no.trt)
    #colnames(trt.stg.0)<-paste("trt.0.",1:no.trt,sep = "")
    trt[[paste("p.",i-1,sep = "")]]             <-matrix(nrow = no.rows, ncol=no.trt)
    trt[[paste("marg.p.",i-1,sep = "")]]        <-vector()
    trt[[paste("cor.trt.stg.",i-1,sep = "")]]   <-matrix(nrow = no.trt, ncol=no.trt)
    trt[[paste("coef.trt.stg.",i-1,sep = "")]]  <-matrix(nrow = no.trt, ncol=no.trt)
}


## Pre evaluation trt
for (i in 1:no.trt){
    trt[["marg.p.0"]][i]<- round(runif(1,min = 0.1,max = 0.8),1)
    trt[["stg.0"]] [,i] <- rbinom(no.rows,1,trt[["marg.p.0"]][i])
    colnames(trt[["stg.0"]]) <- paste("t", 1:no.trt ,".stg.0",sep = "")
    trt[["cor.trt.stg.0"]] <- round(cor(trt[["stg.0"]]),1)
}

## Define coefficient of the treatment depending on the treatments of the previous stage.


for (i in 1:(no.stage-1)){
    for (j in 1:no.trt){
        trt[[paste("coef.trt.stg.",i-1,sep = "")]][j,] = round(runif(no.trt,min = -0.9,max=0.9),digits = 1)
    }
    #coef.trt.stage[[i]]= coef.trt.trt
}


#############################################################################################################################################################
## Creating dependency of trts to linear combination of covariates
## The coefficients created and the linear ombination will remain the same through all stages
## meaning that dependenct of trts to the covariates will not change

## d.cov.trt is the maximum number of covariates that linearly will predict the trts
#cov.trt.all<-cbind(cov.2,cov.3)
d.cov.trt<-(dim(cov.2)[2]+dim(cov.3)[2])
## we set min.num.cov as the minimum number of covariates that linearly will contribute to the trts
## for now we set it to be %50 of the maximum number of covariates
## This can be put as a fraction input in the function later
min.num.cov=floor(quantile(1:d.cov.trt,0.5))
no.select.cov.trt = sample(min.num.cov : d.cov.trt, no.trt, replace = TRUE)
#trt-outcome predictors should be in the treatment linear model
no.select.cov.grp1.y = floor(quantile(1:dim(cov.3)[2],2*true.cov.prop))
which.cov.grp1.y = dim(cov.2)[2]+sample(1:dim(cov.3)[2], no.select.cov.grp1.y, replace = FALSE)




select.cov.trt<-list()
name.cov.trt<-list()
cov.trt<-list()
coef.cov.trt<-list()
no.coef<-vector()

for (i in 1:no.trt){
    select.cov.trt[[i]]    = sample(1:d.cov.trt, no.select.cov.trt[i], replace = FALSE)
    select.cov.trt[[i]]    = c(select.cov.trt[[i]],which.cov.grp1.y)
    select.cov.trt[[i]]    = unique(select.cov.trt[[i]])
    select.cov.trt[[i]]    = dim(cov.1)[2]+select.cov.trt[[i]]
    no.coef[i]=length(select.cov.trt[[i]])
    name.cov.trt[[i]]= paste("x",select.cov.trt[[i]],sep="")
    cov.trt[[i]] =cov.data[,name.cov.trt[[i]]]
    coef.cov.trt[[i]] = round(runif(no.coef[i],min = -0.85,max=0.90),digits = 1)#positive and negative coef
    coef.cov.trt[[i]][coef.cov.trt[[i]]==0]=0.5 # replace zero coef with 0.5 
}

#############################################################################################################################################################
## make the pool of predictors of outcome:  cov.3 and cov.4
## We create linear dependency of y on twice the true.cov.prop of cov.3 and
## and twice the true.cov.prop of cov.4, so that in total we have the true.cov.prop of all covariates 
## predicting the outcome Y.


#dim.cov.y=dim(cov.3)[2]

# number of selected covariated from each covariate group of 3 or 4
no.select.cov.grp1.y = floor(quantile(1:dim(cov.3)[2],2*true.cov.prop))
no.select.cov.grp2.y = floor(quantile(1:dim(cov.4)[2],2*true.cov.prop))
which.cov.grp1.y = sample(1:dim(cov.3)[2], no.select.cov.grp1.y, replace = FALSE)
which.cov.grp2.y = sample(1:dim(cov.4)[2], no.select.cov.grp2.y, replace = FALSE)


## keeping one normal(continuous) cov (last cov from cov.4)
#  as a predictor of y
last.var<-dim(cov.4)[2]
#which.cov.grp2.last.y <- ifelse(last.var %in% which.cov.grp2.y , which.cov.grp2.y , c(which.cov.grp2.y,last.var))
if (last.var %in% which.cov.grp2.y){
    which.cov.grp2.y
} else {which.cov.grp2.y<-c(which.cov.grp2.y,last.var) }
#cov.3[,which.cov.grp1.y]
#cov.4[,which.cov.grp2.y]

## Define a list that contains trts as predictors of y "trt.y"
trt.y<-list()
for (i in 1:(no.stage-1)){
    trt.y[[paste("stg.",i, sep="")]] = matrix(nrow = no.rows,ncol =length(1: floor(quantile(1:no.trt, true.trt.prop))))
}

#############################################################################################################################################################
## BINDING COVARIATE-TREATMENT OF THE PREVIOUS STAGE AND RESPONSE OF THE PREVIOUS STAGE and
## name it t.state

t.state<-list()
y.state<-list()
w<-list()
dim.y<-vector()
#index<-matrix(nrow = no.stage-1 , ncol = no.trt)

for (i in 1:(no.stage-1)){
    for (j in 1:no.trt){
        
 #no.trt.dep[i,]<- sample(1:(no.trt/4),no.trt,replace = T)
 true.trt= true.trt.prop*no.trt
 no.trt.dep1<-matrix(nrow = no.stage-1,ncol =true.trt )
 no.trt.dep1[i,]<- sample(1:true.trt,true.trt,replace = T)
 no.trt.dep2<-matrix(nrow = no.stage-1,ncol =no.trt-true.trt )
 no.trt.dep2[i,]<- sample(1:(no.trt/4), no.trt-true.trt, replace = T)
 no.trt.dep<-cbind(no.trt.dep1,no.trt.dep2)
 
 if (j <= true.trt){
     
     #no.trt.dep1[i,]<- sample(1:true.trt,true.trt,replace = T)
     
     t.state[[paste("stg",i,".trt",j,sep = "")]] = cbind (cov.trt[[j]] , trt[[paste("stg.",i-1,sep = "")]][,1:true.trt], y[,paste("y.",i-1,sep = "")])
     t.state[[paste("coef.stg",i,".trt",j,sep = "")]] = c(coef.cov.trt[[j]], trt[[paste("coef.trt.stg.",i-1,sep = "")]][j,][1:true.trt],
                                                          coef.y.trt[[paste("trt.",j,sep="")]][,i])
 }
 
 else {
     #no.trt.dep2<-matrix(nrow = no.stage-1,ncol =no.trt-true.trt )
     #no.trt.dep2[i,]<- sample(1:(no.trt/2), no.trt-true.trt, replace = T)
     
     t.state[[paste("stg",i,".trt",j,sep = "")]] = cbind(cov.trt[[j]] , trt[[paste("stg.",i-1,sep = "")]][,(true.trt+1) : (true.trt+1+no.trt.dep[i,j])], y[,paste("y.",i-1,sep = "")])
     t.state[[paste("coef.stg",i,".trt",j,sep = "")]] = c(coef.cov.trt[[j]], trt[[paste("coef.trt.stg.",i-1,sep = "")]][j,][(true.trt+1): (true.trt+1+no.trt.dep[i,j])],
                                                          coef.y.trt[[paste("trt.",j,sep="")]][,i])
 }
# t.state is the predictors of each trt per stage
 #t.state[[paste("stg",i,".trt",j,sep = "")]] = cbind(cov.trt[[j]] , trt[[paste("stg.",i-1,sep = "")]][,1:no.trt.dep[i,j]], y[,paste("y.",i-1,sep = "")])
 
 #index[i,j]=dim(t.state[[paste("stg",i,".trt",j,sep = "")]])[2]
 
 colnames(t.state[[paste("stg",i,".trt",j,sep = "")]])[ncol(t.state[[paste("stg",i,".trt",j,sep = "")]])]<-paste("y.",i-1,sep = "")
 
 #t.state[[paste("coef.stg",i,".trt",j,sep = "")]] = c(coef.cov.trt[[j]], trt[[paste("coef.trt.stg.",i-1,sep = "")]][j,][1:no.trt.dep[i,j]],
  #                                                   coef.y.trt[[paste("trt.",j,sep="")]][,i])

 w[[paste("stg",i,".trt",j,sep = "")]]=rep.row(t.state[[paste("coef.stg",i,".trt",j,sep = "")]],no.rows)*t.state[[paste("stg",i,".trt",j,sep = "")]]




 trt[[paste("p.",i,sep = "")]][,j]=round((exp(sum.row( w[[paste("stg",i,".trt",j,sep = "")]] ))) / (1+exp(sum.row( w[[paste("stg",i,".trt",j,sep = "")]] ))),digits = 1)
 ###### because for creating a multivariate binomial we should only have 1 marginal per trt,
 ###### and not 1 marginal per obs, there should be a way to get the expectation of observational
 ###### p[i,] and have a truly marginal instead of mean.
 trt[[paste("marg.p.",i,sep = "")]][j] = round(mean(trt[[paste("p.",i,sep = "")]][,j]),digits = 1)

    }
 #if (j<no.trt){next}
 trt.list <- mod.generate.jointly.binary(no.rows,no.binary=no.trt,prop.vec.binary=trt[[paste("marg.p.",i,sep = "")]],
                                         corr.vec.binary=bincorr.vec, adjust.corrs = TRUE)
 trt[[paste("stg.",i,sep = "")]] <-trt.list[[1]]
 colnames(trt[[paste("stg.",i,sep = "")]]) <- paste("t", 1:no.trt ,".stg.",i,sep = "")
 trt[[paste("cor.trt.stg.",i,sep = "")]] <-round(cor(trt.list[[1]]),1)
 trt[[paste("uppervec.cor.stg.",i,sep="")]]<- upper_tri_vec(round(cor(trt.list[[1]]),1))

 #           Deviding Trts to two groups
 #            1- Trts which will linearly predict the outcome and their proportion is defined by true.trt.prop
 #            2-      Trts that will not affect the "outcome"
 trt.y[[paste("stg.",i, sep="")]]      <-trt[[paste("stg.",i,sep = "")]][ , 1: floor(quantile(1:no.trt, true.trt.prop))]

 y.state[[paste("stg.",i, sep="")]]  = cbind(cov.3[,which.cov.grp1.y],cov.4[,which.cov.grp2.y], trt.y[[paste("stg.",i, sep="")]], y[,paste("y.",i-1,sep = "")])
 
 if (length(which.cov.grp1.y)<2){
     colnames(y.state[[paste("stg.",i, sep="")]])[1]=paste0("x",dim(cov.1)[2]+dim(cov.2)[2]+which.cov.grp1.y)
 }
 
 if (length(which.cov.grp2.y)<2){
     colnames(y.state[[paste("stg.",i, sep="")]])[length(which.cov.grp1.y)+1]=paste0("x",dim(cov.1)[2]+dim(cov.2)[2]+dim(cov.3)[2]+which.cov.grp1.y)
 }
 
 colnames(y.state[[paste("stg.",i, sep="")]])[ncol(y.state[[paste("stg.",i, sep="")]])] = paste("y.",i-1,sep = "")
 dim.y[i]  =  dim(y.state[[paste("stg.",i, sep="")]])[2]
 y.state[[paste("coef.stg.",i, sep="")]] = round(runif(dim.y[i],min = min.coef, max = max.coef),digits = 1)
 ## replacing 0 coef with 0.2
 y.state[[paste("coef.stg.",i, sep="")]][y.state[[paste("coef.stg.",i, sep="")]]==0]=(min.coef+max.coef)/2
 
 y[,paste("y.",i,sep = "")]    = sum.row((rep.row(y.state[[paste("coef.stg.",i, sep="")]],no.rows))*y.state[[paste("stg.",i, sep="")]]) # no error
 y[,paste("y.",i,sep = "")]    = y[,paste("y.",i,sep = "")]+rnorm(no.rows, 0, 1)

}

# Return the data
True.outcome   <-list()
True.treatment <-list()

for(i in 1:(no.stage-1)){
    True.outcome[[paste("stg.",i,sep = "")]]   <-cbind(colnames(y.state[[paste("stg.",i,sep = "")]]),y.state[[paste("coef.stg.",i,sep = "")]])
    #True.treatment[[paste("stg.",i,sep = "")]] <-
}
#for(i in 1:(no.stage-1)){
#for (j in (1:no.trt)){
#True.treatment[[paste("stg",i,".trt",j,sep = "")]]=cbind(colnames(t.state[[paste("stg",i,".trt",j,sep = "")]]),t.state[[paste("coef.stg",i,".trt",j,sep = "")]])

#}
#}
cov.data_wcor<-list("cov"=cov.data,"cov correlation"=round(cor(cov.data),1), "upper corr vec" = upper_tri_vec(round(cor(cov.data),1)))

return(list("cov"=cov.data_wcor,"trt"=trt ,"y"=y,"t.state"=t.state,"True.outcome"=True.outcome))


#trt for IPTW (MIMIC)
#write.csv(trt[["stg.1"]],file = "Stage1_AllNewTreatments_Binary.csv",row.names=FALSE)
#write.csv(trt[["stg.2"]],file = "Stage2_AllNewTreatments_Binary.csv",row.names=FALSE)


#cov for IPTW (MIMIC)
#write.csv(cbind(cov.data,trt[["stg.0"]],y[,1]),file = "NewCombinedbinaryconfoundingvar_stage1.csv",row.names=FALSE)
#write.csv(cbind(cov.data,trt[["stg.1"]],y[,2]),file = "NewCombinedbinaryconfoundingvar_stage2.csv",row.names=FALSE)


#simulated data set at each stage
#data_stage1<-cbind(cov.data,trt[["stg.0"]],trt[["stg.1"]],y[,1],y[,2])
#colnames(data_stage1)[(ncol(data_stage1)-1):ncol(data_stage1)] = c("y.0","y.1")

#data_stage2<-cbind(cov.data,trt[["stg.1"]],trt[["stg.2"]],y[,2],y[,3])
#colnames(data_stage2)[(ncol(data_stage2)-1):ncol(data_stage2)] <- c("y.1","y.2")

#write.csv(data_stage1,file = "data_stage1.csv",row.names=FALSE)
#write.csv(data_stage2,file = "data_stage2.csv",row.names=FALSE)

# True outcome model at each stage
#write.csv(cbind(colnames(y.state[["stg.1"]]),y.state[["coef.stg.1"]]),file="true_outcome_model_stg1.csv")
#write.csv(cbind(colnames(y.state[["stg.2"]]),y.state[["coef.stg.2"]]),file="true_outcome_model_stg2.csv")


#write.csv(y,file = "y.csv",row.names=FALSE)
}
