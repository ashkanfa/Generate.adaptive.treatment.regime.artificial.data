#MIMIC
"
This code estimates joint probability of 't' treatments.
Assumes an additive model form of main effects only
Should be modified for a different model form

The data should be pre-processed before importing
"
rm(list=ls())
setwd("C:/Users/RandKID/OneDrive/Research/Proposal/Generating Simulated Data")
library(MASS)
library(gtools)
library(nnet)
library(geepack)
library(bindata)
library(combinat)
#Read Data
#Trts_Stage1 <- read.csv("Stage1_AllNewTreatments_Binary_229Obs.csv",header=TRUE)
Trts_Stage2 <- read.csv("Stage2_AllNewTreatments_Binary_229Obs.csv",header=TRUE)
#Cov_Stage2 <- read.csv("NewCombinedbinaryconfoundingvar_stage1_229Obs.csv",header=TRUE)
Cov_Stage2 <- read.csv("NewCombinedbinaryconfoundingvar_stage2_229Obs.csv",header=TRUE)

Data_Corr.Var <- Trts_Stage2
#sample size
n=nrow(Data_Corr.Var)
#number of treatments
t=ncol(Data_Corr.Var)

###MIMIC algorithm using greedy search 

  order_min<-NULL
  order_mean<-NULL
  newdata.mimic <- Data_Corr.Var
  
  predprobs.mimic <- NULL
  #data_mimic <- cbind(newdata.mimic,Cov_Stage2)
  for(j in 1:t)
  {
    data_mimic <- cbind(newdata.mimic[j],Cov_Stage2)
    model.form <- as.formula(paste(colnames(data_mimic)[1],"~","."))
    model.logistic <- multinom(formula=model.form, data = data_mimic)
    summary(model.logistic)
    phat.mimic <- predict(model.logistic, newdata = data_mimic, "probs")#obtain predicted probabilities and store it under phat column for each multinomial model
    phatcol.mimic <- cbind(phat.mimic)
    prob.mimic<-NULL
    for(k in 1:n)
    {
      idx <- data_mimic[k,j]
      if(idx==0)
      {
        probtemp <- 1-phat.mimic[k]
        prob.mimic <- rbind(probtemp,prob.mimic,deparse.level = 0)
      }else{
        probtemp <- phat.mimic[k]
        prob.mimic <- rbind(probtemp,prob.mimic,deparse.level = 0)
      } 
    }
    predprobs.mimic <- cbind(predprobs.mimic,prob.mimic,deparse.level = 0)#bind the probabilities in pred probs for IPTW calculations
  }
  predprobs.mimic<-as.data.frame((log(predprobs.mimic)))
  names(predprobs.mimic)<-c(names(newdata.mimic))
  #min.prob<-NULL
  mean.prob<-NULL
  for(j in 1:ncol(predprobs.mimic))
  {
    #min.p <- min(predprobs.mimic[,j])
    #min.prob <- cbind(min.prob,min.p,deparse.level = 0)
    mean.p<-mean(predprobs.mimic[,j])
    mean.prob<-cbind(mean.prob,mean.p,deparse.level = 0)
  }
  #min.prob<-as.data.frame(min.prob)
  mean.prob <- as.data.frame(mean.prob)
  #names(min.prob)<-c(names(data_mimic))
  names(mean.prob)<-c(names(newdata.mimic))
  
  #first.treat_min<-colnames(min.prob)[which(min.prob==min(min.prob))]
  #order_min<-cbind(order_min,first.treat_min)
  
  first.treat_minmean<-colnames(mean.prob)[which(mean.prob==min(mean.prob))]
  order_mean<-cbind(order_mean,first.treat_minmean)
  
  newdata.mimic<-Data_Corr.Var[ , -which(names(Data_Corr.Var) %in% c(first.treat_minmean))]
  data_itr_mimic <- as.data.frame(Data_Corr.Var[,c(first.treat_minmean)])
  names(data_itr_mimic) <- c(first.treat_minmean)
  newdata.mimic <- cbind(data_itr_mimic,newdata.mimic)


for(i in 2:(t-1))
{
  predprobs.mimic <- NULL
  #data_mimic <- cbind(newdata.mimic,Cov_Stage2)
  for(j in i:t)
  {
    data_mimic <- cbind(newdata.mimic[j],newdata.mimic[1],Cov_Stage2)
    model.form <- as.formula(paste(colnames(data_mimic)[1],"~","."))
    model.logistic <- multinom(formula=model.form, data = data_mimic)
    summary(model.logistic)
    phat.mimic <- predict(model.logistic, newdata = data_mimic, "probs")#obtain predicted probabilities and store it under phat column for each multinomial model
    phatcol.mimic <- cbind(phat.mimic)
    prob.mimic<-NULL
    for(k in 1:n)
    {
      idx <- data_mimic[k,j]
      if(idx==0)
      {
        probtemp <- 1-phat.mimic[k]
        prob.mimic <- rbind(probtemp,prob.mimic,deparse.level = 0)
      }else{
        probtemp <- phat.mimic[k]
        prob.mimic <- rbind(probtemp,prob.mimic,deparse.level = 0)
      } 
    }
    predprobs.mimic <- cbind(predprobs.mimic,prob.mimic,deparse.level = 0)#bind the probabilities in pred probs for IPTW calculations
  }
  predprobs.mimic<-as.data.frame((log(predprobs.mimic)))
  names(predprobs.mimic)<-c(names(newdata.mimic[,i:t]))
  #min.prob<-NULL
  mean.prob<-NULL
  for(j in 1:ncol(predprobs.mimic))
  {
    #min.p <- min(predprobs.mimic[,j])
    #min.prob <- cbind(min.prob,min.p,deparse.level = 0)
    mean.p<-mean(predprobs.mimic[,j])
    mean.prob<-cbind(mean.prob,mean.p,deparse.level = 0)
  }
  #min.prob<-as.data.frame(min.prob)
  mean.prob <- as.data.frame(mean.prob)
  #names(min.prob)<-c(names(data_mimic))
  names(mean.prob)<-c(names(newdata.mimic[,i:t]))
  
  #first.treat_min<-colnames(min.prob)[which(min.prob==min(min.prob))]
  #order_min<-cbind(order_min,first.treat_min)
  
  first.treat_minmean<-colnames(mean.prob)[which(mean.prob==min(mean.prob))]
  order_mean<-cbind(order_mean,first.treat_minmean)
  
  newdata.mimic<-newdata.mimic[ , -which(names(newdata.mimic) %in% c(first.treat_minmean))]
  data_itr_mimic <- as.data.frame(Data_Corr.Var[,c(first.treat_minmean)])
  names(data_itr_mimic) <- c(first.treat_minmean)
  newdata.mimic <- cbind(data_itr_mimic,newdata.mimic)
}
final.trt <- colnames(Data_Corr.Var)[-which(names(Data_Corr.Var) %in% c(order_mean))]
order_mean<-cbind(order_mean,final.trt)



