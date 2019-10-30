library("devtools")
library("githubinstall")
library("mvtnorm")
library("BinNor")
library("PReMiuM")
library("corpcor")
library("psych")
library("Matrix")
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
#library("boot")

###################################################################################################
## Generating upper diagnol correlation vector by tru.p
## three correlation structure of 1-within important (cor1), 2- btw important and unimportant (cor2)
## and 3- within unimportant (cor3)

"Levels: a) Low [0,0.2], b) High [0.7,0.9], c) Medium [0.3,0.6] ."  

cor.vec.generate=function(p, true.p, cor1="med",cor2="med",cor3="med"){
  #n: number of variables
  # p: fraction/proportion of important variables in n.
  #define necessary functions
  low=function(no.low){
    round(runif(no.low,min=0,max=0.2),digits = 1)
  }
  
  med=function(no.med){
    round(runif(no.med,min=0.3,max=0.5),digits = 1)
  }
  
  high=function(no.high){
    round(runif(no.high,min=0.7,max=0.9),digits = 1)
  }
  
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
#cor.vec<-cor.vec.generate(n=p,true.p,cor1="high",cor2="low")

###################################################################################################
#      FUNCTION : UPPER TRI VEC: 
#      turns matrix into vector of upper-triangular elements

upper_tri_vec = function(m) {
  v1 = as.vector( t(m) )
  keepElement = as.vector( t(upper.tri(m) ) ) #use transpose to avoid going by columns
  v2 = as.numeric( v1[keepElement] )
  return(v2)
}
###################################################################################################


sum.row<-function(x){ 
  v<-matrix(nrow = dim(x)[1],ncol=1)
  for (i in 1:dim(x)[1])
    v[i,]=sum(x[i,])
  return(v)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}