
#C-1 still in progress
#B-4 still in progress
#B-5 still in progress
#B-1 independency should be checked

#------------------------------------------------------------------------------------------------------------------
#This simulated Data is able to produce the following structure so far (soft coded):
#A -   Continuous, Discrete Structure in covariates and trts:

#A-1      Optional mixture of Continuous, categorical and binary covariates.
#A-2      Optional trts (decision variables or exposures or policies).
#The default here is Binary.

#A-3      Response: continuous

#------------------------------------------------------------------------------------------------------------------
#B -  Dependency Structure "across" variables at time t:

#B-1      full/partial Independency or full/partial Multicollinearity across Trts
#B-2      full/partial Independency or full/partial Multicollinearity across covariates
#B-3      full/partial dependency btw covariates and trts.

#B-4      Partial dependency of covariates and response: to discover performance of
#the proposed models in exploring true underlying model and also performance of dimension
#reduction techniques in selecting the true variables in the model.

#B-5      Partial dependency of Trts and response: to explore the performance of modeling techniques
#in keeping the true trts in the model.
#------------------------------------------------------------------------------------------------------------------
#C -  Dependency Structure "within" subject evolving in time:

#C-1      Time-varying covariates
#C-1-1         Time-varying covariates that do not vary as a function of time(e.g. temperature)
#C-1-2         Time-varying covariates that does vary as a function of time
#C-2      Dependency of an trt at stage t on previous trts
#C-3      Optional Time-varying confounding
#(response is an intermediate variable effecting trt in next stage and affected by trt from the last stage )
#to explore best model that can incorporate IPTW weights that has been developed by Nilabh et al.
#------------------------------------------------------------------------------------------------------------------
#D       Simulate High-Dimensional data
#D1      The algorithm is capable of creating high dimensional covariates or treatments. (High dimensionality                      "across" variables at time t)
#D2      The algorithm is capable of creating determined number of stages. (High dimensionality                                             "within" variables at time t)



#STEP 1-------Creating (binary and continuous) covariates & binary trts and the correlation structure across these variables
#USING THE modified version of PACKAGE BinNor as the basis for creating mixture of Binary and continuous trts and covariates
#  which are partially/fully correlated. by running the function "mod.jointly.generate.binary.normal"

#????????????????????Explain how it works based on the paper here!!!!!!!!!!!!!
#??????? also explain that this step creates what A to D in explained before!!!
#????? also give an example how to work with this algorithm.


#------------------------------------------------------------------------------------------------------------------
#Step PRE: Some functions that is used later on in the mod.jointly.generate.binary.normal function algorithm 
#and so on

#         Function1: closest element
#         returns whichever element of candidates vector is closest in abs value to x

closest = function(x, candidates) {
    return( candidates[ which.min( abs(candidates - x) ) ] )
}


#        Function2: binary-normal correlation bound ######
#        given parameters for 1 binary and 1 normal RV, return the maximum correlation

BN.rBound = function(p) {
    q = 1-p
    
    #compute upper bound
    hiBound = dnorm( qnorm(p) ) / sqrt(p*q)
    return( round(hiBound, 2) )
}


#       Function3: expand a matrix of slot_0 to the number of specified stages for time_invariant covariate setting
#       repeat each subject's entry in each row for obs number of times
expand_matrix = function(.matrix, .n.stage) {
    library(plyr)
    cat("Expanding a time invariant covariates to the number of specified stages\n")
    
    .n = nrow(.matrix)
    .expanded = matrix(c(NA), nrow = .n*.n.stage, ncol = ncol(.matrix) )
    
    adply(.matrix, 1, function(..subject, ..obs) {
        matrix(rep(..subject, ..obs), nrow = ..obs, byrow = TRUE)
    }, .n.stage)
    
}

#      FUNCTION4: UPPER TRI VEC: 
#      turns matrix into vector of upper-triangular elements

upper_tri_vec = function(m) {
    v1 = as.vector( t(m) )
    keepElement = as.vector( t(upper.tri(m) ) ) #use transpose to avoid going by columns
    v2 = as.numeric( v1[keepElement] )
    return(v2)
}
#------------------------------------------------------------------------------------------------------------------
#Step 1: Creating the "Across" data

mod.jointly.generate.binary.normal=function(no.rows,no.bin,no.nor,
                                            prop.vec.bin,mean.vec.nor,var.nor,corr.vec, adjust.corrs = TRUE){
    
    require(mvtnorm)
    require(corpcor)
    require(psych)
    require(Matrix)
    require(BinNor)
    require(ICC)
    require(miscTools)
    require(car)
    require(plyr)
    #Definition of the arguments are as follows:
    #no.rows=Number of subjects
    #no.bin=Number of binary variables (among these )
    #no.nor=Number of normally distributed variables
    #prop.vec.bin=Vector of marginal proportions for binary variables
    #mean.vec.nor=Vector of means for normal variables
    #var.nor=Vector of variances for normal variables
    #corr.vec=Specified correlations among all variables
    b=2*no.bin #b is double number of the total binary vars
    d=no.bin+no.nor #d is the total dimension 
    #adjust.corrs=T/F. If a correlation is out of bounds, should it
    # be adjusted to the closest feasible values?
    
    ############################################################
    #Important note 1: For convenience, binary variables are assumed
    #to come first, then normal variables follow
    #Important note 2: Correlations are specified in vector form,
    #rather than a matrix form. If the dimension is d, d*(d-1)/2
    #non-redundant correlation terms must be specified. The order
    #in which correlations are specified is based on the upper diagonal
    #elements. For example, if there are four variables (X1,X2,X3,X4),
    #corr.vec is specified in the following form:
    #c(Corr(X1,X2),Corr(X1,X3),Corr(X1,X4),Corr(X2,X3),Corr(X2,X4),
    #Corr(X3,X4))
    #------------------------------------------------------------------------------------------------------------------
    #Step 1-1 Series of control statements to prevent obvious ARGUMENT
    #SPECIFICATION ERRORS:
    
    if ((no.rows<1)|(floor(no.rows)!=no.rows)){stop("Number of rows must be
                                                    an integer whose value is at least 1!\n")}
    if ((no.bin<1)|(floor(no.bin)!=no.bin)){stop("Number of binary variables
                                                 must be an integer whose value is at least 1!\n")}
    if ((no.nor<1)|(floor(no.nor)!=no.nor)){stop("Number of normal variables
                                                 must be an integer whose value is at least 1!\n")}
    
    if ((min(prop.vec.bin)<=0)|(max(prop.vec.bin)>=1)){
        stop("Proportions for binary variables must be between 0 and 1!\n")}
    if (length(prop.vec.bin)!=no.bin){stop("Proportion vector
                                           is misspecified, dimension is wrong!\n")}
    
    if (length(mean.vec.nor)!=no.nor){
        stop("Mean vector for the normal part is misspecified,
             dimension is wrong!\n")}
    if (length(var.nor)!=no.nor){
        stop("Vector of variances for the normal part is misspecified,
             dimension is wrong!\n")}
    if (min(var.nor<=0)){stop("Variances must be positive!\n")}
    
    if(length(corr.vec)!=(d*(d-1)/2)){
        stop("Vector of correlations is misspecified, dimension is wrong!\n")}
    
    #------------------------------------------------------------------------------------------------------------------
    #Step 1-2:   Statements to check CORRELATION BOUND VIOLATIONS
    
    #Step 1-2-1: Form a correlation matrix from the specified correlation vector
    # this is equivalent to lower.tri.to.corr.mat function from BinNor package
    # that Converts a lower vector to a full correlation matrix
    sigma=diag(d)
    temp=1
    for(i in 1:(d-1)){
        for(j in (i+1):d){
            sigma[i,j]=sigma[j,i]=corr.vec[temp]
            temp=temp+1
        }
    }
    #------------------------------------------------------------------------------------------------------------------
    
    #Step 1-2-1: Check if the specified correlation matrix is positive definite, if not
    #find the nearest positive definite matrix
    #This is equivalent to the function: compute.sigma.star from BinNor package that
    #Computes intermediate (tetrachoric) correlation matrix
    
    if(is.positive.definite(sigma)==FALSE)
    {sigma=make.positive.definite(sigma)
    print("Specified correlation matrix is not positive definite,")
    print("Algorithm will be using the closest positive definite matrix!")}
    
    diag(sigma)=1
    
    p=prop.vec.bin
    q=1-p
    #------------------------------------------------------------------------------------------------------------------
    
    #Step 1-2-2: Checks if the correlation terms are within the feasible range
    #This is equivalent to the function: validation.range from BinNor package that
    #Checks if the correlation terms among binary-binary, binary-normal
    #and normal-normal combinations.are within the feasible range
    
    #Check if the correlations for binary-binary combinations are
    #in the feasible range
    
    #Boundaries for BB =[max(-sqrt((pi*pj)/(qi*qj)),-sqrt((qi*qj)/(pi*pj))),
    #min(sqrt((pi*qj)/(qi*pj)),sqrt((qi*pj)/(pi*qj)))]
    L_BB=diag(no.bin)
    U_BB=diag(no.bin)
    
    for(i in 1:no.bin){
        for(j in 1:no.bin){
            if (i!=j) L_BB[i,j]=L_BB[j,i]=max(-sqrt((p[i]*p[j])/(q[i]*q[j])),
                                              -sqrt((q[i]*q[j])/(p[i]*p[j])))
            if (i!=j) U_BB[i,j]=U_BB[j,i]=min(sqrt((p[i]*q[j])/(q[i]*p[j])),
                                              sqrt((q[i]*p[j])/(p[i]*q[j])))
        }
    }
    
    
    for(i in 1:no.bin){
        for(j in 1:no.bin){
            if(sigma[i,j]<L_BB[i,j] | sigma[i,j]>U_BB[i,j]) {
                if (!adjust.corrs) {
                    # if user does not want to adjust correlations, give error
                    stop("BB corrrelation [", i,",",j,"] is out of range! Specify a feasible number!")
                } else {
                    #adjust correlation to the closest feasible value
                    cat( c("BB corrrelation [", i,",",j,"],", sigma[i,j], ", is out of range! Used closest feasible correlation instead\n"))
                    sigma[i,j] = sigma[j,i] = closest(sigma[i,j], c( L_BB[i,j], U_BB[i,j] ) )
                }
            }
        }
    }
    
    
    #Compute the biserial correlations for binary-normal combinations and
    #check if they are in the feasible range
    
    #temporary matrix
    BN_temp=sigma
    
    # replace the BN values in BN_temp with the corresponding phi values
    for(i in (no.bin+1):d){
        for(j in 1:no.bin){
            BN_temp[i,j]=BN_temp[i,j]/(dnorm(qnorm(p[j]))/sqrt(p[j]*q[j]))
        }
    }
    
    for(i in (no.bin+1):d){
        for(j in 1:no.bin){
            if (BN_temp[i,j]< -1 | BN_temp[i,j]> 1) {
                
                if (!adjust.corrs) {q
                    # if user does not want to adjust correlations, give error
                    stop("BN correlation [", i,",",j,"] is out of range! Specify a feasible number!")
                } else {
                    #adjust correlation to the closest feasible value
                    BN_temp[i,j] = closest(BN_temp[i,j], c(-1, 1))
                    
                }
            }
        }
    }
    
    #keep the BN part of BN_temp matrix
    BN=BN_temp[(no.bin+1):d,1:no.bin]
    
    #Compute the tetrachoric correlations for binary-binary combinations
    
    
    # create sigmaBB matrix by converting BB part of sigma into polychoric correlations
    sigmaBB=diag(no.bin)
    for(i in 1:no.bin){
        for(j in 1:no.bin){
            if (i > j) {
                sigmaBB[i,j] = sigmaBB[j,i] = phi2poly( sigma[i,j] ,p[i],p[j])
                #force symmetry because phi2poly is an imperfect optimization process with rounding error
            }
            #########################################################################
            ###### NOTE: ABOVE ROUNDING OF SIGMA ENTRY IS A LITTLE SKETCH!!!!! ######
            #########################################################################
            
        }
    }
    
    
    #Biserial correlations for binary-normal combinations
    sigmaBN=BN
    
    #Combine all three types (binary-binary, binary-normal, normal-normal)
    #of correlations to construct the overall correlation matrix
    #(Step 7 in the algorithm)
    sigma_new=sigma
    sigma_new[1:no.bin,1:no.bin]=sigmaBB
    sigma_new[(no.bin+1):d,1:no.bin]=sigmaBN
    sigma_new[1:no.bin,(no.bin+1):d]=t(sigmaBN)
    
    
    #Check if the final correlation matrix is positive definite, if not
    #find the nearest positive definite matrix (Step 8 in the algorithm)
    
    if(is.positive.definite(sigma_new)==FALSE) {
        sigma_new=make.positive.definite(sigma_new)
        print("Final correlation matrix is not positive definite,")
        print("Algorithm will be using the closest positive definite matrix!")
    }
    #------------------------------------------------------------------------------------------------------------------
    
    
    #Step 1-3: Generate multivariate normal data (Step 9 in the algorithm)
    data=rmvnorm(no.rows,mean=rep(0,d), sigma=sigma_new)
    
    #------------------------------------------------------------------------------------------------------------------
    #
    ##
    ###
    ####
    ##### shouldn't be abs(data[i,j])<=qnorm(1-p[j])) ???? check BinNOR 
    #Step 1-3-1: Obtain binary variables by the thresholds determined by marginal proportions
    
    for(i in 1:no.rows){
        for(j in 1:no.bin){
            if(data[i,j]<=qnorm(1-p[j])) data[i,j]=0 else data[i,j]=1
        }
    }
    
    #------------------------------------------------------------------------------------------------------------------
    
    #------------------------------------------------------------------------------------------------------------------
    #Step 1-3-3:    Go back to the original scale for normal variables by reverse centering and
    
    for(i in 1:no.rows){
        temp=1
        for(j in (no.bin+1):d){
            data[i,j]=mean.vec.nor[temp]+(data[i,j]*sqrt(var.nor[temp]))
            temp=temp+1
        }
    }
    #------------------------------------------------------------------------------------------------------------------
    
    
    #step 1-4: Output is the data matrix!
    diag(sigma_new)=1
    return(list(data,sigma_new))
    
    }
#------------------------------------------------------------------------------------------------------------------
