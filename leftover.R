#### Step 2 Creatin "Within Subject structure"
### Step 2.1: Generate one slice of data at the begining and 
#Prepration before defining a function to expand the variables longitudinally in step 2.2
multi_stage_data = function(no.rows, n.stage, no.bin, no.nor, no.trt, prop.vec.bin, mean.vec.nor, var.nor, pcor,
                            no.bin.tv=0, no.nor.tv=0, no.nor.tvt=0, wcor) {
    
    # n: number of subjects to generate    
    # n.stage:         number of stages or repeated observations per subject to generate
    # pcor:            across-subject correlation matrix for time slot 0
    # time.invariant   if false,then there are time variant covariates for which the number of them
    #in each category (normal or binary) is determined by no.bin.tv and no.nor.tv 
    # no.bin.tv        number of time-variant binary covariates that is not function of time
    # no.nor.tv        number of time-variant normal covariates that is not function of time
    # no.nor.tvt       number of time-variant normal covariates that is a function of time
    #                  so far there is no binary variable that is a function of time. only for normal
    no.bin.ti = no.bin-no.bin.tv 
    #no.bin.ti         number of time invariant binary covariates
    no.nor.ti = no.nor - no.nor.tv - nor.tvt
    #no.nor.ti         number of time invariant normal covariates
    # wcor:           within-subject correlation matrix
    
    
    #Series of control statements to prevent obvious ARGUMENT SPECIFICATION ERRORS:
    
    if ((n.stage<1)|(floor(n.stage)!=n.stage)){stop("Number of stages must be
                                                    an integer whose value is at least 1!\n")}
    
    if ((no.bin.tv!=0)|(no.bin.tv<0)|(floor(no.bin.tv)!=no.bin.tv)|(no.bin.tv>no.bin)){stop("no.bin.tv is whether 0 when there is no time-variant binary,
                                                                                            or an integer between 1 and total number of binary variables when there is at least one time variant binary variable!\n")}
    
    if ((no.nor.tv!=0)|(no.nor.tv<0)|(floor(no.nor.tv)!=no.nor.tv)|(no.nor.tv>no.nor)){stop("no.nor.tv is whether 0 when there is no time-variant normal var,
                                                                                            or an integer between 1 and total number of normal variables when there is at least one time variant normal variable!\n")}
    
    if ((no.nor.tvt!=0)|(no.nor.tvt<0)|(floor(no.nor.tvt)!=no.nor.tvt)|(no.nor.tvt>no.nor)){stop("no.nor.tvt is whether 0 when there is no time-variant normal var,
                                                                                                 or an integer between 1 and total number of normal variables when there is at least one time variant normal variable that changes as a function of time!\n")}
    
    if((no.nor.ti+no.nor.tv+no.nor.tvt)!=no.nor){stop("total number of normal variables are not sum of no.nor.ti + no.nor.tv + no.nor.tvt")}
    
    if (no.nor.tvt==0){
        if ((no.bin.tv<0)|(floor(no.bin.tv)!=no.bin.tv)|(no.bin.tv>no.bin)){stop("when there is time variant binary, the Number of time-variant binary variables
                                                                                 must be an integer whose value is at least 0 and at most total number of binary variables!\n")}
        if ((no.nor.tv<1)|(floor(no.nor.tv)!=no.nor.tv)|(no.nor.tv>no.nor)){stop("when there is time variant normal, theNumber of time-variant normal variables
                                                                                 must be an integer whose value is at least 1 and at most total number of normal variables !\n")}
        #dimension of time-variant 
        #  d_tv = no.nor.tv + no.bin.tv
        #  if(length(wcor)!=(d*(d-1)/2)){
        # stop("Vector of correlations is misspecified, dimension is wrong!\n")}
        
        #------------------------------------------------------------------------------------------------------------
        pcor=upper_tri_vec(pcor)
        
        
        
        slot_0 = mod.jointly.generate.binary.normal( no.rows, no.bin, no.nor, prop.vec.bin, mean.vec.nor, var.nor, corr.vec=pcor)
        #no.rows=Number of subjects
        #no.bin=Number of binary variables (trt vars and non-trt covariates e.g. sex )
        #no.nor=Number of normally distributed variables
        #prop.vec.bin=Vector of marginal proportions for binary variables
        #mean.vec.nor=Vector of means for normal variables
        #var.nor=Vector of variances for normal variables
        #corr.vec=Specified correlations among all variables 
        colname(slot_0) <- paste("x", 1:d ,sep = "")  
        
        #for now, I devide covariates into 4 categories. 
        # 1- only outcome related covariates,           cov_1
        # 2- only trt related covariates                cov_2
        # 3- both outcome and trt related covariates    cov_3
        # 4- covariates of no effect.                   cov_4
        
        cov_2    <-slot_0 [ , 1:floor(quantile(1:d,.25))]
        cov_1    <-slot_0 [ , floor(quantile(1:d,.25))+1 : floor(quantile(1:d,.5))]
        cov_3    <-slot_0 [ , floor(quantile(1:d,.5)+1) : floor(quantile(1:d,.75))]
        cov_4    <-slot_0 [ , floor(quantile(1:d,.75)+1) : d]
        
        #we'll plan to get half the trts from class2 and second half from class 3
        #First half:
        th1=floor(quantile(1:no.trt),.5)
        # 70 percentile of the total number of covariates of class 2
        min_cov_2=floor(quantile(1:dim(cov_2)[2],0.7))
        # sample from 70 percentile to all number of covariates of class 2 generated in half number os trts
        num_cov_2=sample(min_cov_2 : dim(cov_2)[2], th1, replace = TRUE)
        #now that we have the number of covariates
        for (i in th1){
            for (j in sample(1:dim(cov_2)[2],num_cov_2[i],replace = FALSE)){
                
            }
        }
        
        paste(colnames(cov_2))[]
        
        p1 <- data.frame(exp(0+0.3*x.2+0.85*x.3)/(1+exp(0+0.3*x.2+0.85*x.3)))
        
        
        # FORM HALF THE TRTS, FROM only_trt_COV
        for (i in 1:th1){
            p[i]=data.frame(exp(0+0.3*x.2+0.85*x.3)/(1+exp(0+0.3*x.2+0.85*x.3)))
        }
        only_trt_cov[,]
        # we would like to have trts of this class (be dependent on more than 70 percent OF THE CLASS 2 COVARIATES
        
        
        
        if   
        
        
        
        
        for (i in 1:no.trt){
            p[i]=data.frame(exp(0+0.3*x.2+0.85*x.3)/(1+exp(0+0.3*x.2+0.85*x.3)))
        }
        
        
        
        }
    
    
    
    
    
    
    
    
    