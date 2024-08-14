library(multilevLCA)

multiLCA_adj3step = function(data, Y, iT, id_high, iM, Z, Zh = NULL,
                             extout = FALSE, dataout = TRUE, kmea = TRUE,
                             maxIter = 1e3, tol = 1e-8, reord = TRUE,
                             NRmaxit = 100, NRtol = 1e-6, verbose = TRUE){
  
  # Step 1: single-level measurement model
  step1 = multiLCA(data,Y,iT,extout=TRUE)
  
  # Step 2: class assignment
  data = as.data.frame(merge(data,unique(step1$mU_modal),Y))
  
  # Step 3: structural model
  step3_init = multiLCA(data,Y,iT,id_high,iM,extout=TRUE)
  
  Y = paste0("C",1:iT)
  if(verbose){
    pb=txtProgressBar(char="Cleaning data...",width=1)
    setTxtProgressBar(pb,1)
    close(pb)
  }
  multilevLCA:::check_inputs1(data,Y,iT,id_high,iM,Z,Zh)
  nrow_data = nrow(data)
  data      = data[complete.cases(data[,Y]),c(Y,id_high,Z,Zh)]
  approach  = multilevLCA:::check_inputs2(data,Y,iT,id_high,iM,Z,Zh)
  if(approach != "direct"){
    stop("This function does not perform model selection",call.=FALSE)
  }
  reord = as.numeric(reord)
  if(!is.null(id_high)){
    id_high_levs    = sort(unique(data[,id_high]))
    data[,id_high]  = as.numeric(factor(data[,id_high]))
    data            = data[order(data[,id_high]),]
    id_high_name    = id_high
    id_high         = data[,id_high]
    vNj             = table(id_high)
  }
  mY        = as.matrix(data[,Y])
  ivItemcat = apply(mY,2,function(x){length(unique(x))})
  itemchar  = any(apply(data[,Y],2,function(x){!is.numeric(x)}))
  if(any(ivItemcat>2)&!itemchar){
    mY  = multilevLCA:::update_YmY(mY,Y,ivItemcat)
    Y   = mY$Y
    mY  = mY$mY
  } else if(itemchar){
    mY  = multilevLCA:::update_YmY_nonnum(mY,Y,ivItemcat)
    Y   = mY$Y
    mY  = mY$mY
  }
  if(!is.null(Z)){
    mZ  = multilevLCA:::clean_cov(data,Z)
    Z   = mZ$covnam
    mZ  = mZ$cov
  }
  if(!is.null(Zh)){
    mZh = multilevLCA:::clean_cov(data,Zh)
    Zh  = mZh$covnam
    mZh = mZh$cov
  }
  
  if(is.null(Zh)){
    
    P             = ncol(mZ)
    vOmega_start  = step3_init$vOmega
    cGamma_start  = array(c(rbind(step3_init$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
    mPhi_start    = t(step1$mClassErrProb)
    mStep1Var     = step3_init$Varmat
    #
    mY      = mY[complete.cases(mZ),]
    id_high = id_high[complete.cases(mZ)]
    vNj     = table(id_high)
    mZ      = mZ[complete.cases(mZ),]
    #
    if(verbose){
      pb=txtProgressBar(char="Fitting structural model...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
    
    step3 = multilevLCA:::MLTLCA_cov_poly(mY,
                                          mZ,
                                          vNj,
                                          vOmega_start,
                                          cGamma_start,
                                          mPhi_start,
                                          mStep1Var,
                                          ivItemcat=iT,
                                          maxIter,
                                          tol,
                                          fixedpars=1,
                                          nsteps=3,
                                          NRtol,
                                          NRmaxit)
    step3 = multilevLCA:::clean_output4(step3,
                                        Y,
                                        iT,
                                        iM,
                                        c("Intercept",Z),
                                        mY,
                                        mZ,
                                        id_high,
                                        length(Y),
                                        P,
                                        id_high_levs,
                                        id_high_name,
                                        extout,
                                        dataout,
                                        ivItemcat)
    
  } else{
    
    P                 = ncol(mZ)
    P_high            = ncol(mZh)
    cGamma_start      = array(c(rbind(step3_init$mGamma,matrix(0,(iT-1)*(P-1),iM))),c(iT-1,P,iM))
    mPhi_start        = t(step1$mClassErrProb)
    mStep1Var         = step3_init$Varmat
    mDelta_start      = matrix(0,iM-1,P_high)
    mDelta_start[,1]  = step3_init$vAlpha
    #
    nomissing = complete.cases(mZ)&complete.cases(mZh)
    mY        = mY[nomissing,]
    id_high   = id_high[nomissing]
    vNj       = table(id_high)
    mZ        = mZ[nomissing,]
    mZh       = mZh[nomissing,]
    mZh       = mZh[!duplicated(id_high),]
    #
    if(verbose){
      pb=txtProgressBar(char="Fitting structural model...",width=1)
      setTxtProgressBar(pb,1)
      close(pb)
    }
    
    step3 = multilevLCA:::MLTLCA_covlowhigh_poly(mY,
                                                 mZ,
                                                 mZh,
                                                 vNj,
                                                 mDelta_start,
                                                 cGamma_start,
                                                 mPhi_start,
                                                 mStep1Var,
                                                 ivItemcat=iT,
                                                 maxIter,
                                                 tol,
                                                 fixedpars=1,
                                                 nsteps=3,
                                                 NRtol,
                                                 NRmaxit)
    step3 = multilevLCA:::clean_output5(step3,
                                        Y,
                                        iT,
                                        iM,
                                        c("Intercept",Z),
                                        c("Intercept",Zh),
                                        mY,
                                        mZ,
                                        mZh,
                                        id_high,
                                        length(Y),
                                        P,
                                        P_high,
                                        id_high_levs,
                                        id_high_name,
                                        extout,
                                        dataout,
                                        ivItemcat)
    
  }
  
  if(nrow(data)<nrow_data){
    warning(paste("Missing values in columns for indicators,",
                  "sample size for estimation of measurement model:",
                  nrow(data)),call.=FALSE)
  }
  if(!is.null(Z)){
    if(nrow(na.omit(mZ))<nrow(data)){
      warning(paste("Missing values in columns for covariates,",
                    "sample size for estimation of structural model:",
                    nrow(na.omit(mZ))),call.=FALSE)
    }
  }
  step3$call = match.call()
  class(step3) = "multiLCA"
  
  return(list(step1=step1,step3=step3))
  
}

# Example 1: multilevel model with lower-level covariates
ex1 = multiLCA_adj3step(data    = dataTOY,
                        Y       = colnames(dataTOY)[2:11],
                        iT      = 3,
                        id_high = "id_high",
                        iM      = 2,
                        Z       = "Z_low")
ex1$step1 # Step-1 model
ex1$step3 # Step-3 model (the model of interest)

# Example 2: multilevel model with lower-level and higher-level covariates
ex2 = multiLCA_adj3step(data    = dataTOY,
                        Y       = colnames(dataTOY)[2:11],
                        iT      = 3,
                        id_high = "id_high",
                        iM      = 2,
                        Z       = "Z_low",
                        Zh      = "Z_high")
ex2$step1 # Step-1 model
ex2$step3 # Step-3 model (the model of interest)
