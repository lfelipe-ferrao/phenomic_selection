env_kernel <-function(env.data,Y=NULL, is.scaled=TRUE, sd.tol = 1, digits=5,
                      tol=1E-3, merge=FALSE,Z_E = NULL, stages=NULL,
                      env.id='env',gaussian2=FALSE, h.gaussian2=NULL){
  
  if(is.null(stages)){
    K_E <- env_kernel_0(env.data=env.data,Y=Y, is.scaled=is.scaled, sd.tol = sd.tol,
                        tol= tol, merge=merge,Z_E = Z_E,
                        env.id=env.id,gaussian2=gaussian2, h.gaussian2=h.gaussian2)
    
    return(list(varCov=round(K_E[[1]],digits),envCov=round(K_E[[2]],digits)))
  }
  
  if(!is.null(stages)){
    
    K_E <- list()
    K_W <- list()
    for(i in 1:length(stages))
    {
      id <- grep(colnames(env.data),pattern = stages[i])
      K_E[[i]] <- round(env_kernel_0(env.data=env.data[,id],Y=Y, is.scaled=is.scaled, sd.tol = sd.tol,
                                     tol= tol, merge=merge,Z_E = Z_E,
                                     env.id=env.id,gaussian2=gaussian2, h.gaussian2=h.gaussian2)[[2]],digits)
      K_W[[i]] <- round(env_kernel_0(env.data=env.data[,id],Y=Y, is.scaled=is.scaled, sd.tol = sd.tol,
                                     tol= tol, merge=merge,Z_E = Z_E,
                                     env.id=env.id,gaussian2=gaussian2, h.gaussian2=h.gaussian2)[[1]],digits)
      
    }
    names(K_E) = names(K_W) = stages
    
    return(list(varCov=K_W,envCov=K_E))
    
  }
  
}


env_kernel_0 <-function(env.data,Y=NULL, is.scaled=TRUE, sd.tol = 1,
                        tol=1E-3, merge=FALSE,Z_E = NULL,
                        env.id='env',gaussian2=FALSE, h.gaussian2=NULL){
  
  nr<-nrow(env.data)
  nc <-ncol(env.data)
  
  GB_Kernel <-function(X,is.center=FALSE){
    if(isFALSE(is.center)) X = scale(x = X,center = T,scale = F)
    XXl <- X %*% t(X)
    K_G <- XXl/(sum(diag(XXl))/nrow(X)) + diag(1e-6, nrow(XXl))
    return(K_G)
  }
  
  if(!is.matrix(env.data)){stop('env.data must be a matrix')}
  if(isFALSE(is.scaled)){
    Amean <- env.data-apply(env.data,2,mean)+tol
    sdA   <- apply(Amean,2,sd)
    A. <- Amean/sdA
    removed <- names(sdA[sdA < sd.tol])
    env.data <- A.[,!colnames(A.) %in% removed]
    t <- ncol(env.data)
    r <- length(removed)
    cat(paste0('------------------------------------------------','\n'))
    cat(paste0(' Removed envirotype markers:','\n'))
    cat(paste0(r,' from ',t,'\n'))
    cat(paste0(removed,'\n'))
    cat(paste0('------------------------------------------------','\n'))
    
  }
  
  if(isTRUE(merge)){
    if(is.null(Z_E)){
      .DF = data.frame(Y[,env.id])
      names(.DF) ='env'
      Z_E = model.matrix(~0+env,.DF)
    }
    env.data = Z_E %*% env.data
  }
  
  if(isTRUE(gaussian2)){
    O <- gaussian2(x = env.data,h=h.gaussian2)
    H <- gaussian2(x = t(env.data),h=h.gaussian2)
  }
  if(isFALSE(gaussian2)){
    
    O = GB_Kernel(env.data)
    H = GB_Kernel(t(env.data))
    
  }
  
  
  return(list(varCov=H,envCov=O))
}

gaussian2 <- function(x,h=NULL){
  d<-as.matrix(dist(x,upper = TRUE,diag = TRUE))^2
  q <- median(d)
  if(is.null(h)) h <- 1
  
  return(exp(-h*d/q))
}

envK = function(env.data,df.pheno,skip=3,env.id){
  df.pheno <-data.frame(df.pheno)
  env.data <-data.frame(env.data)
  env.data$env <- as.factor(rownames(env.data))
  W <- as.matrix(merge(df.pheno,env.data, by=env.id)[,-c(1:skip)])
  return(W)
  
}
