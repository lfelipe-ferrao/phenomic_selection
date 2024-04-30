
#Name: Paul Adunola
#Date: Apr 29, 2024
#Title: A COMPARISON OF GENOMIC AND PHENOMIC SELECTION METHODS FOR YIELD PREDICTION IN Coffea canephora

library(dplyr)
library(AGHmatrix)
library(stringr)
library(asreml)

#Select data in FES2021 location
fes21 = agro1 %>% filter(envi %in% "FES2021") %>%data.frame()
fes21$RG_Sample_Code = as.factor(fes21$RG_Sample_Code)

# Relationship Matrices
# Genomic Kernel Matrix
Gmat = Gmatrix(as.matrix(X1),ploidy=2)
Gmat[1:4,1:4]
diag(Gmat) = diag(Gmat) + 0.0001

# NIR Kernel Matrix
source('env_kernel.R') #This function is taken from envRtype package for creating relationship matrices
Gnir = env_kernel(as.matrix(mean.cof.nir3),is.scaled=T,gaussian2=T)
Gnir = Gnir$envCov
image(Gnir)
Gnir[1:4,1:4]

#Prediction
set.seed(123)
gs_sing = function(dat,G,rep=5){
  rep_list = list()
  cv = 10
  n = nrow(dat)
  for (j in 1:rep) {
    folds = sample(1:cv, size=n,replace=T)  # Creating the folds 
    corr = vector()
    for(i in 1:max(folds)){              # loop across the folds
      pred.gen=which(folds==i)
      train.set = dat
      # Store the information of individuals that we will predict
      test.set = dat[pred.gen,]
      # Assuming NA in the data for the individuals that we will predict
      train.set$yNA = train.set$mean_yield
      train.set$yNA[pred.gen] = NA
      
      # GBLUP model 
      mod <- asreml(
        fixed = yNA ~ 1, 
        random = ~ vm(RG_Sample_Code, G), 
        na.action=na.method(y="include"),
        data = train.set
      )
      mod = update(mod)
      #Extract blups
      yHat = summary(mod, coef=T)$coef.random
      yHat = data.frame(id = rownames(yHat),yHat)
      rownames(yHat) = NULL
      #split RG_Sample_Code name
      yHat[c("o1","o2","o3","RG_Sample_Code")] = str_split_fixed(yHat$id, '_', 4)
      # Creating the id for later combining
      yHat.set = yHat %>%
        filter(RG_Sample_Code %in% test.set$RG_Sample_Code) %>%
        data.frame()
      # Combining the datasets
      comb_data = merge(yHat.set[,c(8,2)], test.set[,c(10,9)], by="RG_Sample_Code")
      # Prediction accuracy (correlation)
      corr[i] = cor(comb_data$solution,comb_data$mean_yield, use="complete")
    }
    rep_list[[j]] = corr
  }
  return(unlist(rep_list))
}

gs_multi = function(dat,G,K,rep=5){

  rep_list = list()
  cv = 10
  n = nrow(dat)
  for (j in 1:rep) {                        # loop across the reps
    folds = sample(1:cv, size=n,replace=T)  # Creating the folds 
    
    corr = vector()
    for(i in 1:max(folds)){                 # loop across the folds
      pred.gen=which(folds==i)
      train.set = dat
      # Store the information of individuals that we will predict
      test.set = dat[pred.gen,]
      # Assuming NA in the data for the individuals that we will predict
      train.set$yNA = train.set$mean_yield
      train.set$yNA[pred.gen] = NA
      
      # GBLUP model 
      mod <- asreml(
        fixed = yNA ~ 1, 
        random = ~ vm(RG_Sample_Code, G) + vm(RG_Sample_Code, K), 
        na.action=na.method(y="include"),
        data = train.set
      )
      mod = update(mod)
      mod = update(mod)
      #Extract blups
      yHat = summary(mod, coef=T)$coef.random
      yHat = data.frame(id = rownames(yHat),yHat)
      rownames(yHat) = NULL
      #split RG_Sample_Code name
      yHat[c("o1","o2","o3","RG_Sample_Code")] = str_split_fixed(yHat$id, '_', 4)
      # Creating the id for later combining
      yHat.set = yHat %>%
        filter(RG_Sample_Code %in% test.set$RG_Sample_Code) %>%
        group_by(RG_Sample_Code) %>%
        summarise(solution = mean(solution)) %>%
        data.frame()
      # Combining the datasets
      comb_data = merge(yHat.set, test.set[,c(10,9)], by="RG_Sample_Code")
      # Prediction accuracy (correlation)
      corr[i] = cor(comb_data$solution,comb_data$mean_yield, use="complete")
    }
    rep_list[[j]] = corr
  }
  return(unlist(rep_list))
}

#Genomic prediction
gs_fem21 = gs_sing(fes21,Gmat)
mean(gs_fem21)

#Phenomic prediction
nir_fem21 = gs_sing(fes21,Gnir)
mean(nir_fem21)

#Genomic + Phenomic prediction
multi_fem21 = gs_multi(fes21,Gmat,Gnir)
mean(multi_fem21)

