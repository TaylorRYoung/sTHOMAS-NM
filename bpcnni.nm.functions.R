################################################################################
#### Functions for training, predicting, and getting z-scores by iterating 
#### over a vector of regions, the formula is specified in the function but
#### could also be passed as a variable to the function
library(data.table)
library(mfp)

## Train a model using control data only and save it to the specified directory
mfp_train_save_bpcnni <- function(region,dat,out_dir,cv=c("age","male","female","eTIV","site")) {
  
  print(region)

  ## Subset the dataset
  dat  <- dat[phenotype=="Control", c(region,cv), with = FALSE]
  
  ## Train the model
  mfp.fit <- mfp(as.formula(
    sprintf("%s ~ fp(age*female,df=4) + fp(age*male,df=4) + fp(eTIV*female,df=4) + fp(eTIV*male,df=4) + site",region)),
    data = dat,
    maxits = 10000,
    x=FALSE,
    y=FALSE,
    rescale = FALSE,
    verbose=FALSE)
  
  ## Remove ancillary data stored in the model fit object
  mfp.fit$fit$data <- NULL
  mfp.fit$fit$model <- NULL
  mfp.fit$y <- NULL
  mfp.fit$x <- NULL
  mfp.fit$weights <- NULL
  mfp.fit$prior.weights <- NULL
  mfp.fit$fit$weights <- NULL
  mfp.fit$fit$prior.weights <- NULL
  mfp.fit$fitted.values <- NULL
  mfp.fit$effects <- NULL
  mfp.fit$qr <- NULL
  mfp.fit$fit$residuals <- NULL
  mfp.fit$fit$fitted.values <- NULL
  mfp.fit$fit$effects <- NULL
  mfp.fit$linear.predictors <- NULL
  mfp.fit$fit$linear.predictors <- NULL

  ## The items in the object need to be restored when loading the model, see mfp_predict_bpcnni
  attr(mfp.fit$terms, ".Environment") <- NULL
  attr(mfp.fit$fit$terms, ".Environment") <- NULL
  
  ## Save the model object
  saveRDS(mfp.fit,paste0(out_dir,region,"_mfp_bpcnni.rds"), compress = "xz")
  
  return(NULL)
}

## Read in the saved model, make predictions, and calculate z-scores
mfp_predict_bpcnni <- function(region,dt,mod_dir,cv=c("age","male","female","eTIV","site")) {

  ## Read in the model object
  fit <- readRDS(paste0(mod_dir,region,"_mfp_bpcnni.rds"))

  ## These items need to be added back to the model fit object
  attr(fit$terms, ".Environment") <- globalenv()
  attr(fit$fit$terms, ".Environment") <- globalenv()
  
  ## Calculate the standard deviation of the residuals
  sg  <- sd(fit$residuals)
  
  ## Subset the dataset to only variables in the model and remove NAs
  dt <- dt[, c(region,cv), with = FALSE]
  #dt <- dt[!is.na(unlist(dt[, 1]))]
  
  ## Predict on the new data, returns a vector
  ft <- predict(fit,dt)
  
  ## Calculate z-scores, subsets the dataset to the first variable, which is the 
  ## observed value in the dataset so you have (observed-expected)/standard deviation
  z  <- (dt[,1][[1]]-ft)/sg
  
  ## Returns a vector only
  return(z)
}

## Wrapper function to get z-scores using mfp_predict_bpcnni and merge with original data
mfp_get_z_bpcnni <- function(region_lst,dt,mod_dir) {
  
  z <- sapply(region_lst, mfp_predict_bpcnni, dt=dt, mod_dir=mod_dir, simplify="matrix", USE.NAMES=TRUE)
  rownames(z) <- dt$id
  z <- as.data.table(z, keep.rownames = "id")
  z <- merge(dt[, names(dt)[!names(dt) %in% region_lst], with=FALSE],z)
  
  return(z)
}

## Function to shift z-scores of new data, this is equivalent to retraining the model
## Site in the new data needs to match and existing site in the trained models
mfp_shift_z_bpcnni <- function(roi_lst,dt,mod_dir) {
  
  ## Predict and get z
  dz <- mfp_get_z_bpcnni(roi_lst,dt,mod_dir)
  
  ## Shift - calculate mean and sd, then z'=(z-mean)/sd
  dz <- melt(dz, id.vars = names(dz)[!names(dz) %in% roi_lst],
              value.name = "z", variable.name = "roi")
  dm <- dz[phenotype=="Control", .(mean=mean(z), sd=sd(z)), by=roi]
  dz <- merge(dz,dm)[, z := (z-mean)/sd]
  dz <- dcast(dz[, -c("mean","sd"), with=FALSE], ... ~ roi, value.var = "z")
  
  return(dz)
}

## GAMLSS - This still requires the original data , only include variables in the model (e.g. leave out id and phenotype)
## Function to train and save GAMLSS models, needs complete records so na.omit() is used
## The function can get stuck in local minima
gamlss_train_save_bpcnni <- function(region,dt,out_dir) {
  
  print(region)
  
  dt  <- dt[phenotype=="Control", c(region,cv), with = FALSE]
  dt  <- na.omit(dt)
  
  gml <- gamlss(formula = as.formula(
    paste0(sprintf("%s ~ cs(age*female) + cs(age*male) + cs(eTIV*female) + cs(eTIV*male) + random(as.factor(site))",region))
  ),
  sigma.formula = as.formula("~cs(age)"),
  tau.formula = ~1,
  nu.formula = ~1,
  family = "SHASHo",
  control = gamlss.control(n.cyc = 10000, trace = FALSE),
  data = dt)
  
  saveRDS(gml,paste0(out_dir,region,"_gamlss_bpcnni.rds"))
  
}

gamlss_train_save_bpcnni_v2 <- function(region,dt,out_dir,n_iter=10000,num_repeat=10) {
  
  print(region)
  
  dt  <- dt[phenotype=="Control", c(region,cv), with = FALSE]
  dt  <- na.omit(dt)
  
  ## Use a repeat loop to re-fit the model if it reaches the maximum number of iterations
  i <- 0
  
  repeat{
    
    gml <- gamlss(formula = as.formula(
      paste0(sprintf("%s ~ cs(age*female) + cs(age*male) + cs(eTIV*female) + cs(eTIV*male) + random(as.factor(site))",region))
    ),
    sigma.formula = as.formula("~cs(age)"),
    tau.formula = ~1,
    nu.formula = ~1,
    family = "SHASHo",
    control = gamlss.control(n.cyc = n_iter, trace = FALSE),
    data = dt)
    
    i <- i+1
    
    if (gml$iter < n_iter) {
      
      print(paste(region,"converged in",gml$iter,"iterations."))
      break
      
    } else if (i==num_repeat) {
      print(paste(region,"did not converge. Consider increasing the number of iterations."))
      break
    }
  }
  
  saveRDS(gml,paste0(out_dir,region,"_gamlss_bpcnni.rds"))
  
}

## Predict
gamlss_predict_bpcnni <- function(region,dt_new,dt_prd,mod_dir) {
  
  fit <- readRDS(paste0(mod_dir,region,"_gamlss_bpcnni.rds"))
  
  dt_new <- dt_new[, c(region,cv,"id"), with = FALSE]
  dt_new <- na.omit(dt_new)
  dt_id  <- dt_new$id
  dt_new <- na.omit(dt_new[, c(region,cv), with=FALSE])
  dt_prd <- dt_prd[phenotype=="Control", c(region,cv), with = FALSE]
  dt_prd <- na.omit(dt_prd)
  
  prd <- as.data.table(predictAll(fit, newdata = dt_new, output = "data.frame", data=dt_prd))
  prd <- prd[, z := (y-mu)/sigma]
  prd <- prd[, id := dt_id]
  #prd <- prd[, roi := region]
  prd <- prd[, c("id","z"), with=FALSE]
  names(prd)[2] <- region
  
  return(prd)
}

## Predict and get z-scores
gamlss_get_z_bpcnni <- function(region_lst,dt_new,dt_prd,mod_dir) {
  
  # z <- sapply(region_lst, gamlss_predict_bpcnni, dt_new=dt_new, dt_prd=dt_prd, mod_dir=mod_dir, simplify="matrix", USE.NAMES=TRUE)
  # rownames(z) <- dt_new$id
  # z <- as.data.table(z, keep.rownames = "id")
  z <- lapply(region_lst, gamlss_predict_bpcnni, dt_new=dt_new, dt_prd=dt_prd, mod_dir=mod_dir)
  # z <- rbindlist(z)
  # z <- dcast(z,id~roi, value.var = "z")
  # z <- merge(dt[, names(dt)[!names(dt) %in% region_lst], with=FALSE],z, all=TRUE)
  
  return(z)
}
