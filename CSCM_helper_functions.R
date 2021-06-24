# ----------------------------------------------------------- #
# This file contains the helper functions needed to replicate #
# "A (Flexible) Synthetic Control Method for Count Data ... " #
# ----------------------------------------------------------- #

# Note: This file is called automatically by the other R scripts and
#       does not have to be opened manually.

# ------------ FUNCTIONS BELOW ------------- #

# Variable importance functions

.g0_varimp = function(X, y) {
  
  fit = suppressWarnings(glmnet::cv.glmnet(x=as.matrix(scale(X)),
                                           y=as.matrix(y),
                                           nfolds=nrow(X),
                                           family="poisson",
                                           alpha=0))
  
  pars = coef(fit, s="lambda.1se")[-1]
  results = abs(pars)/sum(abs(pars))
  
  return(results)
}

.g0_vmat = function(X,y) {
  return(vmat = diag(x = .g0_varimp(X,y), nrow = ncol(X), ncol = ncol(X)))
}

# Loss functions (including gradient and hessian)

.cscm_loss = function(par,X1,X0,V,lambda,synth.w) {
  W <- par
  loss <- t(X1 - X0 %*% W) %*% V %*% (X1 - X0 %*% W)
  penalty <- lambda*sum((synth.w-W)^2)
  return(as.numeric(loss+penalty))
}

.cscm_gradient = function(par, X1,X0,V,lambda,synth.w) {
  W <- par
  loss <- -2 * t(X0) %*% V %*% (X1-X0%*%W)
  penalty <- 2 * lambda * diag(ncol(X0)) %*% (par - synth.w)
  return(loss+penalty)
}

.cscm_hess = function(par, X1,X0,V,lambda,synth.w) {
  W <- par
  loss <- 2 * t(X0) %*% V %*% X0
  penalty <- 2 * lambda * diag(ncol(X0))
  return(loss+penalty)
}

# Standard synth estimation (fast)

.cscm_init_W = function(X1,X0,V) { #Code from augsynth, credit goes to Ben-Michael et al
  
  #Code from augsynth
  
  X0 = t(X0)
  
  Pmat <- X0 %*% V %*% t(X0)
  qvec <- - t(X1) %*% V %*% t(X0)
  
  n0 <- nrow(X0)
  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))
  
  settings = osqp::osqpSettings(verbose = FALSE,
                                eps_rel = 1e-8,
                                eps_abs = 1e-8)
  sol <- osqp::solve_osqp(P = Pmat, q = qvec,
                          A = A, l = l, u = u, 
                          pars = settings)
  
  return(abs(sol$x))
  
}

# CountSynth minimization

.cscm_solve_W = function(X1,X0,Y0,V=NULL,lambda,synth.w, init=NULL) {
  
  if (is.null(V)==TRUE) {V=.g0_vmat(t(X0),Y0)}
  
  if (is.null(init)==TRUE) {init <- synth.w}
  
  solution.w <- optimx::optimr(par=init,
                               fn=.cscm_loss,
                               gr=.cscm_gradient,
                               hess=.cscm_hess,
                               X1=X1,
                               X0=X0,
                               V=V,
                               lambda=lambda,
                               synth.w=synth.w,
                               method="nlminb",
                               lower=rep(0,ncol(X0)),
                               control=list(kkt=FALSE,
                                            starttests=FALSE))
  return(solution.w)
}


# Helper functions to prepare data for optimization

.cscm_scaleX = function(X1,X0,Z1,Z0) { #Code from Synth, credit goes to Abadie et al
  X0 <- rbind(Z0,X0)
  X1 <- rbind(Z1,X1)
  nvarsV <- dim(X0)[1]
  big.dataframe <- cbind(X0, X1)
  divisor <- sqrt(apply(big.dataframe, 1, var))
  scaled.matrix <- t(t(big.dataframe) %*% (1/(divisor) * diag(rep(dim(big.dataframe)[1], 
                                                                  1))))
  X0.scaled <- scaled.matrix[, c(1:(dim(X0)[2]))]
  if (is.vector(X0.scaled) == TRUE) {
    X0.scaled <- t(as.matrix(X0.scaled))
  }
  X1.scaled <- scaled.matrix[, dim(scaled.matrix)[2]]
  res = list()
  res[["X1"]] <- X1.scaled
  res[["X0"]] <- X0.scaled
  return(res)
}

.cscm_postY0 = function(Y0, t_int,Tmax=NULL) {
  if (is.null(Tmax)==TRUE) {Tmax<-nrow(Y0)}
  y<-as.vector(rowMeans(t(Y0[t_int:Tmax,])))
  return(y)
}

# Find lambda sequence

.cscm_lseq = function(X1,X0,V,synth.w,lambda_min_ratio=1E-16, n_lambda=100) {
  
  n = ncol(X0)
  #Setup n-length coefficient vector
  beta <- synth.w+as.vector(rep(0.0001,n))
  
  target <- .cscm_loss(par=synth.w,X1,X0,V,lambda=0,synth.w=synth.w)
  
  f <- function(x,y,beta) {y-x*sum((synth.w-beta)^2)}
  max.lambda = uniroot(f, y=target, beta=beta, lower=0.0000000000001, upper=10000000000000000000)$root
  
  
  scaler = (lambda_min_ratio) ^ (1/n_lambda)
  lambdasq = max.lambda * (scaler ^ (seq(0:n_lambda)-1))
  
  return(lambdasq)
}

# Fit cscm on lambda sequence

csynth_inner = function(X1, X0, V, synth.w, lambdas) {
  iter = length(lambdas)
  init <- synth.w
  results = list()
  run = 1
  for (j in 1:iter) {
    
    if (isTRUE(run)==1) {
      
      init = synth.w
      
    } 
    
    if (isTRUE(run)>1) {
      init = temp$par}
    
    #temp = stats::optim(par=init,
    #                    fn=.cscm_loss,
    #                    gr=.cscm_gradient,
    #                    X1=X1,
    #                    X0=X0,
    #                    V=V,
    #                    synth.w=synth.w,
    #                    lambda=lambdas[j],
    #                    method="L-BFGS-B",
    #                    lower=rep(0,ncol(X0)))
    
    temp <- optimx::optimr(par=init,
                           fn=.cscm_loss,
                           gr=.cscm_gradient,
                           hess=.cscm_hess,
                           X1=X1,
                           X0=X0,
                           V=V,
                           lambda=lambdas[j],
                           synth.w=synth.w,
                           method="nlminb",
                           lower=rep(0,ncol(X0)),
                           control=list(kkt=FALSE,
                                        starttests=FALSE))
    
    results[[j]] = temp
    run = run+1
  }
  return(results)
}

# Cross-validated CSCM

holdout.csynth = function(Z1,Z0,X1,X0,Y0plot,t_int,Tmax,lambda_min_ratio=1E-16, n_lambda=20, min_1se=FALSE) {
  
  
  # Scale data and prepare data
  sf <- .cscm_scaleX(Z1,Z0,X1,X0)
  y0post <- .cscm_postY0(Y0plot,t_int,Tmax)
  
  # Find lambda sequence for the entire sample
  V.full <- .g0_vmat(t(sf$X0), y0post)
  synth.W.full <- .cscm_init_W(sf$X1,sf$X0,V.full)
  lambda.seq <- .cscm_lseq(sf$X1,
                           sf$X0,
                           V=V.full,
                           synth.w=synth.W.full)
  
  # Create storage
  
  holdoutlist = list()
  
  for (t in 1:nrow(Z0)) {
    
    Z1.holdout = as.matrix(Z1[-t])
    Z0.holdout = Z0[-t,]
    
    sm <- .cscm_scaleX(Z1.holdout,Z0.holdout,X1,X0) #To do: add holdout functionality for X1,X0
    V.h <- .g0_vmat(t(sm$X0), y0post)
    w.init <- .cscm_init_W(sm$X1,sm$X0,V.h)
    
    holdout.models = csynth_inner(X1=sm$X1, 
                                  X0=sm$X0, 
                                  V=V.h,
                                  synth.w=w.init,
                                  lambdas=lambda.seq)
    
    Z1.test = as.matrix(Z1[t])
    Z0.test = Z0[t,]
    
    error_holdout = function(Z1.test, Z0.test, W) {
      
      err <- (Z1.test - Z0.test %*% W)^2
      
      
    }
    
    
    error_csynth = do.call("rbind",lapply(holdout.models, function(p) error_holdout(Z1.test, Z0.test, p$par)))
    
    holdoutlist[[t]] = error_csynth
    
  }
  
  lambda_errors <- rowMeans(do.call("cbind",holdoutlist), na.rm=T)
  lambda_errors_se <- apply(do.call("cbind",holdoutlist),1,function(x)  sd(x) / sqrt(length(x)))
  min_idx = which.min(lambda_errors)
  min_error <- lambda_errors[min_idx]
  min_se <- lambda_errors_se[min_idx]
  lambda_min <- lambda.seq[min_idx]
  lambda_1se <- max(lambda.seq[lambda_errors <= min_error + min_se])
  
  lambda.choose <- ifelse(isTRUE(min_1se)==TRUE, lambda_1se, lambda_min)
  
  
  # Fit full model
  
  final.model = .cscm_solve_W(X1=sf$X1,
                              X0=sf$X0,
                              Y0=y0post,
                              V=V.full,
                              synth.w=synth.W.full,
                              lambda=lambda.choose)
  
  results <- list()
  results[["fit"]] <- final.model
  results[["synth.orig"]] <- synth.W.full
  
  return(results)
  
}

# Main wrapper

countSynth = function(data,
                      predictors=NULL,
                      dependent,
                      unit.variable,
                      time.variable,
                      treatment.identifier,
                      controls.identifier,
                      unit.names.variable=NULL,
                      t_int,
                      K=2,
                      n.lambda=20,
                      lambda_min_ratio=1E-16,
                      min_1se=FALSE,
                      ci.level=0.95,
                      full.model=TRUE,
                      ...) {
  
  # To do: add preliminary checks
  # - Are outcomes non-negative?
  # - Is the data balanced?
  # - Are there missing values?
  # - Etc
  
  
  # No predictors? Supply mean of dependent as predictor for compatibility with Synth::dataprep.
  if (is.null(predictors)==TRUE) {predictors <- dependent}
  
  
  # Helper function to keep Synth quiet
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  
  
  # Preliminaries
  
  t_unique = unique(data[,time.variable])
  t_length = length(t_unique)
  times <- 1:t_length
  min_t = min(times)
  max_t = max(times)
  time.plot = c(min_t:max_t)
  pre.period = min_t:(t_int-1)
  df = data
  df[,"timeindex"] <- data[,time.variable]-(min(t_unique)-1)
  
  t_unique <- unique(df[,"timeindex"])
  
  
  # Define length of the pre (T0) and post (T1) periods
  T0 = length(min_t:(t_int-1))
  T1 = length(T0+1:max_t)
  
  # Length of the holdout samples (either T0/K or T1 depending on which is smaller)
  r = floor(min(c(T0/K,T1)))
  
  
  # Define the minimum and maximum points in each holdout sample
  Holdout = list()
  HL = list()
  for (k in 1:K) {
    HL[[k]] = c((k-1)*r+1, k*r)
    Holdout[[k]] = c(t_unique[min(HL[[k]])], t_unique[max(HL[[k]])])
  }
  
  
  # Dataprep objects
  datapreplist = list()
  message("Preparing data for cross-fitted synthetic control estimation...")
  for (k in 1:K) {
    
    datapreplist[[k]] = quiet(Synth::dataprep(foo=df,
                                              predictors=predictors,
                                              dependent=dependent,
                                              unit.variable=unit.variable,
                                              time.variable="timeindex",
                                              treatment.identifier=treatment.identifier,
                                              controls.identifier=controls.identifier,
                                              time.predictors.prior = pre.period[!(pre.period %in% min(Holdout[[k]]):max(Holdout[[k]]))],
                                              time.optimize.ssr = pre.period[!(pre.period %in% min(Holdout[[k]]):max(Holdout[[k]]))],
                                              #unit.names.variable=unit.names.variable,
                                              time.plot=times))
    
  }
  
  augatts = augests = cscmlist = cscm.weights = list()
  message("Estimating unit weights in each training sample...")
  for (k in 1:K) {
    
    
    fit <- holdout.csynth(Z1=datapreplist[[k]]$Z1,
                          Z0=datapreplist[[k]]$Z0,
                          X1=datapreplist[[k]]$X1,
                          X0=datapreplist[[k]]$X0,
                          Y0plot=datapreplist[[k]]$Y0plot,
                          lambda_min_ratio = lambda_min_ratio,
                          n_lambda=n_lambda,
                          t_int=t_int,
                          Tmax=max_t,
                          min_1se=min_1se)
    
    
    #Store objects and info
    cscmlist[[k]] = fit[["fit"]]
    
    #Predict at all time points
    
    opt.w <- fit[["fit"]]$par
    cf.out = datapreplist[[k]]$Y0plot %*% opt.w
    obs.out = datapreplist[[k]]$Y1plot
    
    att.raw = log(sum(obs.out[(T0+1):t_length]))-log(sum(cf.out[(T0+1):t_length]))
    bias.est = log(sum(obs.out[min(HL[[k]]):max(HL[[k]])]))-log(sum(cf.out[min(HL[[k]]):max(HL[[k]])]))
    att = att.raw-bias.est
    
    augatts[[k]] = att
    augests[[k]] = data.frame(obs.out=obs.out, cf.out=cf.out, diff = obs.out-cf.out, 
                              time=min_t:max_t)
    cscm.weights[[k]] = opt.w
  }
  
  #Average the debiased holdout estimates
  att.mean = mean(do.call("rbind", augatts))
  
  # Estimate the standard deviation
  att.vec = as.vector(do.call("rbind", augatts))
  squared.err = (att.vec-att.mean)^2
  att.sd = sqrt(1+(K*r)/T1) * sqrt( (1/(K-1)) * sum(squared.err) ) 
  
  # Confidence intervals based on t-distribution with K-1 degrees of freedom
  
  alf = (1-ci.level)/2
  
  att.lower = att.mean-qt(1-alf, K-1)*(att.sd/sqrt(K))
  att.upper = att.mean+qt(1-alf, K-1)*(att.sd/sqrt(K))
  
  # Store RR estimates
  
  RR = exp(att.mean)
  RR.lower = exp(att.lower)
  RR.upper = exp(att.upper)
  
  
  # Run Csynth on full sample
  
  if (isTRUE(full.model)==TRUE) {
  
  message("Estimating unit weights on the full sample...")
  
  dataprep.full <- Synth::dataprep(foo=df,
                                   predictors=predictors,
                                   dependent=dependent,
                                   unit.variable=unit.variable,
                                   time.variable="timeindex",
                                   treatment.identifier=treatment.identifier,
                                   controls.identifier=controls.identifier,
                                   time.predictors.prior = pre.period,
                                   time.optimize.ssr = pre.period,
                                   #unit.names.variable=unit.names.variable,
                                   time.plot=times)
  
  
  csynth.full <- holdout.csynth(Z1=dataprep.full$Z1,
                                Z0=dataprep.full$Z0,
                                X1=dataprep.full$X1,
                                X0=dataprep.full$X0,
                                Y0plot=dataprep.full$Y0plot,
                                lambda_min_ratio = lambda_min_ratio,
                                n_lambda=n_lambda,
                                t_int=t_int,
                                Tmax=max_t,
                                min_1se=min_1se)
  
  
  solution.w <- csynth.full[["fit"]]$par
  synth.w.main <- csynth.full[["synth.orig"]]
  
  }
  
  message("Estimating effects and finishing up... Done.")
  
  # Store estimates
  require(dplyr)
  estdata = invisible(do.call("rbind", augests) %>% group_by(time) %>% summarise_all(mean))
  colnames(estdata) = c("time",
                        "observed",
                        "est.counterfactual.cross.fitted",
                        "pointwise.effect.cross.fitted")
  
  if (isTRUE(full.model)==TRUE) {
  estdata$est.counterfactual <- dataprep.full$Y0plot %*% solution.w
  estdata$pointwise.effect <- dataprep.full$Y1plot - estdata$est.counterfactual
  }
  
  add <- function(x) Reduce("+", x)
  cross.fitted.w = add(cscm.weights)/K
  
  #Store results
  
  results = list()
  results[["ATT"]] = data.frame(RR=RR, RR.lower=RR.lower, RR.upper=RR.upper)
  results[["Estimates"]] = estdata
  results[["cscm.optim.res.cross.fitted"]] = cscmlist
  results[["dataprep.cross.fitted"]] = datapreplist
  results[["unit.weight.cross.fitted"]] = cross.fitted.w
  results[["ATTs.cross.fitted"]] = augatts
  if (isTRUE(full.model)==TRUE) {
    results[["unit.weight.full.sample"]] = solution.w
    results[["unit.weight.SCM"]] <- synth.w.main
    results[["dataprep.main"]] = dataprep.full
  }
  return(results)
}

# Helper functions for simulation study

factor.init <- function(data,t_int,id,time,outcome,trvar="D",r=3) {
  
  mod <- gsynth::gsynth(data=data,
                        Y=outcome,
                        D=trvar,
                        X=NULL,
                        r=r,
                        index=c(id,time),
                        force=c("two-way"),
                        se=FALSE,
                        normalize=FALSE,
                        CV=FALSE)
  
  
  simparams = list()
  
  simparams[["factors"]] <- mod$factor
  simparams[["residuals"]] <- mod$res.co
  simparams[["loadings"]] <- rbind(mod$lambda.co, mod$lambda.tr)
  simparams[["unitfe"]] <- rbind(mod$alpha.co, mod$alpha.tr)
  simparams[["timefe"]] <- mod$xi
  
  return(simparams)
  
}

factor.sim <- function(factor.init, n=NULL, ar_param=NULL) {
  
  require(MASS)
  
  if (is.null(n)) {n <- nrow(factor.init$unitfe)}
  
  # Prepare fixed parameters
  
  total_t <- nrow(factor.init$timefe)
  cov_loadings <- cov(factor.init$loadings)
  meanLoadings <- colMeans(factor.init$loadings)
  sd_resid <- sd(as.vector(factor.init$residuals))
  timefe <- factor.init$timefe
  factors <- factor.init$factors
  mean_alpha <- mean(factor.init$unitfe)
  sd_alpha <- sd(factor.init$unitfe)
  TIME <- 1:total_t
  
  resid_transp <- t(factor.init$residuals)
  cov_resids <- cov(resid_transp)
  
  datalist = list()
  treatlist = list()
  # Generate outcomes per unit
  
  for (i in 1:n) {
    
    ID <- rep(i, total_t)
    alpha <- rnorm(total_t, mean_alpha, sd_alpha)
    lambda <- mvrnorm(n = 1, meanLoadings, cov_loadings)
    
    if(!is.null(ar_param)==TRUE) {
      tshock <- arima.sim(list(ar=ar_param), n=total_t,
                                        sd=sd_resid)
    } else {
      tshock <- mvrnorm(n = 1, rep(0, ncol(resid_transp)), cov_resids)
    }
    
    #tshock <- rnorm(total_t, 0, sd_resid)
    
    Y <- alpha + timefe + (factors %*% lambda) + tshock
    
    pr_tr <- plogis((1/2) * (scale(alpha) + sum(scale(lambda))))
    
    treatlist[[i]] <- as.data.frame(cbind(i, pr_tr))
    
    datalist[[i]] <- as.data.frame(cbind(ID,TIME,tshock,exp(Y)))
    colnames(datalist[[i]]) <- c("ID","TIME","tshock","Y")
    
  }
  
  # Select one unit for "treatment"
  
  tr.out <- as.data.frame(do.call("rbind",treatlist))
  tr.out[,2] <- tr.out[,2] / sum(tr.out[,2])
  tr_id <- as.numeric(tr.out[which.max(tr.out[,2]),][,1])
  
  data.out <- as.data.frame(do.call("rbind",datalist))
  
  results <- list()
  results[["simdata"]] <- data.out
  results[["treated.id"]] <- tr_id
  results[["controls.id"]] <- unique(data.out$ID)[!(unique(data.out$ID) %in% tr_id)]
  return(results)
}

