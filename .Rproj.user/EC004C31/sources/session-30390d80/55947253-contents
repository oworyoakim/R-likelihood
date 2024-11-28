# Dependencies
library(ggplot2)
library(dplyr)

#---- MSingarch: Specify an MS-INGARCH model
# m: number of regimes
# par: list with parameters. Must contain vectors a and b of length m, and an mxm TPM named gamma

MSingarch <- function(m,
                      par,
                      mean_spec = "linear"){

  obj <- list(m = m,
              par = par,
              mean_spec = mean_spec)
  class(obj) <- "MSingarch"
  return(obj)

}

#---- rINGARCH Generates a random sample from a MS-INGARCH model

# n: The sample size
# model: An MSingarch object with model specifications
# init: list of initial values of the timeseries. Must contain y0 and lambda0
# model: which model to be used. Either "linear" (Default) or "log-linear".
# stationary_mc: Logical value: Should the stationary distribution of the Markov-Chain be used as initial distribution?

rMSingarch  <- function(n,
                        model,
                        init_mean = NULL,
                        init_mc = NULL){

  # Extract needed objects from lists
  m <- model$m
  a <- model$par$a
  b <- model$par$b
  d <- model$par$d
  gamma <- model$par$gamma


  # Create empty objects
  mvec <- seq(m)
  y <- regime <- lambda <- rep(NA, n)
  lambdavec <- matrix(nrow = n, ncol = m)

  # Stationary distribution of MC or specified P(S1)
  if(is.null(init_mc)){
    delta <- solve(t(diag(m) - gamma + 1), rep(1,m))
  }else{delta = init_mc$delta}

  # Linear model
  if(model$mean_spec == "linear"){

    # How should the process be initiated?
    if(is.null(init_mean)){
      y0 <- sum((delta*d)/(1 - a))/(1 - sum((delta*b)/(1 - a))) #Update when stationary mean is found
      lambdavec0 <- rep(y0, m)
      init_mean <- list(y0 = y0,
                        lambdavec0 = lambdavec0)
    }else{
      y0 <- init_mean$y0
      lambdavec0 <- init_mean$lambdavec0
    }

    # for t = 1
    regime[1] <- sample(mvec, 1, prob = delta)
    lambdavec[1, ] <- d + a*lambdavec0 + b*y0
    lambda[1] <- lambdavec[1, regime[1]]
    y[1] <- rpois(1, lambda = lambda[1])

    # for t = 2, ..., n
    for(t in 2:n){
      regime[t] <- sample(m, 1, prob = gamma[regime[t - 1], ])
      # regime[t] <- sample(m, 1, prob = delta)
      lambdavec[t, ] <- d + a*lambdavec[t - 1, ] + b*y[t - 1]
      lambda[t] <- lambdavec[t, regime[t]]
      y[t] <- rpois(1, lambda = lambda[t])
    }
  }

  # Log-linear model
  if(model$mean_spec == "log-linear"){

    # How should the process be initiated?
    if(is.null(init_mean)){
      y0 <- sum(delta*exp(d/(1 - a - b)))  #Update when stationary mean is found
      lambdavec0 <- rep(y0, m)
      init_mean <- list(y0 = y0,
                        lambdavec0 = lambdavec0)
    }else{
      y0 <- init_mean$y0
      lambdavec0 <- init_mean$lambdavec0
    }

    eta <- matrix(NA, nrow = n, ncol = m)

    # for t = 1
    regime[1] <- sample(mvec, 1, prob = delta)
    eta[1, ] <- d + a*log(lambdavec0) + b*log(y0 + 1)
    lambdavec[1, ] <- exp(eta[1, ])
    lambda[1] <- lambdavec[1, regime[1]]
    y[1] <- rpois(1, lambda = lambda[1])

    # for t = 2, ...,n
    for(t in 2:n){
      regime[t] <- sample(m, 1, prob = gamma[regime[t - 1], ])
      eta[t, ] <- d + a*eta[t - 1, ] + b*log(y[t - 1] + 1)
      lambdavec[t, ] <- exp(eta[t, ])
      lambda[t] <- lambdavec[t, regime[t]]
      y[t] <- rpois(1, lambda = lambda[t])
    }
  }

  # Gather results
  res <- list(y = y,
              lambda = lambda,
              regime = regime,
              lambdavec = lambdavec,
              n = n,
              delta = delta,
              init_mean = init_mean,
              init_mc = init_mc,
              model = model)

  class(res) <- "rMSingarch"

  return(res)

  }



## ---- fitMSingarch (fix y0, lambda0 when we no stationary mean)
# Estimation using TMB
fitMSingarch <- function(data,
                         model,
                         init_par = NULL,
                         init_mean = list(y0 = NULL, lambda0 = NULL),
                         ctrl = list(map = list(),
                                     gradient = FALSE,
                                     hessian = FALSE)){

  # TMB data-input
  ldata <- list(y = data,
                m = model$m,
                lambda0 = init_mean$lambdavec0,
                y0 = init_mean$y0)

  # Take init_par values from model spec if it is not specified
  if(is.null(init_par)){
    init_par = list()
    init_par$a = model$par$a
    init_par$b = model$par$b
    init_par$d = model$par$d
    init_par$gamma = model$par$gamma
  }

  # Define MakeADfun object according to model

if(model$mean_spec == "linear"){

    parameters <- list(ta = log(init_par$a),
                       tb = log(init_par$b),
                       td = log(init_par$d),
                       tgamma = Gamma_n2w(m = model$m, gamma = init_par$gamma))
    obj <- MakeADFun(ldata,
                     parameters,
                     DLL = "msingarch_linear",
                     silent = TRUE,
                     map = ctrl$map)
}
  if(model$mean_spec == "log-linear"){

    parameters <- list(a = init_par$a,
                       b = init_par$b,
                       d = init_par$d,
                       tgamma = Gamma_n2w(m = model$m, gamma = init_par$gamma),
                       Beta = init_par$Beta)
    obj <- MakeADFun(ldata,
                     parameters,
                     DLL = "msingarch_log_linear",
                     silent = TRUE,
                     map = ctrl$map)
  }

  # Note: should write something that throws an error if non of the
  # above mean_specs are specified

  # The function ifelse cannot return a NULL value, the function switch can
  # If gradient is FALSE, then gradient + 1 is 1 and the switch returns NULL
  # If gradient is TRUE, then gradient + 1 is 2 and the switch returns mod$gr
  gr <- switch(ctrl$gradient + 1, NULL, obj$gr)
  he <- switch(ctrl$hessian + 1, NULL, obj$he)


  # Optimizing while handling errors
  mod <- tryCatch({
    nlminb(start = obj$par, objective = obj$fn, gradient = gr, hessian = he)
  },
  error = function(e) {
    message("nlminb error:")
    message("m = ", model$m)
    message("gradient = ", ctrl$gradient)
    message("hessian = ", ctrl$hessian)
    message("The original error message is:\n", e)
    return()
  })
  convergence <- mod$convergence
  if (is.null(mod)) {
    return()
  }
  if (convergence != 0) {
    w <- paste0("nlminb didn't succesfully converge:\n",
                mod$message,
                "\nm = ", m,
                "\ngradient = ", gradient,
                "\nhessian = ", hessian, "\n")
    warning(w)
  }



  # Return standard errors
    adrep <- summary(sdreport(obj), "report")

  # Decoding
    # Retrieve the objects at ML value
    delta <- obj$env$report()$delta
    gamma <- obj$env$report()$gamma
    emission_probs <- obj$env$report()$emission_probs
    n <- length(data$y)
    m <- length(delta)
    mllk <- obj$env$report()$mllk

    # Compute log-forward probabilities (scaling used)
    lalpha <- matrix(NA, m, n)
    foo <- delta * emission_probs[1, ]
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, 1] <- log(foo) + lscale
    for (i in 2:n) {
      foo <- foo %*% gamma * emission_probs[i, ]
      sumfoo <- sum(foo)
      lscale <- lscale + log(sumfoo)
      foo <- foo / sumfoo
      lalpha[, i] <- log(foo) + lscale
    }

    # Compute log-backwards probabilities (scaling used)
    lbeta <- matrix(NA, m, n)
    lbeta[, n] <- rep(0, m)
    foo <- rep (1 / m, m)
    lscale <- log(m)
    for (i in (n - 1):1) {
      foo <- gamma %*% (emission_probs[i + 1, ] * foo)
      lbeta[, i] <- log(foo) + lscale
      sumfoo <- sum(foo)
      foo <- foo / sumfoo
      lscale <- lscale + log(sumfoo)
    }

    # Compute conditional state probabilities, smoothing probabilities
    stateprobs <- matrix(NA, ncol = n, nrow = m)
    llk <- - mllk
    for(i in 1:n) {
      stateprobs[, i] <- exp(lalpha[, i] + lbeta[, i] - llk)
    }

    # Most probable states
    ldecode <- rep(NA, n)
    for (i in 1:n) {
      ldecode[i] <- which.max(stateprobs[, i])
    }

    # In-sample predictions using regime-probabilities (not really correct)
    lambda_mat <- obj$env$report()$lambda_mat
    lambda_pred <- rowSums(t(stateprobs)*lambda_mat)


    # In-sample fit statistics
    mllk <- mod$objective
    np <- length(unlist(parameters))
    n <- sum(!is.na(data))
    residuals <- data$y - lambda_pred
    residuals_sc <- (data$y - lambda_pred)

    AIC <- 2 * (mllk + np)
    BIC <- 2 * mllk + np * log(n)
    MSE <- mean(residuals^2)
    MSE_sc <- sum(residuals_sc)/(n - np)


  # Gather results
    ret <- list(data = data,
                model = model,
                par = adrep,
                TMBobj = obj,
                nlminobj = mod,
                convergence = convergence,
                loglik = -mllk,
                AIC = AIC,
                BIC = BIC,
                MSE = MSE,
                m = m,
                MSE_sc = MSE_sc,
                lambda_pred = lambda_pred,
                lalpha = lalpha,
                lbeta = lbeta,
                regimeprobs = t(stateprobs),
                ldecode = ldecode,
                lambda_mat = lambda_mat
                )

    class(ret) <- "fitMSingarch"

    return(ret)

}



#-- Gamma_w2n: Function to transform working TPM parameters to natural
# m: Number of regimes
# tgamma: Working parameters

Gamma_w2n <- function(m, tgamma){

  gamma <- diag(m)
  if (m == 1) return(gamma)

  gamma[!gamma] <- exp(tgamma)
  gamma <- gamma/apply(gamma, 1, sum)

  return(gamma)
}

#-- Gamma_n2w: Function to transform natural TPM parameters to working parameters
# m: Number of regimes
# gamma: Transition probability matrix

Gamma_n2w <- function(m, gamma){

  foo <- log(gamma / diag(gamma))

  tgamma <- as.vector(foo[!diag(m)])

  return(tgamma)
}



#--- plot.rMSingarch: Plot function of an object of class rMSingarch

plot.rMSingarch <- function(obj, ...){

  # Organize data in dataframe
  df <- data.frame(y = obj$y, t = seq(obj$y),
                   lambda = obj$lambda, regime = obj$regime)

  p <- ggplot(df, aes(x = t, y = y)) +
    geom_rect(aes(xmin = t, xmax = dplyr::lead(t),
                  ymin = -Inf, ymax = Inf, fill = factor(regime)), alpha = 0.3) +
    geom_line(alpha = 0.5) +
    geom_line(aes(x = t, y = lambda), color = "red", size = 0.5)

  p$labels$fill <- "Underlying regime"

  p

}


#--- plot.rMSingarch: Plot function of an object of class rMSingarch

plot.fitMSingarch <- function(obj, ...){

  # Organize data in dataframe
  df <- data.frame(y = obj$data$y, t = seq(obj$data$y),
                   lambda = obj$lambda_pred, regime = obj$ldecode)

  p <- ggplot(df, aes(x = t, y = y)) +
    geom_rect(aes(xmin = t, xmax = dplyr::lead(t),
                  ymin = -Inf, ymax = Inf, fill = factor(regime)), alpha = 0.3) +
    geom_line(alpha = 0.5) +
    geom_line(aes(x = t, y = lambda), color = "red", size = 0.5)

  p$labels$fill <- "Most probable regime"

  p

}

#--- print.fitMSingarch

print.fitMSingarch <- function(x, ...){



  cat("Model:\n")
  cat(paste(x$m, "state MS-INGARCH model with", x$model$mean_spec, "mean specification"))
  cat("\nCoefficients:\n")
  print(x$par)
  cat("\n In-sample fit:\n")
  cat("AIC:", x$AIC, "BIC:", x$BIC, "MSE:", x$MSE, "MSE_sc:", x$MSE_sc)

}
