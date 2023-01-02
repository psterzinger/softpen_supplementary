## Distributed as part of the supplementary material for the manuscript
## "Maximum softly-penalized likelihood for mixed effects logistic regression"
##
## Authors: Philipp Sterzinger, Ioannis Kosmidis
## Date: 1 June 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!
mv_loglik <- function(par, data, nAGQ = 1, lfun = NULL) {
  ## glmer parameterization: c(L,beta), L: lower triangular cholesky factor of sigma
  if (is.null(lfun))
    lfun <- mv_infer_loglik(data, nAGQ)
  out <- try(lfun(par), silent = TRUE)
  if (is(out, "try-error")) NA else out
}

mv_infer_loglik <- function(data, nAGQ = 1) {
  data$quadtable <- NULL
  p <- ncol(data$X)
  q <- max(ncol(data$Z), 1)
  dfun <- glmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
  ## Change to logsigma parameterization and fixef first
  function(par) -dfun(pars_w2n(par, p, q)) / 2
}

pars_w2n <- function(par, p, q) {
  par_fe <- par[1:p]
  par_re <- par[(p + 1):length(par)]
  inds <- vapply(1:q, function(i) {2 + (q + 2) * (i - 1) - i * (i + 1) / 2}, 1)
  par_re[inds] <- exp(par_re[inds])
  return(c(par_re, par_fe))
}

pars_n2w <- function(beta, par_re, q) {
  inds <- vapply(1:q, function(i) { 2 + (q + 2) * (i - 1) - i * (i + 1) / 2}, 1)
  par_re[inds] <- log(par_re[inds])
  return(c(beta, par_re))
}

## Numerically evaluate score function
mv_scores <- function(par, data, nAGQ = 1, lfun) {
  numDeriv::grad(mv_loglik, par, data = data, nAGQ = nAGQ, lfun = lfun)
}

nHuber <- function(x, a) {
  ifelse(abs(x) > a,
         - a * (abs(x) - a / 2),
         - x^2 / 2)
}

mv_penalty <- function(par, data, mult = c(0, 0), pen_log_sigma = NULL) {
  X <- data$X
  Y <- data$Y
  M <- data$M
  grouping <- data$grouping
  npar <- length(par)
  p <- ncol(X)
  q <- max(ncol(data$Z), 1)
  beta <- par[1:p]
  par_re <- par[(p + 1):npar]
  fetas <- drop(X %*% beta)
  probs <- plogis(fetas)
  w <- M * probs * (1 - probs)
  pen <- mult[1] * log(det(crossprod(X * sqrt(w)))) / 2 # a multiple of jeffreys
  if (!is.null(pen_log_sigma)) {
    pen <- pen + mult[2] * pen_log_sigma(par_re)
  }
  pen
}

mv_glmm_fit <- function(start, data,
                        mult = c(0, 0),
                        nAGQ = 1,
                        pen_log_sigma = NULL, ...) {
  ## starting vals are in working parameterization: c(beta,L_tilde), L_tilde is lower triangular cholesky of sigma with log-transform of diagonal entries 
  ## null multiplier means we get default multiplier
  if (is.null(mult)) {
    mult <- rep(mv_multiplier(data), 2)
  }
  lfun <- mv_infer_loglik(data, nAGQ = nAGQ)
  obj <- function(pars) {
    pen <- mv_penalty(pars, data,
                      mult = mult,
                      pen_log_sigma = pen_log_sigma)
    -mv_loglik(pars, data, nAGQ = nAGQ, lfun = lfun) - pen
  }
  out <- optimx(par = start,
                fn = obj,
                ...)
  out <- summary(out, order = value)
  attr(out, "data") <- data
  attr(out, "nAGQ") <- nAGQ
  attr(out, "lfun") <- if (is.null(lfun)) NULL else lfun
  attr(out, "mult") <- mult
  class(out) <- c("glmm_fit", class(out))
  attr(out, "call") <- match.call()
  out
}

mv_multiplier <- function(data) {
  n <- nrow(data$X)
  p <- ncol(data$X)
  return(2 * sqrt(p / n))
}

mv_get_glmer <- function(start = NULL, data, nAGQ = 1, method_name = NULL) {
  p <- ncol(data$X)
  q <- max(ncol(data$Z), 1)
  npar <- p + 0.5 * q * (q + 1)
  data$quadtable <- NULL
  out <- try({
    if (is.null(start)) {
      mod <- glmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                   family = binomial(), nAGQ = nAGQ)
    }
    else {
      ## Start value parameterization: c(beta,L_tilde),
      ## where L_tilde: lower triangular entries of lower Cholesky factor of sigma, diag is in logs 
      fixef <- start[1:p]
      theta <- start[(p + 1):npar]
      inds <- vapply(1:q, function(i) {2 + (q + 2) * (i - 1) - i * (i + 1) / 2}, 1)
      theta[inds] <- exp(theta[inds])
      start_vals <- list(theta = theta, fixef = fixef)
      mod <- glmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                   family = binomial(), nAGQ = nAGQ, start = start_vals)
    }
    cfe <- pars_n2w(mod@beta, mod@theta, q)
    bd_diagnostics <- transform_bd_grad(mod,p,q) 
    bd <- bd_diagnostics$bd
    grad <- bd_diagnostics$grad
    ## Get deviance for unpenalized model to compute standard errors
    dfun <- glmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                  family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
    ## Change to logsigma parameterization and fixef first
    nll <- function(par) dfun(pars_w2n(par, p, q)) / 2
    sfe <- sqrtdiaginv(numDeriv::hessian(nll, cfe))
  }, silent = TRUE)
  if (is(out, "try-error")) {
    grad <- bd <- cfe <- sfe <- rep(NA, npar)
  }
  out <- data.frame(estimate = cfe,
                    se = sfe,
                    grad = grad,
                    bd = bd,
                    method = ifelse(is.null(method_name), "glmer", method_name),
                    parameter = factor(c(colnames(data$X), ind_mat(q)), levels = c(colnames(data$X), ind_mat(q)), ordered = TRUE))
  class(out) <- c("get_results_output", class(out))
  attr(out, "call") <- match.call()
  out
}

sqrtdiaginv <- function(mat) {
  out <- try(sqrt(diag(solve(mat))), silent = TRUE)
  if (inherits(out, "try-error")) {
    out <- rep(NA, nrow(mat))
  }
  out
}

hess <- function(object, best = TRUE) {
  details <- attr(object, "details")
  if (is.null(dim(details))) {
    out <- details["nhatend"]
    names(out) <- details$method
    if (best) out[[1]] else out
  }
  else {
    if (best) details[, "nhatend"][[1]] else details[, "nhatend"]
  }
}

transform_bd_grad <- function(mod,p,q){
  nrep <- q*(q+1)/2
  npar <- nrep + p
  
  H <- mod@optinfo$derivs$Hessian/2
  grad <- mod@optinfo$derivs$gradient/2
  
  ## get grad/hessian wrt log transformed pars 
  log_sigma_inds <- vapply(1:q, function(i) {2 + (q + 2) * (i - 1) - i * (i + 1) / 2}, 1)
  grad[log_sigma_inds] <- grad[log_sigma_inds] *  mod@theta[log_sigma_inds]
  
  H[,log_sigma_inds] <-  H[,log_sigma_inds] *  mod@theta[log_sigma_inds]
  H[log_sigma_inds,] <-  H[log_sigma_inds,] *  mod@theta[log_sigma_inds]
  H[cbind(log_sigma_inds,log_sigma_inds)] <- H[cbind(log_sigma_inds,log_sigma_inds)] + grad[log_sigma_inds]
  
  
  
  ## reorder hessian grad from c(RE,FE) to c(FE,RE) 
  grad_reorder <- c(grad[(nrep+1):npar],grad[1:nrep])
  
  H_reorder <- H
  H_reorder[1:p,1:p] <- H[(nrep+1):npar,(nrep+1):npar]
  H_reorder[(p+1):npar,(p+1):npar] <- H[1:nrep,1:nrep]
  H_reorder[1:p,(p+1):npar] <- H[(nrep+1):npar,1:nrep]
  H_reorder[(p+1):npar,1:p] <- H[1:nrep,(nrep+1):npar]
  
  return(list(bd = sqrtdiaginv(H_reorder), grad = grad_reorder))
}

mv_get_MSPAL <- function(start, data, nAGQ = 1,
                         mult = NULL, pen_log_sigma = NULL,
                         optimization_methods = c("BFGS", "CG"),
                         method_suffix = "", method_name = NULL, ...) {
  p <- ncol(data$X)
  q <- max(ncol(data$Z), 1)
  npar <- p + 0.5 * q * (q + 1)
  out <- try({
    fit0 <- mv_glmm_fit(start,
                        data = data, nAGQ = nAGQ,
                        method = optimization_methods,
                        mult = mult,
                        pen_log_sigma = pen_log_sigma,
                        ...)
    mult <- attr(fit0, "mult")
    cfe <- coef.glmm_fit(fit0, best = TRUE)
    sfe <- se(fit0, best = TRUE)
    grad <- gradient(fit0, best = TRUE)
    bd <- boundary_diagnostics(fit0, best = TRUE)
  }, silent = TRUE)
  if (is(out, "try-error")) {
    grad <- bd <- cfe <- sfe <- rep(NA, npar)
  }
  if (is.null(method_name)) {
    method_name <- paste0("MSPAL_", method_suffix, "[", paste(format_dec(mult, 1), collapse = ","), "]")
  }
  out <- data.frame(estimate = cfe,
                    se =  sfe,
                    grad = grad,
                    bd = bd,
                    method = method_name,
                    parameter = factor(c(colnames(data$X), ind_mat(q)), levels = c(colnames(data$X), ind_mat(q)), ordered = TRUE))
  attr(out, "mult") <- mult
  class(out) <- c("get_results_output", class(out))
  attr(out, "call") <- match.call()
  out
}

se <- function(object, best = TRUE) {
  if (best) {
    sqrtdiaginv(info(object, best = TRUE))
  }
  else {
    do.call("rbind", lapply(info(object, best = FALSE), sqrtdiaginv))
  }
}

info <- function(object, best = TRUE) {
  data <- attr(object, "data")
  nAGQ <- attr(object, "nAGQ")
  lfun <- attr(object, "lfun")
  info_fun <- function(x) {
    -numDeriv::hessian(mv_loglik, x, data = data, nAGQ = nAGQ, lfun = lfun)
  }
  if (best) {
    info_fun(coef(object, best = TRUE))
  }
  else {
    apply(coef(object, best = FALSE), 1, info_fun, simplify = FALSE)
  }
}

gradient <- function(object, best = TRUE) {
  details <- attr(object, "details")
  if (is.null(dim(details))) {
    out <- details["ngatend"]
    names(out) <- details$method
    out[[1]]
  }
  else {
    if (best) details[, "ngatend"][[1]] else do.call("rbind", details[, "ngatend"])
  }
}

boundary_diagnostics <- function(object, best = TRUE) {
  if (best) {
    sqrtdiaginv(hess(object, best = TRUE))
  }
  else {
    lapply(hess(object, best = FALSE), sqrtdiaginv)
  }
}

format_dec <- function(x, k) format(round(x, k), trim = TRUE, nsmall = k)

mv_get_bglmer <- function(start = NULL, data, nAGQ = 1,
                          c_prior, f_prior, method_name = NULL) {
  p <- ncol(data$X)
  q <- max(ncol(data$Z), 1)
  npar <-  p + 0.5 * q * (q + 1)
  out <- try({
    if (is.null(start)) {
      mod <- bglmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                    family = binomial, nAGQ = nAGQ,
                    cov.prior = c_prior, fixef.prior = f_prior)
    }
    else {
      ## starting vals are in working parameterization: c(beta,L_tilde), L_tilde is lower triangular cholesky of sigma with log-transform of diagonal entries 
      fixef <- start[1:p]
      theta <- start[(p + 1):npar]
      inds <- vapply(1:q, function(i) {2 + (q + 2) * (i - 1) - i * (i + 1) / 2}, 1)
      theta[inds] <- exp(theta[inds])
      start <- list(fixef = fixef, theta = theta)
      mod <- bglmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                    family = binomial, nAGQ = nAGQ, start = start,
                    cov.prior = c_prior, fixef.prior = f_prior)
    }
    cfe <- pars_n2w(mod@beta, mod@theta, q)
    bd_diagnostics <- transform_bd_grad(mod,p,q) 
    bd <- bd_diagnostics$bd
    grad <- bd_diagnostics$grad
    ## Get deviance for unpenalized model to compute standard errors
    dfun <- glmer(Y ~ -1 + X + (-1 + Z | grouping), weights = M, data = data,
                  family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
    ## Change to logsigma parameterization and fixef first
    nll <- function(pars) dfun(pars_w2n(pars, p, q)) / 2
    sfe <- sqrtdiaginv(numDeriv::hessian(nll, cfe))
  }, silent = FALSE)
  if (is(out, "try-error")) {
    grad <- bd <- cfe <- sfe <- rep(NA, npar)
  }
  out <- data.frame(estimate = cfe,
                    se = sfe,
                    grad = grad,
                    bd = bd,
                    method = ifelse(is.null(method_name), paste0("BGLMER[", as.character(f_prior), ",", as.character(c_prior), "]"), method_name),
                    parameter = factor(c(colnames(data$X), ind_mat(q)), levels = c(colnames(data$X), ind_mat(q)), ordered = TRUE))
  class(out) <- c("get_results_output", class(out))
  attr(out, "call") <- match.call()
  out
}

coef.glmm_fit <- function(object, best = TRUE) {
  class(object) <- class(object)[-1]
  cc <- coef(object)
  if (best) drop(cc[1, ]) else  cc
}

## Helper for progress reporting in simulations
report_res <- function(out, digits = 2) {
  est <- format_dec(out$estimate, digits)
  bd <- format_dec(out$bd, digits)
  mm <- as.character(unique(out$method))
  nc <- nchar(mm)
  cat(paste0(mm, ":"), "\t\t")
  cat(paste0(est, " [", bd, "] "), "\n")
}

## getSummary helper functions to interface memisc::mtable
getSummary.get_results_output <- function(x) {
  coefs <- x[, "estimate"]
  ses <- x[, "se"]
  stat <- coefs / ses
  p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  cc <- cbind(coefs, ses, stat, p)
  colnames(cc) <- c("est", "se", "stat", "p")
  rownames(cc) <- x[, "parameter"]
  list(coef = cc, call = attr(x, "call"))
}

getSummary.glmm_fit <- function(x) {
  coefs <- coef(x)
  ses <- se(x)
  stat <- coefs / ses
  p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  cc <- cbind(coefs, ses, stat, p)
  colnames(cc) <- c("est", "se", "stat", "p")
  rownames(cc) <- c(colnames(attr(x, "data")$X), ind_mat(max(ncol(attr(x, "data")$Z), 1)))
  list(coef = cc, call = attr(x, "call"))
}

## Functions for simulation experiments

mv_simufun <- function(par, data, nsimu, seed = NULL) {
  # par takes logl parameterization of sigma
  X <- data$X
  M <- data$M
  Z <- data$Z
  grouping <- data$grouping
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  p <- ncol(X)
  q <- max(ncol(Z), 1)
  npar <- p + 0.5 * q * (q + 1)
  n <- nrow(X)
  ## Number of groups
  ng <- length(unique(grouping))
  beta <- par[1:p]
  if (q > 1) {
    sigma <- matrix(0, nrow = q, ncol = q)
    sigma[lower.tri(sigma, diag = TRUE)] <- par[(p + 1):npar]
    diag(sigma) <- exp(diag(sigma))
    sigma <- tcrossprod(sigma)
    reff <- vapply(1:nsimu, function(i, ng, sigma, Z, grouping) {u <- MASS::mvrnorm(ng, rep(0, q), sigma);
    return(rowSums(Z * u[as.numeric(grouping), ]))}, rep(1, n),ng = ng, sigma = sigma, Z = Z, grouping = grouping)
  }else{
    sigma <- exp(par[(p + 1)])
    reff <- matrix(rnorm(ng * nsimu, 0, sigma), nrow = ng)
    reff <- reff[as.numeric(grouping), ]
  }
  fetas <- drop(X %*% beta)
  etas <- fetas + reff
  probs <- plogis(etas)
  if (nsimu > 1) {
    val <- apply(probs, 2, function(prob) rbinom(n, M, prob))
  }else{
    val <- rbinom(n, M, probs)
  }
  val
}

mv_fit_all_methods <- function(startl, data, s, trace = TRUE, nAGQ = 1, optimization_methods, c_prior=NULL) {
  for (st in startl) {
    res_mal <- mv_get_MSPAL(start = st, data = data, nAGQ = nAGQ,
                            mult = c(0, 0),
                            optimization_methods = optimization_methods,
                            method_name = "ML")
    if (isTRUE(all(abs(res_mal$grad) < 1e-03))) break
  }
  startl <- c(list(res_mal$estimate), startl)
  for (st in startl) {
    res_mspal <- mv_get_MSPAL(start = st, data = data, nAGQ = nAGQ,
                              pen_log_sigma = function(x) sum(nHuber(x, 1)),
                              optimization_methods = optimization_methods,
                              method_name = "MSPL")
    if (isTRUE(all(abs(res_mspal$grad) < 1e-03))) break
  }
  startl <- c(list(res_mspal$estimate), startl)
  for (st in startl) {
    if (is.null(c_prior)) {
      res_bglmer_t <- mv_get_bglmer(start = st, data = data,  nAGQ = nAGQ,
                                    c_prior = bquote(wishart), f_prior = bquote(t),
                                    method_name = "bglmer[t]")
    }else{
      res_bglmer_t <- mv_get_bglmer(start = st, data = data,  nAGQ = nAGQ,
                                    c_prior = c_prior, f_prior = bquote(t),
                                    method_name = "bglmer[t]")
    }
    
    if (isTRUE(all(abs(res_bglmer_t$grad) < 1e-03))) break
  }
  startl <- c(list(res_bglmer_t$estimate), startl)
  for (st in startl) {
    if (is.null(c_prior)) {
      res_bglmer_n <- mv_get_bglmer(start = st, data = data, nAGQ = nAGQ,
                                    c_prior = bquote(wishart), f_prior = bquote(normal),
                                    method_name = "bglmer[n]")
    }else{
      res_bglmer_n <- mv_get_bglmer(start = st, data = data, nAGQ = nAGQ,
                                    c_prior = c_prior, f_prior = bquote(normal),
                                    method_name = "bglmer[n]")
    }
    if (isTRUE(all(abs(res_bglmer_n$grad) < 1e-03))) break
  }
  if (trace) {
    cat("===================================\n")
    cat("Sample:\t", s, "\n")
    report_res(res_mal)
    report_res(res_mspal)
    report_res(res_bglmer_t)
    report_res(res_bglmer_n)
  }
  rbind(res_mal,
        res_bglmer_n,
        res_bglmer_t,
        res_mspal)
}

mv_perform_experiment <- function(truth, data,
                                  nsimu = 100, seed = NULL,
                                  alt_start = rep(0, length(truth)), nAGQ = 1,
                                  optimization_methods = c("BFGS", "CG"),
                                  ncores = 4,
                                  mathpar = NULL,
                                  c_prior = NULL) {
  registerDoMC(ncores)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  simulated_data <- mv_simufun(truth, data = data, nsimu = nsimu)
  base_startl <- list(truth, alt_start)
  out <- mclapply(1:nsimu, function(j) {
    res <- mv_fit_all_methods(startl = base_startl,
                              data = within(data, Y <- simulated_data[, j]),
                              s = j,
                              trace = TRUE,
                              optimization_methods = optimization_methods,
                              nAGQ = nAGQ,
                              c_prior = c_prior)
    res$sample <- j
    res
  }, mc.cores = ncores)
  #out <- lapply(1:nsimu, function(j) {
  #  res <- mv_fit_all_methods(startl = base_startl,
  #                         data = within(data, Y <- simulated_data[, j]),
  #                         s = j,
  #                         trace = TRUE,
  #                         optimization_methods = optimization_methods,
  #                         nAGQ = nAGQ)
  #  res$sample <- j
  #  res
  #})
  out <- do.call("rbind", out)
  out$true <- truth
  parnames <- c(colnames(data$X), ind_mat(max(ncol(data$Z), 1)))
  if (is.null(mathpar)) {
    mathpar <- parnames
    out$mathpar <- mathpar
  }
  else {
    names(mathpar) <- parnames
    out$mathpar <-  dplyr::recode(out$parameter, !!!mathpar)
  }
  out$mathpar <- factor(out$mathpar, levels = rev(mathpar),
                        ordered = TRUE)
  out$parameter <- factor(out$parameter, levels = rev(parnames),
                          ordered = TRUE)
  out$method <- factor(out$method, levels = rev(c("ML", "bglmer[n]", "bglmer[t]", "MSPL")),
                       ordered = TRUE)
  out
}

mv_summarize_experiment <- function(out, se_threshold = 40, grad_threshold = 1e-03,
                                    estimate_threshold = 50) {
  ## Exclusion
  out <- within(out, estimate[(se > se_threshold) |
                                (abs(grad) > grad_threshold) |
                                abs(estimate) > estimate_threshold] <- NA)
  out %>%
    group_by(method, parameter, mathpar) %>%
    summarize(
      pu = mean(estimate < true, na.rm = TRUE),
      bias = mean(estimate - true, na.rm = TRUE),
      var = var(estimate, na.rm = TRUE),
      cover95 = mean(abs(estimate - true) < qnorm(0.975) * se, na.rm = TRUE),
      mse = var + bias^2,
      mad = mean(abs(estimate - true), na.rm = TRUE),
      mean = mean(estimate, na.rm = TRUE),
      nused = sum(!is.na(estimate)),
      ntried = length(estimate),
      perc_used = nused / ntried * 100,
      min = min(estimate, na.rm = TRUE),
      max = max(estimate, na.rm = TRUE),
      cmin = min(estimate - true, na.rm = TRUE),
      cmax = max(estimate - true, na.rm = TRUE))
}

mv_plot_experiment <- function(out, out_summary = NULL, xlims = c(-10, 10),
                               max_var = 10, max_abs_bias = 3, ...) {
  if (is.null(out_summary)) out_summary <- mv_summarize_experiment(out, ...)
  ## Summary plots
  stats <- c("nused", "bias", "var", "mse", "pu", "cover95")
  xlabs <- c("R", "Bias", "Variance", "MSE", "PU", "Coverage")
  plots <- as.list(numeric(length(stats)))
  names(xlabs) <- names(plots) <- stats
  for (stat in stats) {
    plots[[stat]] <- ggplot(out_summary) +
      geom_col(aes_string(y = "mathpar", x = stat,
                          col = "method", fill = "method"),
               alpha = 0.7, position = "dodge") +
      labs(x = NULL, y = NULL, title = xlabs[stat]) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_y_discrete(labels = label_parse()) +
      scale_color_grey(start = 0.1, end = 0.8) +
      scale_fill_grey(start = 0.2, end = 0.9)
    if (stat %in% c("nused", "var", "mse", "pu", "cover95"))
      plots[[stat]] <- plots[[stat]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank())
    if (stat == "cover95")
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0.95),
                   lty = 2, alpha = 0.7) +
        coord_cartesian(x = c(0.6, 1.00))
    if (stat %in% c("bias", "var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0),
                   lty = 2, alpha = 0.7)
    if (stat %in% c("var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(0, max_var))
    if (stat %in% c("bias"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(-1,  1) * max_abs_bias)
    if (stat %in% "pu")
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0.5),
                   lty = 2, alpha = 0.7)
    if (stat == "nused") {
      ntried <- unique(out_summary$ntried)
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = ntried),
                   lty = 2, alpha = 0.7) +
        scale_x_continuous(breaks = ntried / 5 * (0:5),
                           labels = c("0", paste0(seq(0, ntried / 1000, length = 6), "K")[-1])) +
        theme(axis.text.x = element_text(angle = 45))
    }
  }
  ## boxplots
  bpp <- ggplot(out) +
    geom_boxplot(aes(y = mathpar, x = estimate - true, colour = method, fill = method), alpha = 0.5,
                 outlier.alpha = 0.1, outlier.size = 0.3) +
    geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.7) +
    coord_cartesian(x = xlims)  +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = NULL, title = "Centred distributions") +
    scale_y_discrete(labels = label_parse()) +
    scale_color_grey(start = 0.1, end = 0.8) +
    scale_fill_grey(start = 0.2, end = 0.9)
  ## patchwork
  ((bpp + plots[[1]] + plot_layout(widths = c(4, 1))) /
      (Reduce("+", plots[-1]) + plot_layout(ncol = length(stats) - 1)))
}

mv_boxplot_simul <- function(dat,
                             par_ind = 3,
                             exclude = FALSE,
                             ylims = c(-30, 30),
                             xlab = NULL,
                             ylab = NULL,
                             truelab = NULL,
                             se_threshold = 40,
                             estimate_threshold = 50,
                             grad_threshold = 1e-3,
                             x_breaks = 1,
                             centered = TRUE) {
  if (exclude) {
    simul_dat <- within(dat,
                        estimate[(se > se_threshold) |
                                   (abs(grad) > grad_threshold) |
                                   abs(estimate) > estimate_threshold] <- NA)
    split_list <- split(simul_dat, simul_dat$method)
  }else{
    split_list <- split(dat, dat$method)
  }
  methods <- names(split_list)
  plot_list <- list()
  scaleFUN <- function(x) sprintf("%.1f", as.numeric(x))
  every_nth <- function(n) {
    return(function(x) x[c(TRUE, rep(FALSE, n - 1))])
  }
  for (i in seq_len(length(split_list))) {
    dat <- split_list[[i]]
    dat_select <- dat %>%
      group_by(n, lambda) %>%
      slice(seq(par_ind, n(), by = length(unique(parameter)))) %>%
      dplyr::select(estimate, n, lambda)
    if (centered) {
      df <- data.frame(lambda = factor(round(dat_select$lambda, digits = 2)),
                       estimate = dat_select$estimate - dat_select$lambda,
                       n = factor(dat_select$n),
                       lambdas = dat_select$lambda)
    }else{
      df <- data.frame(lambda = factor(round(dat_select$lambda, digits = 2)),
                       estimate = dat_select$estimate,
                       n = factor(dat_select$n),
                       lambdas = dat_select$lambda)
    }
    bp <- ggplot(df, aes(x = lambda, y = estimate, fill = n)) +
      geom_boxplot() +
      theme_minimal() +
      scale_color_manual(values = "black",
                         name = truelab, labels = "",
                         position = "bottom") +
      scale_fill_grey(start = 0.5, end = 1) +
      labs(title = methods[i], x = xlab, y = ylab) +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(fill = guide_legend(override.aes = list(size = .5), order = 1)) +
      scale_x_discrete(labels = scaleFUN, breaks = every_nth(n = x_breaks)) +
      coord_cartesian(ylim = ylims)
    if (!centered) {
      bp <- bp +
        geom_line(data = df,
                  mapping = aes(x = lambda, y = lambdas, group = 1, col = "black")) +
        geom_point(data = df,
                   aes(x = lambda, y = lambdas, col = "black"),
                   size = 2.5, shape = 23, show.legend = TRUE)
    }else{
      bp <- bp + geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)
    }
    plot_list[[i]] <- bp
  }
  ((plot_list[[1]] + plot_list[[4]]) / (plot_list[[2]] + plot_list[[3]]))+ plot_layout(guides = "collect")
}

mv_simul_table <- function(dat,
                           par_ind = 3,
                           se_threshold = 40,
                           estimate_threshold = 50,
                           grad_threshold = 1e-3,
                           log_estimate_lower_threshold = log(1e-2),
                           log_mask = NULL) {
  out <- dat
  out$method <- ordered(out$method, unique(levels(out$method))[c(1, 4, 2, 3)])
  if (!is.null(log_mask)) {
    log_mask <- rep(log_mask, nrow(out) / length(log_mask))
  }else{
    log_mask <- rep(FALSE, nrow(out))
  }
  out$log_mask <- log_mask
  inds <-  which((out$se > se_threshold) |
                   (abs(out$grad) > grad_threshold) |
                   ((abs(out$estimate) > estimate_threshold) & !out$log_mask) |
                   ((out$estimate > log(estimate_threshold)) & out$log_mask) |
                   ((out$estimate < log_estimate_lower_threshold) & out$log_mask))
  out$estimate[inds] <- NA
  table_dat <- out %>%
    group_by(method, n, lambda) %>%
    slice(seq(par_ind, n(), by = length(unique(parameter)))) %>%
    dplyr::select(estimate, n, lambda) %>%
    summarize(
      nused = sum(!is.na(estimate)),
      ntried = length(estimate),
      perc_not_used = (1 - nused / ntried) * 100) %>%
    dplyr::select(method, n, lambda, perc_not_used)
  table_dat <- tidyr::pivot_wider(table_dat,
                                  id_cols = 1:2,
                                  names_from = 3,
                                  values_from = 4)
  steps <- nrow(table_dat) / length(unique(table_dat$method))
  inds <- rep(FALSE, steps)
  inds[ceiling(steps / 2)] <- TRUE
  inds <- rep(inds, length(unique(table_dat$method)))
  table_dat$method[!inds] <- NA
  names(table_dat) <- c("", "", round(unique(out$lambda), digits = 2))
  table_dat[, 2] <- rep(c("n=50", "n=100", "n=200"), 4)
  table_dat[, 3:ncol(table_dat)] <- t(apply(
    table_dat[, 3:ncol(table_dat)], 1, function(vec) {
      apply(as.matrix(vec), 1, function(entry) {
        if (entry == 0) {
          round(entry, digits = 0)
        }
        else{
          paste("\\textbf{", round(entry, digits = 0), "}", sep = "")
        }
      })
    }))
  return(xtable(table_dat,
                align = c("l", "l", "r",
                          rep("c", ncol(table_dat) - 2)),
                digits = 2))
}

mv_report_outliers <- function(simul_data, ylims, par_ind) {
  simul_data %>% group_by(method) %>%    slice(seq(par_ind, n(), by = length(unique(parameter)))) %>% summarize(outliers = sum(estimate>ylims[2] | estimate<ylims[1],na.rm=TRUE),
                                                                                                                nas = sum(is.na(estimate)) )
}

mv_percentile_table <- function(simul_dat, p, q, var_only = TRUE,
                                perc = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
                                centered = TRUE,
                                se_threshold = 40,
                                grad_threshold = 1e-03,
                                estimate_threshold = 50) {
  if (is.null(simul_dat$mathpar)) {
    simul_dat$mathpar <- simul_dat$parameter
  }
  simul_dat_clean <- within(simul_dat, estimate[(se > se_threshold) |
                                                  (abs(grad) > grad_threshold) |
                                                  abs(estimate) > estimate_threshold] <- NA)
  npar <- p + 0.5 * q * (q + 1)
  if (var_only) {
    inds <- c(npar - 2, npar - 1, npar, 1:(npar - 3))
  }else{
    inds <- c(1:(npar - 3), npar - 2, npar - 1, npar)
  }
  simul_dat_clean$mathpar <- ordered(simul_dat_clean$mathpar,
                                     levels = unique(simul_dat_clean$mathpar)[inds])
  if (centered) {
    tab <- simul_dat_clean %>%
      group_by(method, mathpar, parameter) %>%
      summarise(percentiles = quantile(estimate - true, perc, na.rm = TRUE),
                nused = 1 - sum(is.na(estimate)) / length(estimate))
  }else{
    tab <- simul_dat_clean %>%
      group_by(method, mathpar, true) %>%
      summarise(percentiles = quantile(estimate, perc, na.rm = TRUE),
                nused = 1 - sum(is.na(estimate)) / length(estimate))
  }
  tab$perc <- rep(perc, npar * length(unique(tab$method)))
  
  #nused <- tab$nused[seq.int(from=1,by=length(perc),to = length(tab$nused))]
  tab <- tidyr::pivot_wider(tab,
                            id_cols = c(method, mathpar),
                            names_from = perc,
                            values_from = percentiles)
  #tab$nused <- nused
  
  if (var_only) {
    inds <- rep(c(rep(TRUE, 0.5 * q * (q + 1)),
                  rep(FALSE, p)),
                length(unique(tab$method)))
    tab <- tab[inds, ]
    par_step <- 0.5 * q * (q + 1)
  }else{
    par_step <- p + 0.5 * q * (q + 1)
  }
  steps <- nrow(tab) / length(unique(tab$method))
  inds <- rep(FALSE, steps)
  inds[ceiling(steps / 2)] <- TRUE
  inds <- rep(inds, length(unique(tab$method)))
  tab$method[!inds] <- NA
  colnames(tab) <- c("a", "b", colnames(tab)[3:ncol(tab)])
  
  
  addtorow <- list()
  m_col <- ncol(tab)
  m_col_w <- ncol(tab) - 2
  addtorow$pos <- as.list(c(-1, -1, -1, seq.int(from = 0, to = nrow(tab), by = par_step)))
  addtorow$command <- c("\\toprule \n",
                        paste0("&& \\multicolumn{",m_col_w,"}{c}{Percentiles} \\\\ \\cmidrule{3-", m_col, "} \n", collapse = ""),
                        paste0("&", paste0("&$", 100 * perc, "\\%$", collapse = ""), "\\\\ \\cmidrule{3-", m_col, "} \n", collapse = ""),
                        rep("\\cmidrule{3-9} \n", length(addtorow$pos) - 4),
                        "\\bottomrule \n")
  print(xtable(tab,
               digits = 2,
               align = c("l", "l", "r", rep("c", ncol(tab) - 2))),
        add.to.row = addtorow,
        include.colnames = FALSE,
        include.rownames = FALSE,
        sanitize.text.function = function(x) x,
        hline.after = NULL)
}

re_names_plot <- function(q) {
  mat <- matrix(NA,nrow=q,ncol=q) 
  t = rep(NA,0.5*q*(q+1))
  count <- 1
  for(i in 1:q){
    for(j in 1:i){
      if(i==j){
        t[count] <- paste0('log(l[',i,'*','","','*',j,'])')
      }else{
        t[count] <- paste0('l[',i,'*','","','*',j,']')
      }
      count <- count+1 
    }
  }
  return(t)
}

re_names <- function(q) {
  mat <- matrix(NA,nrow=q,ncol=q) 
  for(i in 1:q){
    for(j in 1:i){
      mat[i,j] <- paste("l_{",paste(i,j,sep=","),"}",sep="")
    }
  }
  diag(mat) <- vapply(diag(mat),function(entry){paste("\\log ",entry,sep="")},"a")
  return(paste("$",mat[lower.tri(mat,diag=T)],"$",sep=""))
}

make_math <- function(mats, ndigits = 2, inline_math = FALSE) {
  mat <- format(round(mats, digits=ndigits), nsmall = ndigits,trim = TRUE)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(i%%2==0){
        if(inline_math){
          mat[i,j] <- paste("$(",mat[i,j],")$",sep="")
        }else{
          mat[i,j] <- paste("(",mat[i,j],")",sep="")
        }
        
      }else{
        if(inline_math){
          mat[i,j] <- paste("$",mat[i,j],"$",sep="")
        }
      }
    }
  }
  mat
}

mv_bar_plot_experiment <- function(out, out_summary = NULL, xlims = c(-10, 10), var_only=TRUE , p = NULL, q = NULL, 
                                   max_var = 10, max_abs_bias = 3, ...) {
  if (is.null(out_summary)) out_summary <- mv_summarize_experiment(out, ...)
  if (var_only & !is.null(p) &!is.null(q)){
    inds <- rep( c(rep(TRUE,0.5*q*(q+1)),rep(FALSE,p)),length(unique(out_summary$method)))
    out_summary <- out_summary[inds,]
  }
  ## Summary plots
  stats <- c("nused", "bias", "var", "mse", "pu")
  xlabs <- c("R", "Bias", "Variance", "MSE", "PU")
  plots <- as.list(numeric(length(stats)))
  names(xlabs) <- names(plots) <- stats
  for (stat in stats) {
    plots[[stat]] <- ggplot(out_summary) +
      geom_col(aes_string(y = "mathpar", x = stat,
                          col = "method", fill = "method"),
               alpha = 0.7, position = "dodge") +
      labs(x = NULL, y = NULL, title = xlabs[stat]) +
      theme_minimal() +
      scale_y_discrete(labels = label_parse()) +
      scale_color_grey(start = 0.1, end = 0.8) +
      scale_fill_grey(start = 0.2, end = 0.9)
    if (stat %in% c("nused", "var", "mse", "pu"))
      plots[[stat]] <- plots[[stat]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank())
    if (stat %in% c("bias", "var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0),
                   lty = 2, alpha = 0.7)
    if (stat %in% c("var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(0, max_var))
    if (stat %in% c("bias"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(-1,  1) * max_abs_bias)
    if (stat %in% "pu")
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0.5),
                   lty = 2, alpha = 0.7)
    if (stat == "nused") {
      ntried <- unique(out_summary$ntried)
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = ntried),
                   lty = 2, alpha = 0.7) +
        scale_x_continuous(breaks = ntried / 5 * (0:5),
                           labels = c("0", paste0(seq(0, ntried / 1000, length = 6), "K")[-1])) +
        theme(axis.text.x = element_text(angle = 45))
    }
  }
  ## patchwork
  (Reduce("+", plots[-1]) +   plots[[1]] + plot_layout(ncol = length(stats),guides = "collect" )) 
}

mv_plot_experiment_color <- function(out, out_summary = NULL, xlims = c(-10, 10),
                                     max_var = 10, max_abs_bias = 3, color_palette=NULL,outline=FALSE, ...) {
  out$method <- factor(out$method, levels = unique(out$method))
  if (is.null(out_summary)) out_summary <- mv_summarize_experiment(out, ...)
  ## Colors
  if (is.null(color_palette)) {
    color_palette <- brewer.pal(n = length(unique(out$method)), name = "PuBuGn")
  }else if(length(color_palette) < length(unique(out$method))){
    color_palette <- brewer.pal(n = length(unique(out$method)), name = color_palette)
  }
  ## Summary plots
  stats <- c("nused", "bias", "var", "mse", "pu", "cover95")
  xlabs <- c("R", "Bias", "Variance", "MSE", "PU", "Coverage")
  plots <- as.list(numeric(length(stats)))
  names(xlabs) <- names(plots) <- stats
  for (stat in stats) {
    plots[[stat]] <- ggplot(out_summary) +
      geom_col(aes_string(y = "mathpar", x = stat,
                          col = "method", fill = "method"),
               alpha = 0.7, position = "dodge") +
      labs(x = NULL, y = NULL, title = xlabs[stat]) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_y_discrete(labels = label_parse()) +
      scale_fill_manual(values = color_palette)
    if(outline){
      plots[[stat]] <- plots[[stat]] +
        scale_color_manual(values = rep("black", length(color_palette)))
    }else{
      plots[[stat]] <- plots[[stat]] + scale_color_manual(values = color_palette)
    }
    
    #+
    # scale_color_grey(start = 0.1, end = 0.8) + scale_fill_grey(start = 0.2, end = 0.9)
    if (stat %in% c("nused", "var", "mse", "pu", "cover95"))
      plots[[stat]] <- plots[[stat]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank())
    if (stat == "cover95")
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0.95),
                   lty = 2, alpha = 0.7) +
        coord_cartesian(x = c(0.6, 1.00))
    if (stat %in% c("bias", "var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0),
                   lty = 2, alpha = 0.7)
    if (stat %in% c("var", "mse"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(0, max_var))
    if (stat %in% c("bias"))
      plots[[stat]] <- plots[[stat]] +
        coord_cartesian(x = c(-1,  1) * max_abs_bias)
    if (stat %in% "pu")
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = 0.5),
                   lty = 2, alpha = 0.7)
    if (stat == "nused") {
      ntried <- unique(out_summary$ntried)
      plots[[stat]] <- plots[[stat]] +
        geom_vline(aes(xintercept = ntried),
                   lty = 2, alpha = 0.7) +
        scale_x_continuous(breaks = ntried / 5 * (0:5),
                           labels = c("0", paste0(seq(0, ntried / 1000, length = 6), "K")[-1])) +
        theme(axis.text.x = element_text(angle = 45))
    }
  }
  ## boxplots
  bpp <- ggplot(out) +
    geom_boxplot(aes(y = mathpar, x = estimate - truth, colour = method, fill = method), alpha = 0.5, 
                 outlier.alpha = 0.1, outlier.size = 0.3) +
    geom_vline(aes(xintercept = 0), lty = 2, alpha = 0.7) +
    coord_cartesian(x = xlims)  +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = NULL, title = "Centred distributions") +
    scale_y_discrete(labels = label_parse()) + 
    scale_fill_manual(values=color_palette) + 
    guides( color= guide_legend(reverse = TRUE), fill=guide_legend(reverse=TRUE) )
  if(outline){
    bpp <- bpp +
      scale_color_manual(values=rep("black",length(color_palette)))
  }else{
    bpp <- bpp + scale_color_manual(values=color_palette)
  }
  #+
  #  scale_color_grey(start = 0.1, end = 0.8) + scale_fill_grey(start = 0.2, end = 0.9)
  ## patchwork
  ((bpp + plots[[1]] + plot_layout(widths = c(4, 1))) /
      (Reduce("+", plots[-1]) + plot_layout(ncol = length(stats) - 1)))
}

ind_mat <- function(q){
  if(q>1){
    mat <- matrix(NA,nrow=q,ncol=q) 
    for(i in 1:q){
      for(j in 1:i){
        mat[i,j] <- paste("l[",paste(i,j,sep=","),"]",sep="")
      }
    }
    
    diag(mat) <- vapply(diag(mat),function(entry){paste("log(",entry,")",sep="")},"a")
    return(mat[lower.tri(mat,diag=T)])
  }else{
    return("log(sigma)")
  }
}

mv_fit_all <- function(start, data, nAGQ = 1, results.only = TRUE) {
  m_MAL <- mv_glmm_fit(start, mult = c(0, 0), data = data, nAGQ = nAGQ,
                       method = c("BFGS", "CG"))
  m_bglmer_t <- mv_get_bglmer(start, data = data, nAGQ = nAGQ,
                              f_prior = bquote(t), c_prior = bquote(wishart))
  m_bglmer_n <- mv_get_bglmer(start, data = data, nAGQ = nAGQ,
                              f_prior = bquote(normal), c_prior = bquote(wishart))
  m_MSPAL <- mv_get_MSPAL(start, data = data, nAGQ = nAGQ,
                          pen_log_sigma = function(x) sum(nHuber(x, 1)))
  return( list(ML = m_MAL, bglmer_t = m_bglmer_t, bglmer_n = m_bglmer_n, MSPL = m_MSPAL) )
  if (results.only) {
    return( list(ML = m_MAL, bglmer_t = m_bglmer_t, bglmer_n = m_bglmer_n, MSPL = m_MSPAL) )
  }else{
    tab <- memisc::mtable("BFGS" = m_MAL["BFGS", ],  "CG" = m_MAL["CG", ],
                          "BGLMER(t)" = m_bglmer_t, "BGLMER(n)" = m_bglmer_n,
                          "MSPL" = m_MSPAL, digits = 2, #coef.style = "horizontal",
                          signif.symbols = NULL, summary.stats = NULL)
    print(tab)
    return( list(mtable = tab, ML = m_MAL, bglmer_t = m_bglmer_t, bglmer_n = m_bglmer_n, MSPL = m_MSPAL) ) 
  }
}

mv_make_table <- function(cond_inf_results, p, q, npar) {
  fe_names <-  rep(NA, p)
  for (i in 1:p) {
    if (i == 1) {
      fe_names[i] <- "$\\beta_0$"
    }else {
      fe_names[i] <- paste("$\\beta_", i, "$", sep = "")
    }
  }
  
  re_name <- re_names(q)
  row_names <-  c(rbind(c(fe_names, re_name), rep(" ", npar)))
  col_names <- c("ML(BFGS)", "ML(CG)", "bglmer[t]", "bglmer[n]", "MSPL")
  
  MAL_BFGS_col <-  c(t(getSummary(cond_inf_results$ML["BFGS", ])$coef[, 1:2]))
  MAL_CG_col <-  c(t(getSummary(cond_inf_results$ML["CG", ])$coef[, 1:2]))
  bglmer_t_col <-  c(t(getSummary(cond_inf_results$bglmer_t)$coef[, 1:2]))
  bglmer_n_col <-   c(t(getSummary(cond_inf_results$bglmer_n)$coef[, 1:2]))
  MSPAL_col <-   c(t(getSummary(cond_inf_results$MSPL)$coef[, 1:2]))
  
  table_mat <- cbind(MAL_BFGS_col, MAL_CG_col, bglmer_t_col, bglmer_n_col, MSPAL_col)
  table_mat <- make_math(table_mat, ndigits = 2, inline_math = FALSE)
  
  table_mat <- cbind(row_names, table_mat)
  table_mat <- gsub("NaN", "$-$", table_mat)
  colnames(table_mat) <- c(" ", col_names)
  xt <- xtable(table_mat,
               type = "latex",
               align = c("l", "l", rep("c", 5)))
  return(xt)
}