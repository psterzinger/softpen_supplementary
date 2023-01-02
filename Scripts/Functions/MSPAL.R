## Distributed as part of the supplementary material for the manuscript
## "Maximum softly-penalized likelihood for mixed effects logistic regression"
##
## Authors: Philipp Sterzinger, Ioannis Kosmidis
## Date: 1 June 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!
infer_loglik <- function(data, nAGQ = 100) {
    data$quadtable <- NULL
    p <- ncol(data$X)
    dfun <- glmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                  family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
    ## Change to logsigma parameterization and fixef first
    function(par) -dfun(c(exp(par[p + 1]), par[1:p])) / 2
}


loglik <- function(par, data, nAGQ = 100, lfun = NULL) {
    if (is.null(lfun))
        lfun <- infer_loglik(data, nAGQ)
    out <- try(lfun(par), silent = TRUE)
    if (is(out, "try-error")) NA else out
}

## Numerically evaluate score function
scores <- function(par, data, nAGQ, lfun) {
    grad(loglik, par, data = data, nAGQ = nAGQ, lfun = lfun)
}

## Options for pen_log_sigm
linear_LR <- function(x, c, d, e, centering = 0) {
    x <- x - centering
    ifelse(x <= c, - e * c * x + e * c^2 / 2,
           ifelse(x >= d, - e * d * x + e * d^2 / 2,
                  - e * x^2 / 2))
}

nHuber <- function(x, a) {
    ifelse(abs(x) > a,
           - a * (abs(x) - a / 2),
           - x^2 / 2)
}

chungetal <- function(x) {
    1.5 * x
}

## Overall penalty
penalty <- function(par, data, mult = c(0, 0), pen_log_sigma = NULL) {
    X <- data$X
    Y <- data$Y
    M <- data$M
    grouping <- data$grouping
    npar <- length(par)
    p <- ncol(X)
    beta <- par[seq.int(p)]
    logsigma <- par[p + 1]
    psi <- exp(logsigma)
    fetas <- drop(X%*%beta)
    probs <- plogis(fetas)
    w <- M * probs * (1 - probs)
    pen <- mult[1] * log(det(crossprod(X * sqrt(w)))) / 2 # a multiple of jeffreys
    if (!is.null(pen_log_sigma)) {
        pen <- pen + mult[2] * pen_log_sigma(logsigma)
    }
    pen
}

glmm_fit <- function(start, data,
                     mult = NULL,
                     nAGQ = 100,
                     pen_log_sigma = NULL, ...) {
    if (is.null(mult)) {
        mult <- rep(multiplier(data), 2)
    }
    penalize <- all(mult == 0)
    lfun <- infer_loglik(data, nAGQ = nAGQ)
    obj <- function(pars) {
        pen <- penalty(pars, data,
                       mult = mult,
                       pen_log_sigma = pen_log_sigma)
        -loglik(pars, data, nAGQ = nAGQ, lfun = lfun) - pen
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


coef.glmm_fit <- function(object, best = TRUE) {
    class(object) <- class(object)[-1]
    cc <- coef(object)
    if (best) drop(cc[1, ]) else  cc
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

sqrtdiaginv <- function(mat) {
    out <- try(sqrt(diag(solve(mat))), silent = TRUE)
    if (inherits(out, "try-error")) {
        out <- rep(NA, nrow(mat))
    }
    out
}

boundary_diagnostics <- function(object, best = TRUE) {
    if (best) {
        sqrtdiaginv(hess(object, best = TRUE))
    }
    else {
        lapply(hess(object, best = FALSE), sqrtdiaginv)
    }
}


info <- function(object, best = TRUE) {
    data <- attr(object, "data")
    nAGQ <- attr(object, "nAGQ")
    lfun <- attr(object, "lfun")
    info_fun <- function(x) {
        -hessian(loglik, x, data = data, nAGQ = nAGQ, lfun = lfun)
    }
    if (best) {
        info_fun(coef(object, best = TRUE))
    }
    else {
        apply(coef(object, best = FALSE), 1, info_fun, simplify = FALSE)
    }
}

se <- function(object, best = TRUE) {
    if (best) {
        sqrtdiaginv(info(object, best = TRUE))
    }
    else {
        do.call("rbind", lapply(info(object, best = FALSE), sqrtdiaginv))
    }
}

multiplier <- function(data) {
    vc_eta <- with(data, tcrossprod(X %*% solve(crossprod(X / 2)), X))
    mean(diag(vc_eta))^{0.5}
}


format_dec <- function(x, k) format(round(x, k), trim = TRUE, nsmall = k)

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

## Functions for simulation study
get_MSPAL <- function(start, data, nAGQ = 100,
                      mult = NULL, pen_log_sigma = NULL,
                      optimization_methods = c("BFGS", "CG"),
                      method_suffix = "", method_name = NULL) {
    pp <- ncol(data$X) + 1
    out <- try({
        fit0 <- glmm_fit(start,
                         data = data, nAGQ = nAGQ,
                         method = optimization_methods,
                         mult = mult,
                         pen_log_sigma = pen_log_sigma)
        mult <- attr(fit0, "mult")
        cfe <- coef(fit0, best = TRUE)
        sfe <- se(fit0, best = TRUE)
        grad <- gradient(fit0, best = TRUE)
        bd <- boundary_diagnostics(fit0, best = TRUE)
    }, silent = TRUE)
    if (is(out, "try-error")) {
        grad <- bd <- cfe <- sfe <- rep(NA, pp)
    }
    if (is.null(method_name)) {
        method_name <- paste0("MSPAL_", method_suffix, "[", paste(format_dec(mult, 1), collapse =","), "]")
    }
    out <- data.frame(estimate = cfe,
                      se = sfe,
                      grad = grad,
                      bd = bd,
                      method = method_name,
                      parameter = factor(c(colnames(data$X), "log_sigma"), levels = c(colnames(data$X), "log_sigma"), ordered = TRUE))
    attr(out, "mult") <- mult
    class(out) <- c("get_results_output", class(out))
    attr(out, "call") <- match.call()
    out
}

## starting values only for theta
get_glmer <- function(start = NULL, data, nAGQ = 100, method_name = NULL) {
    pp <- ncol(data$X) + 1
    data$quadtable <- NULL
    out <- try({
        if (is.null(start)) {
            mod <- glmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                         family = binomial(), nAGQ = nAGQ)
        }
        else {
            start <- list(fixef = start[-pp], theta = start[pp])
            mod <- glmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                         family = binomial(), nAGQ = nAGQ, start = start)
        }
        cfe <- c(mod@beta, log(mod@theta))
        bd_diagnostics <- transform_bd_grad(mod,pp,1) 
        bd <- bd_diagnostics$bd
        grad <- bd_diagnostics$grad
        ## Get deviance for unpenalized model to compute standard errors
        dfun <- glmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                      family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
        ## Change to logsigma parameterization and fixef first
        nll <- function(pars) dfun(c(exp(pars[pp]), pars[-pp])) / 2
        sfe <- sqrtdiaginv(hessian(nll, cfe))
    }, silent = TRUE)
    if (is(out, "try-error")) {
        grad <- bd <- cfe <- sfe <- rep(NA, pp)
    }
    out <- data.frame(estimate = cfe,
                      se = sfe,
                      grad = grad,
                      bd = bd,
                      method = ifelse(is.null(method_name), "glmer", method_name),
                      parameter = factor(c(colnames(data$X), "log_sigma"), levels = c(colnames(data$X), "log_sigma"), ordered = TRUE))
    class(out) <- c("get_results_output", class(out))
    attr(out, "call") <- match.call()
    out
}


## starting values only for theta
get_bglmer <- function(start = NULL, data, nAGQ = 20,
                       c_prior, f_prior, method_name = NULL) {
    pp <- ncol(data$X) + 1
    data$quadtable <- NULL
    out <- try({
        if (is.null(start)) {
            mod <- bglmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                          family = binomial, nAGQ = nAGQ,
                          cov.prior = c_prior, fixef.prior = f_prior)
        }
        else {
            start <- list(fixef = start[-pp], theta = start[pp])
            mod <- bglmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                          family = binomial, nAGQ = nAGQ, start = start,
                          cov.prior = c_prior, fixef.prior = f_prior)
        }
        cfe <- c(mod@beta, log(mod@theta))
        bd_diagnostics <- transform_bd_grad(mod,pp,1) 
        bd <- bd_diagnostics$bd
        grad <- bd_diagnostics$grad
        ## Get deviance for unpenalized model to compute standard errors
        dfun <- glmer(Y ~ -1 + X + (1 | grouping), weights = M, data = data,
                      family = binomial, nAGQ = nAGQ, devFunOnly = TRUE)
        ## Change to logsigma parameterization and fixef first
        nll <- function(pars) dfun(c(exp(pars[pp]), pars[-pp])) / 2
        sfe <- sqrtdiaginv(hessian(nll, cfe))
    }, silent = TRUE)
    if (is(out, "try-error")) {
        grad <- bd <- cfe <- sfe <- rep(NA, pp)
    }
    out <- data.frame(estimate = cfe,
                      se = sfe,
                      grad = grad,
                      bd = bd,
                      method = ifelse(is.null(method_name), paste0("BGLMER[", as.character(f_prior), ",", as.character(c_prior), "]"), method_name),
                      parameter = factor(c(colnames(data$X), "log_sigma"), levels = c(colnames(data$X), "log_sigma"), ordered = TRUE))
    class(out) <- c("get_results_output", class(out))
    attr(out, "call") <- match.call()
    out
}


## Helper for progress reporting in simulations
report_res <- function(out, digits = 2) {
    est <- format_dec(out$estimate, digits)
    bd <- format_dec(out$bd, digits)
    mm <- unique(out$method)
    nc <- nchar(mm)
    cat(paste0(mm, ":"), "\t\t")
    cat(paste0(est, " [", bd, "] "), "\n")
}

## getSummary helper functions to interface memisc::mtable
getSummary.get_results_output <- function(x) {
    coefs <- x[, "estimate"]
    ses <- x[, "se"]
    stat <- coefs/ses
    p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    cc <- cbind(coefs, ses, stat, p)
    colnames(cc) <- c("est", "se", "stat", "p")
    rownames(cc) <- x[, "parameter"]
    list(coef = cc, call = attr(x, "call"))
}

getSummary.glmm_fit <- function(x) {
    coefs <- coef(x)
    ses <- se(x)
    stat <- coefs/ses
    p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    cc <- cbind(coefs, ses, stat, p)
    colnames(cc) <- c("est", "se", "stat", "p")
    rownames(cc) <- c(colnames(attr(x, "data")$X), "log_sigma")
    list(coef = cc, call = attr(x, "call"))
}


## Functions for simulation experiments

simufun <- function(par, data, nsimu, seed = NULL) {
    X <- data$X
    Y <- data$Y
    M <- data$M
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
    q <- nrow(X)
    ## Number of groups
    ng <- nlevels(grouping)
    beta <- par[seq.int(p)]
    psi <- exp(par[p + 1])
    fetas <- drop(X%*%beta)
    reff <- matrix(rnorm(ng*nsimu, 0, psi), nrow = ng)
    etas <- fetas + reff[as.numeric(grouping),]
    probs <- plogis(etas)
    val <- apply(probs, 2, function(prob) rbinom(q, M, prob))
    val
}

fit_all_methods <- function(startl, data, s, trace = TRUE, nAGQ = 100, optimization_methods) {
    for (st in startl) {
        res_mal <- get_MSPAL(start = st, data = data, nAGQ = nAGQ,
                             mult = c(0, 0),
                             optimization_methods = optimization_methods,
                             method_name = "MAL")
        if (isTRUE(all(abs(res_mal$grad) < 1e-03))) break
    }
    startl <- c(list(res_mal$estimate), startl)
    for (st in startl) {
        res_mspal <- get_MSPAL(start = st, data = data, nAGQ = nAGQ,
                               pen_log_sigma = function(x) nHuber(x, 1.0),
                               optimization_methods = optimization_methods,
                               method_name = "MSPAL")
        if (isTRUE(all(abs(res_mspal$grad) < 1e-03))) break
    }
    startl <- c(list(res_mspal$estimate), startl)
    for (st in startl) {
        res_bglmer_t <- get_bglmer(start = st, data = data,  nAGQ = nAGQ,
                                   c_prior = bquote(wishart), f_prior = bquote(t),
                                   method_name = "bglmer(t)")
        if (isTRUE(all(abs(res_bglmer_t$grad) < 1e-03))) break
    }
    startl <- c(list(res_bglmer_t$estimate), startl)
    for (st in startl) {
        res_bglmer_n <- get_bglmer(start = res_bglmer_t$estimate, data = data, nAGQ = nAGQ,
                                   c_prior = bquote(wishart), f_prior = bquote(normal),
                                   method_name = "bglmer(n)")
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

perform_experiment <- function(truth, data,
                               nsimu = 100, seed = NULL,
                               alt_start = rep(0, length(truth)), nAGQ = 100,
                               optimization_methods = c("L-BFGS-B", "nlminb"),
                               ncores = 10,
                               mathpar = NULL) {
    registerDoMC(10)
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
    simulated_data <- simufun(truth, data = data, nsimu = nsimu)
    base_startl <- list(truth, alt_start)
    out <- mclapply(1:nsimu, function(j) {
        res <- fit_all_methods(startl = base_startl,
                               data = within(dat, Y <- simulated_data[, j]),
                               s = j,
                               trace = TRUE,
                               optimization_methods = optimization_methods,
                               nAGQ = nAGQ)
        res$sample <- j
        res
    }, mc.cores = ncores)
    out <- do.call("rbind", out)
    out$true <- truth
    parnames <- c(colnames(data$X), "log_sigma")
    if (is.null(mathpar)) {
        out$mathpar <- out$parameter
    }
    else {
        names(mathpar) <- parnames
        out$mathpar <-  recode(out$parameter, !!!mathpar)
    }
    out$mathpar <- factor(out$mathpar, levels = rev(mathpar),
                          ordered = TRUE)
    out$parameter <- factor(out$parameter, levels = rev(parnames),
                            ordered = TRUE)
    out$method <- factor(out$method, levels = rev(c("MAL", "bglmer(n)", "bglmer(t)", "MSPAL")),
                         ordered = TRUE)
    out
}

## parmath should correspond to the order of c(colnames(data$X), "log_sigma")
summarize_experiment <- function(out, se_threshold = 40, grad_threshold = 1e-03,
                                 estimate_threshold = 50) {
    ## Exlusion
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

plot_experiment <- function(out, out_summary = NULL, xlims = c(-10, 10),
                            max_var = 10, max_abs_bias = 3, ...) {
    if (is.null(out_summary)) out_summary <- summarize_experiment(out, ...)
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
            scale_color_grey(start = 0.1, end = 0.8) + scale_fill_grey(start = 0.2, end = 0.9)
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
                scale_x_continuous(breaks = ntried/5 * (0:5),
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
        scale_color_grey(start = 0.1, end = 0.8) +
        scale_fill_grey(start = 0.2, end = 0.9)
    ## patchwork
    ((bpp + plots[[1]] + plot_layout(widths = c(4, 1))) /
            (Reduce("+", plots[-1]) + plot_layout(ncol = length(stats) - 1)))
}
