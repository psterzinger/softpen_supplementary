## Distributed as part of the supplementary material for the manuscript
## "Maximum softly-penalized likelihood for mixed effects logistic regression"
##
## Authors: Philipp Sterzinger, Ioannis Kosmidis
## Date: 2 January 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org"
       options(repos=r)})

pckg <- c("numDeriv", "lme4", "blme", "optimx", "dplyr", "parallel", "doMC", "ggplot2", "patchwork", "scales") 

for (i in seq_len(length(pckg))) {
  if (!is.element(pckg[i], installed.packages()[, 1]))
    install.packages(pckg[i], dep = TRUE)
  require(pckg[i], character.only = TRUE)
}
update.packages(ask = FALSE) 

functions_path <- "./Functions"
data_path <- "../Data"
results_path <- "../Results"
figures_path <- results_path

source(file.path(functions_path, "MSPAL.R"))
load(file.path(data_path, "culcita.RData"))

## Parameterization
estimate_threshold <- 30
se_threshold <- 30
grad_threshold <- 0.001
nAGQ <- 100

## Prepare list with model matrix, response vector, weights, and grouping variable
dat <- list(X = model.matrix(~ ttt, data = culcita_dat),
            Y = culcita_dat$predation,
            M = rep(1, nrow(culcita_dat)),
            grouping = culcita_dat$block)

## Get MAL
fit <- glmm_fit(c(0, 0, 0, 0, 0), mult = c(0, 0), data = dat,
                method = c("L-BFGS-B", "nlminb"), nAGQ = nAGQ)

truth <- drop(coef(fit))

sim_culcita_log <- file(file.path(results_path,"sim_culcita_log.txt"))
tryCatch({
out <- perform_experiment(truth = truth,
                          data = dat,
                          nsimu = 10000,
                          seed = 123,
                          alt_start = rep(0, length(truth)),
                          nAGQ = nAGQ,
                          optimization_methods = c("L-BFGS-B", "nlminb"),
                          ncores = 48,
                          mathpar = c("beta[0]",
                          "beta[1]",
                          "beta[2]",
                          "beta[3]",
                          "log(sigma)"))



out_summary <- summarize_experiment(out,
                                    estimate_threshold = estimate_threshold,
                                    se_threshold = se_threshold,
                                    grad_threshold = grad_threshold)

#load(file.path(results_path, "simulation_study_example1.rda"))
out$method <- ordered(
  rep(c("ML", "bglmer[n]", "bglmer[t]", "MSPL"), each = 5),
  levels = c("MSPL", "bglmer[t]", "bglmer[n]", "ML"))
pdf(file.path(figures_path,"Fig1.pdf"), width = 8, height = 7)
plot_experiment(out,
                estimate_threshold = estimate_threshold,
                se_threshold = se_threshold,
                grad_threshold = grad_threshold,
                xlims = c(-11, 11))
dev.off()

out %>%
        group_by(method, sample) %>%
        summarize(failed = any(is.na(estimate) |
        is.na(se) | isTRUE(abs(estimate) > estimate_threshold) |
        isTRUE(se > se_threshold) |
        isTRUE(abs(grad) > grad_threshold))) %>%
            summarize(mean(failed))


save.image(file.path(results_path, "simulation_study_example1.rda"))

}, error = function(e) {
  writeLines(as.character(e), sim_culcita_log)
})
close(sim_culcita_log)
