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

pckg <- c("lme4","blme","expint","optimx","numDeriv","memisc")

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

## A wrapper function to fit all methods in Example 1 of the main text and pretty-report the output
fit_all <- function(start, data, nAGQ = 100, cv = 1) {
    m_MAL <- glmm_fit(start, mult = c(0, 0), data = data, nAGQ = nAGQ,
                   method = c("BFGS", "CG"))
    m_bglmer_t <- get_bglmer(start, data = data, nAGQ = nAGQ,
                              f_prior = bquote(t), c_prior = bquote(wishart))
    m_bglmer_n <- get_bglmer(start, data = data, nAGQ = nAGQ,
                             f_prior = bquote(normal), c_prior = bquote(wishart))
    m_MSPAL <- get_MSPAL(start, data = data, nAGQ = nAGQ,
                          pen_log_sigma = function(x) nHuber(x, cv))
    mtable("BFGS" = m_MAL["BFGS", ],  "CG" = m_MAL["CG", ],
           "BGLMER(t)" = m_bglmer_t, "BGLMER(n)" = m_bglmer_n,
           "MSPL" = m_MSPAL, digits = 2, #coef.style = "horizontal",
           signif.symbols = NULL, summary.stats = NULL)
}

## Data from the worked examples of Bolker (2015)
## downloaded from https://github.com/bbolker/mixedmodels-misc/raw/master/data/culcita.RData on 20220208
load(file.path(data_path, "culcita.RData"))

atypical_obs <- 20
## Remove atypical observation and fit using treatment contrasts with "none" as baseline
## culcita_dat[extreme_obs,]

## "none" as baseline
culcita_none <- within(culcita_dat, ttt <- relevel(ttt, ref = "none"))
dat_none <- list(X = model.matrix(~ ttt, data = culcita_none)[-atypical_obs, ],
                 Y = culcita_dat$predation[-atypical_obs],
                 M = rep(1, nrow(culcita_dat))[-atypical_obs],
                 grouping = culcita_dat$block[-atypical_obs])
baseline_none <- fit_all(rep(0, 5), dat_none)

## "both" as baseline
culcita_both <- within(culcita_dat, ttt <- relevel(ttt, ref = "both"))
dat_both <- within(dat_none, X <- model.matrix(~ ttt, data = culcita_both)[-atypical_obs, ])
baseline_both <- fit_all(rep(0, 5), dat_both)

## Differences in estimates
Cmat <- matrix(c(1, 0, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1,
                 1, -1, -1, -1), 4, 4)
sink(file=file.path(results_path,"estimate_diffs.txt"))
print("Differences in estimates") 
print("ML BFGS:")
baseline_both[["BFGS"]]$coef[1:4] - Cmat %*% baseline_none[["BFGS"]]$coef[1:4, 1]
print("ML CG:")
baseline_both[["CG"]]$coef[1:4] - Cmat %*% baseline_none[["CG"]]$coef[1:4, 1]

print("MSPL:")
baseline_both[["MSPL"]]$coef[1:4] - Cmat %*% baseline_none[["MSPL"]]$coef[1:4, 1]

print("BGLMER(t):")
baseline_both[["BGLMER(t)"]]$coef[1:4] - Cmat %*% baseline_none[["BGLMER(t)"]]$coef[1:4, 1]
print("BGLMER(n):")
baseline_both[["BGLMER(n)"]]$coef[1:4] - Cmat %*% baseline_none[["BGLMER(n)"]]$coef[1:4, 1]
sink() 


writeLines(toLatex(baseline_none), con = file.path(results_path,"baseline_none.tex"))
writeLines(toLatex(baseline_both), con = file.path(results_path,"baseline_both.tex"))

## Full culcita data, "none" as baseline
dat_none_full <- list(
                X = model.matrix(~ ttt, data = culcita_none),
                Y = culcita_dat$predation,
                M = rep(1, nrow(culcita_dat)),
                grouping = culcita_dat$block)
baseline_none_full <- fit_all(rep(0, 5), dat_none_full)
writeLines(toLatex(baseline_none_full), con = file.path(results_path,"baseline_none_full.tex"))
