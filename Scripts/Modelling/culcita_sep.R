library("lme4")
library("blme")
library("expint")
library("optimx")
library("numDeriv")
library("memisc")
source("../Software/MSPAL.R")

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
           "MSPAL" = m_MSPAL, digits = 2, #coef.style = "horizontal",
           signif.symbols = NULL, summary.stats = NULL)
}

## Data from the worked examples of Bolker (2015)
## downloaded from https://github.com/bbolker/mixedmodels-misc/raw/master/data/culcita.RData on 20220208
load("../../Data/culcita.RData")

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

baseline_both[["BFGS"]]$coef[1:4] - Cmat %*% baseline_none[["BFGS"]]$coef[1:4, 1]
baseline_both[["CG"]]$coef[1:4] - Cmat %*% baseline_none[["CG"]]$coef[1:4, 1]


baseline_both[["MSPAL"]]$coef[1:4] - Cmat %*% baseline_none[["MSPAL"]]$coef[1:4, 1]


baseline_both[["BGLMER(t)"]]$coef[1:4] - Cmat %*% baseline_none[["BGLMER(t)"]]$coef[1:4, 1]
baseline_both[["BGLMER(n)"]]$coef[1:4] - Cmat %*% baseline_none[["BGLMER(n)"]]$coef[1:4, 1]

baseline_none
baseline_both

## Full culcita data, "none" as baseline
dat_none_full <- list(
                X = model.matrix(~ ttt, data = culcita_none),
                Y = culcita_dat$predation,
                M = rep(1, nrow(culcita_dat)),
                grouping = culcita_dat$block)
baseline_none_full <- fit_all(rep(0, 5), dat_none_full)
