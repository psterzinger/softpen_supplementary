## Distributed as part of the supplementary material for the manuscript
## "Maximum softly-penalized likelihood for mixed effects logistic regression"
##
## Authors: Philipp Sterzinger, Ioannis Kosmidis
## Date: 1 June 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org"
       options(repos=r)})

pckg <- c("blme","lme4","doMC","parallel","numDeriv","patchwork","ggplot2","grid","xtable")

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

source(file.path(functions_path, "mv_MSPAL.R"))

sim_1_log <- file(file.path(results_path,"sim_1_log.txt"))
tryCatch({
ns <- c(50, 100, 200)
k <- 5
p <- 5
B <- 100
ln <- 11
lmin <- -10
lmax <- 10
lambdas <- seq(from = lmin, to = lmax, length.out = ln)
sigma <- 3

simul_list <- list()
for (i in seq_len(length(ns))) {
  set.seed(i)
  n <- ns[i]
  x1 <- rep(1, k * n)
  x2 <- rnorm(k * n)
  x3 <- rbinom(k * n, 1, .5)
  x4 <- rbinom(k * n, 1, .2)
  x5 <- rexp(k * n)
  X <- cbind(x1, x2, x3, x4, x5)
  colnames(X) <- c("intercept_FE", "norm", "bern_sep", "bern_noise", "exp")
  data <- list(X = X,
               Y = NULL,
               Z = rep(1, k * n),
               M = rep(1, k * n),
               grouping = rep(1:k, n))
  out_list <- NULL
  if (i == 1) {
    alt_start <-  rep(0, p + 1)
  }else{
    alt_start <- simul_list[[i - 1]] %>%
      group_by(method, parameter) %>%
      filter(lambda == lmin & method == "MSPAL") %>%
      summarise(mean_par = mean(estimate, na.rm = TRUE)) %>%
      dplyr::select(mean_par)
    alt_start <- rev(unlist(alt_start[, 2], use.names = FALSE))
  }
  for (j in 1:ln) {
    beta <- c(1, -0.5, lambdas[j], 0.25, -1)
    par <- c(beta, log(sigma))
    if (j > 1) {
      alt_start <- out %>%
        group_by(method, parameter) %>%
        summarise(mean_par = mean(estimate, na.rm = TRUE)) %>%
        filter(method == "MSPAL") %>%
        dplyr::select(mean_par)
      alt_start <- rev(unlist(alt_start[, 2], use.names = FALSE))
    }
    out <- mv_perform_experiment(truth = par,
                                 data = data,
                                 nsimu = B,
                                 nAGQ = 20,
                                 c_prior = "gamma",
                                 seed = i + j,
                                 alt_start = alt_start,
                                 optimization_methods = c("CG", "nlminb","nlm", "BFGS", "L-BFGS-B"), 
                                 ncores = 48)
    out$lambda <- lambdas[j]
    out_list <- rbind(out_list, out)
  }
  out_list$n <- ns[i]
  simul_list[[i]] <- out_list
}
simul <- do.call(rbind, simul_list)
saveRDS(simul, file.path(results_path,"Sim_1.Rds"))}, error = function(e) {
  writeLines(as.character(e), sim_1_log)
})
close(sim_1_log)


#simul_data <- readRDS(file.path(results_path, "Sim_1.Rds"))
simul_data <- simul 
## Figure
xlab <- expression(lambda)
ylab <- expression(hat(beta)[3])
truelab <- expression(beta[3])
ylims <- c(-20, 20)
sim <- simul_data
sim$method <- ordered(
  rep(c("ML", "bglmer[n]", "bglmer[t]", "MSPL"), each = 6),
  levels = c("MSPL", "bglmer[t]", "bglmer[n]", "ML"))
bp <-  mv_boxplot_simul(sim, truelab = truelab, par_ind = 3, ylims = ylims)

pdf(file.path(figures_path, "Fig4.pdf"), height = 8.27, width = 11.69)
bp
grid.draw(textGrob(ylab, x = .02, rot = 90))
grid.draw(textGrob(xlab, y = 0.01, rot = 0))
dev.off()

sink(file = file.path(results_path,"sim1_outliers.txt"))
mv_report_outliers(simul_data, ylims, 3)
sink()

## Table
sim1_tab <- mv_simul_table(sim)
addtorow <- list()
inds <- c(-1, -1, seq(from = 0, by = 3, length.out = 4), nrow(sim1_tab))
addtorow$pos <- as.list(inds)
addtorow$command <- c("\\toprule \n",
    "&& \\multicolumn{11}{c}{$\\lambda$} \\\\ \\cmidrule{3-13} \n",
    rep("\\cmidrule{3-13} \n", 4),
    "\\bottomrule \n")
print(sim1_tab,
  add.to.row = addtorow,
  include.colnames = TRUE,
  include.rownames = FALSE,
  sanitize.text.function = function(x) x,
  file = file.path(results_path, "Sim_1_tab.tex"),
  hline.after = NULL)
