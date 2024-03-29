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

pckg <- c("lme4",
	"optimx",
	"numDeriv",
	"blme",
	"xtable",
	"dotwhisker",
	"patchwork",
	"scales",
	"cowplot",
	"dplyr",
	"doMC",
	"memisc",
	"ggplot2",
	"parallel") 

for (i in seq_len(length(pckg))) {
  if (!is.element(pckg[i], installed.packages()[, 1]))
    install.packages(pckg[i], dep = TRUE)
  require(pckg[i], character.only = TRUE)
}
update.packages(ask = FALSE) 

ncores <- as.numeric(Sys.getenv("NCORES"))

functions_path <- "./Functions"
data_path <- "../Data"
results_path <- "../Results"
figures_path <- results_path

source(file.path(functions_path, "mv_MSPAL.R"))

# Conditional inference data from "http://pastebin.com/raw.php?i=0yDpEri8"), last accessed: May 31, 2022
d.binom <- read.table(file.path(data_path, "d.binom.csv"),
  header = TRUE,
  sep = ",")

data <- list(Y = d.binom$response,
            X = model.matrix(
              ~ type * p.validity * counterexamples, data = d.binom),
            Z = cbind(
              rep(1, nrow(d.binom)),
              as.numeric(factor(d.binom$counterexamples)) - 1),
            M = rep(1, nrow(d.binom)),
            grouping = as.numeric(factor(d.binom$code)))
p <- ncol(data$X)
q <- max(1, ncol(data$Z))
npar <- p + 0.5 * q * (q + 1)
cond_inf_results <- mv_fit_all(start = rep(0, npar), data = data)

# make table
tab <- mv_make_table(cond_inf_results, p, q, npar)
print(tab,
    sanitize.text.function = function(x) x,
    include.rownames = FALSE,
    file = file.path(results_path,"cond_inf_table.tex"),
    compress = FALSE,
    only.contents = TRUE,
    comment = FALSE)

# simulation
truth <- cond_inf_results$MSPL$estimate
fe_names <-  rep(NA, p)
for (i in 1:p) {
  if (i == 1) {
    fe_names[i] <- "beta[0]"
  }else {
    fe_names[i] <- paste("beta[", i, "]", sep = "")
  }
}
re_names_p <- re_names_plot(q)
mathpar <- c(fe_names, re_names_p)

sim_cond_inf_log <- file(file.path(results_path,"sim_cond_inf_log.txt"))
tryCatch({
simul_data <-  mv_perform_experiment(truth = truth,
                            data = data,
                            nsimu = 2,
                            seed = 0,
                            mathpar = mathpar, 
                            ncores = ncores)

saveRDS(simul_data, file.path(results_path, "cond_inf_sim10000.Rds"))
}, error = function(e) {
  writeLines(as.character(e), sim_cond_inf_log)
})
close(sim_cond_inf_log)

## simulation analysis
#simul_data <- readRDS(file.path(results_path, "cond_inf_sim10000.Rds"))
fe_names <-  rep(NA, p)
    for (i in 1:p) {
        if (i == 1) {
            fe_names[i] <- "$\\beta_0$"
        }else {

            fe_names[i] <- paste("$\\beta_", i, "$", sep = "")
        }
    }
sim <- simul_data
sim$method <- ordered(
  rep(c("ML", "bglmer[n]", "bglmer[t]", "MSPL"), each = npar),
  levels = c("MSPL", "bglmer[t]", "bglmer[n]", "ML"))
re_name <- re_names(q)
sim$mathpar <- c(fe_names, re_name)

sink(file=file.path(results_path,"perc_table.tex"))
mv_percentile_table(sim, 8, 2)
sink()

simul_data$method <- ordered(
  rep(c("ML", "bglmer[n]", "bglmer[t]", "MSPL"), each = npar),
  levels = c("MSPL", "bglmer[t]", "bglmer[n]", "ML"))

pdf(file.path(figures_path, "Fig2.pdf"), width = 10, height = 5)
  mv_bar_plot_experiment(simul_data, p = 8, q = 2)
dev.off()

### Full data
pdf(file.path(figures_path, "Fig3.pdf"), width = 10, height = 5)
  mv_bar_plot_experiment(simul_data, p = 8, q = 2, var_only = FALSE)
dev.off()
sink(file=file.path(results_path,"perc_table_full.tex"))
mv_percentile_table(sim, p = 8, q = 2, var_only = FALSE)
sink()

if (FALSE) {
  fe_names <-  rep(NA, p)
    for (i in 1:p) {
        if (i == 1) {
            fe_names[i] <- "$\\beta_0$"
        }else {
            fe_names[i] <- paste("$\\beta_", i, "$", sep = "")
        }
    }
re_name <- re_names(q)
simul_data$mathpar <- c(fe_names, re_name)
sums <- mv_summarize_experiment(simul_data)
methods <- as.character(unique(simul_data$method))
test4 <- sums %>%
    group_by(method) %>%
    arrange(factor(parameter, levels = rev(levels(parameter)))) %>%
    dplyr::select(method, mathpar, bias, var, mse, pu, nused) %>%
    rename(Bias = bias, Variance = var, MSE = mse, PU = pu, R = nused)

varnames <- c("Bias", "Variance", "MSE", "PU", "R")
minitable2 <- function(var_name) {
    test4 %>%
    dplyr::select(method, mathpar, var_name) %>%
    pivot_wider(names_from = mathpar,
        values_from = var_name)
}
test4 <- do.call(rbind, lapply(varnames, minitable2))

var_col <- do.call(rbind, lapply(varnames,
    function(var_name) {
        out <- matrix(NA, 4, 1);
        out[3] <- var_name;
        return(out)}))
test4$var_col <- var_col[, 1]
test4 <- test4 %>% relocate(var_col, method, .before = "$\\beta_0$")
colnames(test4)[1:2] <- " "

addtorow <- list()
inds <- c(-1, seq(from = 0, by = 4, length.out = length(varnames)), nrow(test4))
addtorow$pos <- as.list(inds)
addtorow$command <- c("\\toprule \n",
    rep("\\cmidrule{3-13} \n", length(varnames)),
    "\\bottomrule \n")

print(xtable(test4, digits = 2),
    add.to.row = addtorow,
    include.colnames = TRUE,
    include.rownames = FALSE,
    sanitize.text.function = function(x) x,
    #file = file.path(results_path, "cond_inf_table_full.tex"),
    hline.after = NULL,
    scalebox = 0.7)
}

seshinfo <- file(file.path(results_path,"seshinfo.tex"))
writeLines(toLatex(sessionInfo()), seshinfo)
close(seshinfo)
