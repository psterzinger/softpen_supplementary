#Required packages 
library(lme4)
library(optimx) 
library(numDeriv)
library(blme)
library(xtable)
library(dotwhisker)
library(patchwork)
library(scales)
library(cowplot)
library(dplyr)
library(doMC)
d.binom <- read.table("http://pastebin.com/raw.php?i=0yDpEri8") 


# run glmer to get model matrices 
out_glmer <- glmer(response ~ type * p.validity * counterexamples + (1+counterexamples|code),data = d.binom, family = binomial)

# define data
data = list(Y=d.binom$response,X=getME(out_glmer,"X"),Z=cbind(rep(1,nrow(d.binom)),as.numeric(d.binom$counterexamples)-1),M=rep(1,nrow(d.binom)),grouping=as.numeric(d.binom$code)) 

p = ncol(data$X) 
q <- ncol(data$Z) 
npar <- p + 0.5*q*(q+1)

# run estimations 
MAL_BFGS <- get_MSPAL(start=rep(0,npar),data=data,nAGQ=1,mult=c(0,0),pen_log_sigma = NULL,optimization_methods ="BFGS",method_name = "MAL_BFGS")

MAL_CG <- get_MSPAL(start=rep(0,npar),data=data,nAGQ=1,mult=c(0,0),pen_log_sigma = NULL,optimization_methods ="CG",method_name = "MAL_CG")

MAL <- get_MSPAL(start=NULL,data=data,nAGQ=1,mult=c(0,0),pen_log_sigma = NULL,method_name = "MAL",control=list(all.methods = TRUE) )

bglmer_t <- get_bglmer(start=NULL,data=data,nAGQ=1,c_prior = "wishart",f_prior = "t")

bglmer_n <- get_bglmer(start=NULL,data=data,nAGQ=1,c_prior = "wishart",f_prior = "normal")

MSPAL <- get_MSPAL(start=rep(0,npar),data=data,nAGQ=1,mult=NULL,pen_log_sigma = function(logsigma){sum(nHuber(logsigma,1))})


# make table
fe_names <-  rep(NA,p) 
for(i in 1:p){
  if(i==1){
    fe_names[i] <-"$\\beta_0$"
  }else{
    fe_names[i] <- paste("$\\beta_",i,"$",sep="")
  }
}

re_names <- re_names(q) 
row_names <-  c(t(cbind(c(fe_names,re_names),rep("",npar))))
col_names <- c("MAL(BFGS)","MAL(CG)","bglmer(t)","bglmer(n)","MSPAL")

MAL_BFGS_col <-  c(t(cbind(MAL_BFGS$estimate,MAL_BFGS$se)))
MAL_CG_col <-  c(t(cbind(MAL_CG$estimate,MAL_CG$se)))
bglmer_t_col <-  c(t(cbind(bglmer_t$estimate,bglmer_t$se)))
bglmer_n_col <-  c(t(cbind(bglmer_n$estimate,bglmer_n$se)))
MSPAL_col <-  c(t(cbind(MSPAL$estimate,MSPAL$se)))

table_mat <- cbind(MAL_BFGS_col,MAL_CG_col,bglmer_t_col,bglmer_n_col,MSPAL_col)
table_mat <- make_math(table_mat,ndigits=2,inline_math=FALSE)

table_mat <- cbind(row_names,table_mat)
colnames(table_mat) <- c(" ",col_names ) 
xt <- xtable(table_mat,type="latex",align=c("l","l",rep("c",5)) ) 
print(xt, sanitize.text.function = function(x){x},include.rownames = FALSE)

# simulation 
MSPAL <- get_MSPAL(start=rep(0,npar),data=data,nAGQ=1,mult=NULL,pen_log_sigma = function(logsigma){sum(nHuber(logsigma,1))})
truth <- MSPAL$estimate 

fe_names <-  rep(NA,p) 
for(i in 1:p){
  if(i==1){
    fe_names[i] <- "beta[0]"
  }else{
    fe_names[i] <- paste("beta[",i,"]",sep="")
  }
}
re_names_p <- re_names_plot(q)
mathpar <- c(fe_names,re_names_p)
simul <-  perform_experiment(truth=truth,data=data,nsimu=1000,seed=0,mathpar = mathpar)
saveRDS(simul,"/Users/Philipp/Repositories/softpen/Data/Logic/cond_inf_sim1000.Rds")
simul_data <- readRDS("/Users/Philipp/Repositories/softpen/Data/Logic/cond_inf_sim1000.Rds")
sim <- simul_data
sim$mathpar <- c(fe_names,re_names(q)) 

percentile_table(sim,8,2) 

pdf("/Users/Philipp/Repositories/softpen/LaTeX/Paper/Figures/cond_inf_simul.pdf",width=10,height=5) 
bar_plot_experiment(simul_data,p=8,q=2)
dev.off() 

