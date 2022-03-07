library(blme)
library(lme4)
library(doMC)
library(parallel)
library(numDeriv)

ns <- c(50,100,200)
k <- 5
p <- 5
B <- 100
sn <- 8
smin <- -5
smax <- 2
lambdas <- seq(from=smin,to=smax,length.out = sn)
beta <- c(1,-0.5,0.5,0.25,-1)

simul_list <- list()
for(i in 1:3){
  set.seed(i)
  n <- ns[i]
  x1 <- rep(1,k*n) 
  x2 <- rnorm(k*n) 
  x3 <- rbinom(k*n,1,.5) 
  x4 <- rbinom(k*n,1,.2) 
  x5 <- rexp(k*n)
  X <- cbind(x1,x2,x3,x4,x5)
  colnames(X) <- c("intercept_FE","norm","bern_sep","bern_noise","exp")
  clusters <- rep(1:k,n)
  z <- rep(1,k*n)
  data <- list(X=X,Y=NULL,Z=rep(1,k*n),M=rep(1,k*n),grouping=clusters)
  out_list <- NULL
  if(i==1){
    alt_start <-  rep(0, p+1)
  }else{
    alt_start <- simul_list[[i-1]] %>%
      group_by(method,parameter) %>%
      filter(lambda==smin & method=="MSPAL") %>%
      summarise(mean_par = mean(estimate,na.rm=TRUE)) %>%
      dplyr::select(mean_par)
    alt_start <- rev(unlist(alt_start[,2],use.names = FALSE)  )
  }
  for(j in 1:sn){
    print(c(i,j)) 
    logsigma <- lambdas[j]
    par <- c(beta,logsigma)
    if(j>1){
      alt_start <- out %>%
        group_by(method,parameter) %>%
        summarise(mean_par = mean(estimate,na.rm=TRUE)) %>%
        filter(method=="MSPAL") %>%
        dplyr::select(mean_par)
      alt_start <- rev(unlist(alt_start[,2],use.names = FALSE)  )
    }
    out <- perform_experiment(truth=par,data=data,nsimu=B,nAGQ=20,c_prior="gamma",seed = i+j,alt_start = alt_start,optimization_methods = c("CG","nlminb","BFGS","L-BFGS-B"))
    out$lambda <- lambdas[j]
    out_list <- rbind(out_list,out)
  }
  out_list$n <- ns[i]
  simul_list[[i]] <- out_list
}
simul <- do.call(rbind,simul_list)
saveRDS(simul,"/Users/Philipp/Repositories/softpen/Data/Simulation_2.Rds")

library(patchwork)
library(ggplot2) 
library(grid)
library(xtable)

simul_data <- readRDS("/Users/Philipp/Repositories/softpen/Data/Simulation_2.Rds")

## Plot
xlab = expression(lambda)
ylab = expression(widehat(log(sigma))) 
truelab <- expression(log(sigma))
bp <-  boxplot_simul(simul_data,truelab = truelab,par_ind = 6,ylim=c(-10,10))

pdf("/Users/Philipp/Repositories/softpen/LaTeX/Paper/Figures/sim2.pdf",height=8.27,width=11.69) 
bp 
grid.draw(textGrob(ylab, x = .02, rot = 90)) 
grid.draw(textGrob(xlab, y = 0.01, rot = 0)) 
dev.off() 

report_outliers(simul_data,ylims,6)
## Table 
simul_table(simul_data,log_mask = c(rep(FALSE,5),TRUE),par_ind = 6)
