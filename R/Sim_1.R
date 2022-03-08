library(blme)
library(lme4)
library(doMC)
library(parallel)
library(numDeriv)

ns <- c(50,100,200)
k <- 5
p <- 5
B <- 100
ln <- 11
lmin <- -10
lmax <- 10
lambdas <- seq(from=lmin,to=lmax,length.out = ln)
sigma <- 3

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
  sigma <- 3
  data <- list(X=X,Y=NULL,Z=rep(1,k*n),M=rep(1,k*n),grouping=clusters)
  out_list <- NULL
  if(i==1){
    alt_start <-  rep(0, p+1)
  }else{
    alt_start <- simul_list[[i-1]] %>% 
      group_by(method,parameter) %>% 
        filter(lambda==lmin & method=="MSPAL") %>%
          summarise(mean_par = mean(estimate,na.rm=TRUE)) %>%
            dplyr::select(mean_par)
    alt_start <- rev(unlist(alt_start[,2],use.names = FALSE)  )
  }
  for(j in 1:ln){
    print(c(i,j)) 
    beta <- c(1,-0.5,lambdas[j],0.25,-1)
    par <- c(beta,log(sigma))
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
simul_list <- do.call(rbind,simul_list)

#saveRDS(simul_list,"../Data/Simulation_1.Rds")
saveRDS(simul_list,"/Users/Philipp/Repositories/softpen/Data/Simulation_1.Rds")

library(patchwork)
library(ggplot2) 
library(grid)
library(xtable)

#simul_data <- readRDS("../Data/Simulation_1.Rds")
simul_data <- readRDS("/Users/Philipp/Repositories/softpen/Data/Simulation_1.Rds")

xlab = expression(lambda)
ylab = expression(hat(beta)[3]) 
truelab <- expression(beta[3])
ylims=c(-20,20)

bp <-  boxplot_simul(simul_data,truelab = truelab,par_ind = 3,ylims=ylims)

#pdf("../LaTeX/Figures/sim1.pdf",height=8.27,width=11.69) 
pdf("/Users/Philipp/Repositories/softpen/LaTeX/Paper/Figures/sim1.pdf",height=8.27,width=11.69) 
bp 
grid.draw(textGrob(ylab, x = .02, rot = 90)) 
grid.draw(textGrob(xlab, y = 0.01, rot = 0)) 
dev.off() 

report_outliers(simul_data,ylims,3)
## Table 
simul_table(simul_data)
