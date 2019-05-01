################################################################
##
## This resource is designed to introduce rstan
## as a tool for Bayesian data analysis
##
##################################################################


## Load rstan libraries
library(rstan)
library(rshinystan)
library("bayesplot")
library("ggplot2")
library(adegenet)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## 1 Simple binomial logistic regressions
## 1.1 Running the stan model

data1 <- read_rdump("Q:\\RProjects\\Teaching_Rstan\\data_IRS.R")

stan_base <- stan(file="Q:\\RProjects\\Teaching_Rstan\\simple_binomial_regression.stan", 
                  data=data1, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)

## extract the model parameters and predictive posterior draws (if run)
base <- extract(stan_base)

# Extract posterior draws for later use
posterior_cp <- as.array(stan_base)

#####################################
## 1.2 Diagnostics
## http://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
## display important HMC/NUTS diagnostic information

##Simplest way is to print the model
print(stan_base)

lp_cp <- log_posterior(stan_base)
np_cp <- nuts_params(stan_base)
## Check the chain acceptance stats
## This pulls out data for the iterations on:
  ## the acceptance_stat (somewhere bewteen 85% and 95% is good)
  ## the step_size (searching for the best parameters)
  ## The tree dephth (which we can alter to help the fit)
  ## The n_leapfrog
  ## The divergent
  ## energy
## We could also use to visualise these diagnostics
stan_diag(stan_base,
          information = c("sample","stepsize", "treedepth","divergence"),
          chain = 0)

## identifying collinearity between variables (which manifests as narrow bivariate plots) 
## as well as the presence of multiplicative non-identifiabilities (banana-like shapes)
color_scheme_set("darkgray") ##any issues will show up as red
mcmc_pairs(posterior_cp, np = np_cp, pars = c("alpha1","alpha2"))

## traceplots - time series plots of the chains converging (ideally!)
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, pars = "alpha1", np = np_cp) +
  xlab("Post-warmup iteration")
## we could also used
traceplot(stan_base, pars = c("alpha1","alpha2"),  inc_warmup = TRUE)
traceplot(stan_base, pars = c("alpha1","alpha2"),  inc_warmup = FALSE)

## overlaid histograms of the (centered) marginal energy distribution πE and the first-differenced distribution πΔE,
## ideally these look the same (Betancourt, 2017)
color_scheme_set("red")
mcmc_nuts_energy(np_cp)

## We could also pull out the histograms for parameters
stan_hist(stan_base, include = TRUE, unconstrain = FALSE,
          inc_warmup = FALSE, pars = c("alpha1","alpha2"))
## And visualise the Rhat stat (tideally 0.99 > Rhat < 1.01)
stan_rhat(stan_base,bins = 10) 
## or
rhats <- rhat(stan_base)
print(rhats)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

## Checking effective sample sizes
## The effective sample size is an estimate of the number of independent 
## draws from the posterior distribution of the estimand of interest
## The larger the ratio of neff to N the better (see Gelman et al. 2013, Stan Development Team 2018 for more details)
ratios_cp <- neff_ratio(stan_base)
print(ratios_cp)
mcmc_neff(ratios_cp, size = 2)
## only need worry about any neff/N less than 0.
## There are different algorithms for MCMC
## these ratios will depend not only on the model being fit but also on the particular MCMC algorithm used
## the one used by stan should give higher ratios than that used by Gibbs for example.

##Autocorrelation
mcmc_acf(posterior_cp, pars = c("alpha1","alpha2"), lags = 10)

#################################
## 1.3 Visuallising the model predictions against data
## For the binomial regression we can simply put in the equation and estimated parameters from the model runs
proportion_Y = data1$d_t/data1$n_t
x = data1$time

plot(proportion_Y ~ x, ylim = c(0,1), ylab = "Proportion of mosquitoes killed by IRS", xlab = "Time in days",
     cex.axis = 1.6,cex.lab = 1.6,cex = 1.6,pch=20,col="aquamarine3")
model_predictions_1 = array(dim=c(365,4000))
time = 1:365
for(i in 1:4000){
  model_predictions_1[,i] = 1/(1 + exp(-base$alpha2[i] * time - base$alpha1[i]))
  
}
for(i in 1:4000){
  lines(model_predictions_1[,i] ~ time,col=transp("grey",0.01))
}

#######################################################
##
## 2 Binomial logistic regressions with random effects
## 2.1 Running the stan model
##
#######################################################

data2 <- read_rdump("Q:\\RProjects\\Teaching_Rstan\\data_IRS_reff.R")

stan2 <- stan(file="Q:\\RProjects\\Teaching_Rstan\\random_effects_binomial_regression.stan", 
              data=data2, 
              warmup=500,
              control = list(adapt_delta = 0.9,
                             max_treedepth = 20),
              iter=1000, chains=4)

base2 = extract(stan2)

## Adjust the plot for the random effects of different studies
time=seq(1,365,length=365)
N_IRS = data2$N_IRS

mean_valssp_checker4 = mean_valssp_checker4u = mean_valssp_checker4l = 
  array(dim=c(length(time), N_IRS))

for(i in 1:N_IRS){
  mean_valssp_checker4[,i] = 1 / (1 + exp(-mean(base2$alpha1[,i]) - mean(base2$alpha2[,i])*time))
  mean_valssp_checker4u[,i] = 1 / (1 + exp(-quantile(base2$alpha1[,i],0.9) - quantile(base2$alpha2[,i],0.9)*time))
  mean_valssp_checker4l[,i] = 1 / (1 + exp(-quantile(base2$alpha1[,i],0.1) - quantile(base2$alpha2[,i],0.1)*time))
}

data2$prop_dead = data2$d_t/data2$n_t
colour_list = c("darkred","orange","aquamarine3","purple","blue","darkgreen")
for(i in 1:N_IRS){
  points(data2$prop_dead[data2$IRS == i]~data2$time[data2$IRS == i],col=colour_list[i],pch=14+i,cex=1.6)
  lines(mean_valssp_checker4[,i] ~ time,col=colour_list[i])
  polygon(c(time,rev(time)),c(mean_valssp_checker4u[,i],rev(mean_valssp_checker4l[,i])),border=FALSE,col=transp(colour_list[i],0.3))
}


## Or you can use the generated quantities
## eg for the 1st study
for(i in 1:2000){
  lines(base2$sp_ppc[i,1,] ~ time,col=transp(colour_list[1],0.01))
}


###########################################################
##
## 3.1 PKPD models example
##
###########################################################

# https://github.com/stan-dev/stancon_talks/tree/master/2017/Contributed-Talks/05_margossian/models/effCpt

RUNDATA <- read_rdump("Q:\\RProjects\\Teaching_Rstan\\template_data.R") ## Coverage data March 2019

## initial estimates will be generated randomly for each chain
init <- function(){
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       ke0 = exp(rnorm(1,log(1),0.2)),
       EC50 = exp(rnorm(1,log(100),0.2)),
       sigma = 0.5,
       sigmaResp = 20)
}

## Specify the variables for which you want history plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

nChains <- 4
nPost <- 1000 ## Number of post-warm-up samples per chain after thinning
nBurn <- 1000 ## Number of warm-up samples per chain after thinning
nThin <- 1
nIter <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = "Q:\\RProjects\\Teaching_Rstan\\Template_ODE_model.stan",
            data = RUNDATA,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))

stan_trace(fit, parametersToPlot)
print(fit)
