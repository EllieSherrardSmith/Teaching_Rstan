// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;                  // rows of data
  
  int<lower=0> n_t[N];             // Total number of mosquitoes entering IRS huts
  int<lower=0> d_t[N];             // Number mosquites dead sprayed hut
  
  vector<lower=0>[N] time;       // predictor

}

parameters {
  //Consider death. This is the proportion of mosquitoes dying (d_t) in treated huts (n_t)
  real alpha1;
  real alpha2;
}

model {
  real sp[N];
  
  alpha1 ~ normal(0,100);
  alpha2 ~ normal(0,100);
   
  for (n in 1:N) {
    sp[n] = alpha1  + alpha2 * time[n];
   }
  
  d_t ~ binomial_logit(n_t, sp);
 }

generated quantities{
  real sp_ppc[365];

    for(t in 1:365){
      sp_ppc[t] = binomial_rng(365, inv_logit(alpha1 + alpha2 * t)) / 365.0;
  }
}