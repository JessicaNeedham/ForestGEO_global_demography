data {
  int<lower=0> N;       // no. of trees 
  real incr[N];         // diameter increment in mm
  real groups[N];   // which growth groups
}

parameters {
  real<lower=0> alpha1;
  real<lower=0> alpha2;
  real<lower=0> beta1;
  real<lower=0> beta2;
}

model {
  for(i in 1:N)
  {
    if(groups[i] == 1)
      target += (gamma_lpdf(incr[i] | alpha1, beta1));
    else
      target += (gamma_lpdf(incr[i] | alpha2, beta2));
  }
}

