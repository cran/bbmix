// fit mixture of beta binomials to allelic counts based on known genotypes

// fit mixture of beta binomials to allelic counts based on known genotypes
data {
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of beta binomials
  int n[N]; // number of reads overlapping SNP alt allele
  int m[N]; // number of total reads overlapping SNP
  vector[K] alpha_p; //prior for alpha component
  vector[K] beta_p; // prior for beta
}

parameters {
  simplex[K] theta; 
  ordered[K] mu; // to avoid labelling switching impose order this way
  vector<lower=0.1>[K] lambda; // from Stan manual
  /* vector<lower=0>[K] alpha; //shape 1 */
  /* vector<lower=0>[K] beta; // shape 2 */
  
}

model {
  vector[K] log_theta = log(theta);
  /* vector[K] alpha_p = [1,10,19]'; */
  /* vector[K] beta_p = [19,10,1]'; */
  vector[K] lps; 
  // priors
  for (k in 1:K){
    mu[k] ~ beta(alpha_p[k],beta_p[k]);
    lambda[k] ~ gamma(alpha_p[k]+beta_p[k],1);  // mean a+b
  }
  
  theta ~ dirichlet(rep_vector(1.0,K));
  // likelihood
  for (i in 1:N) {
    lps = log_theta; 
    for (k in 1:K) {
      lps[k] += beta_binomial_lpmf(n[i] | m[i], lambda[k]*mu[k], lambda[k]*(1-mu[k]));
    }
    target += log_sum_exp(lps);
  }
}
