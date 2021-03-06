#include /Users/sean00002/Dropbox/PhD/Andrew_Rotation/functions.stan

data {
   int<lower=1> N_VARIANTS;             // number of variants in the model
   real<lower=0,upper=1> v[N_VARIANTS]; // alt allele freq in VCF file
   int<lower=1> N_DNA;                  // number of DNA replicates
   int<lower=0> a[N_VARIANTS,N_DNA];   // DNA alt read counts
   int<lower=0> b[N_VARIANTS,N_DNA];   // DNA ref read counts
   int<lower=1> N_RNA;                  // number of RNA replicates
   int<lower=0> k[N_VARIANTS,N_RNA];   // RNA alt read counts
   int<lower=0> m[N_VARIANTS,N_RNA];   // RNA ref read counts
}

parameters {
   real<lower=0,upper=1> p[N_VARIANTS]; // alt allele freq in DNA
   real<lower=0,upper=1> qi[N_VARIANTS,N_RNA]; // alt freqs in RNA
   real<lower=0> theta[N_VARIANTS];
   real<lower=2> c1; // concentration parameter of beta prior for qi
   real<lower=2> c2; // concentration parameter of beta prior for p
   real<lower=0> s; // std.dev. parameter of lognormal prior for theta
}

transformed parameters { // ORDER MATTERS!
   real<lower=0,upper=1> q[N_VARIANTS]; // alt allele freq in RNA
   for(j in 1:N_VARIANTS)
      q[j]=theta[j]*p[j]/(1.0-p[j]+theta[j]*p[j]);
}

model {
   // Parameters
   c1 ~ gamma(1.1, 0.005);
   c2 ~ gamma(1.1, 0.005);
   s ~ gammaModeSD(1,1);
   for(j in 1:N_VARIANTS) {
      p[j] ~ betaModeConc(v[j],c2);
      log(theta[j])/s ~ normal(0,1); // theta~lognormal(0,s)
      target+=-log(theta[j])-log(s); // Jacobian
      for(i in 1:N_RNA)
         qi[j,i] ~ betaModeConc(q[j],c1);
      for(i in 1:N_DNA)
         a[j,i] ~ binomial(a[j,i]+b[j,i],p[j]);
      for(i in 1:N_RNA)
         k[j,i] ~ binomial(k[j,i]+m[j,i],qi[j,i]);
   }
}
