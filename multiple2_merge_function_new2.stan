functions {
   // This function can be used via: g~gammaModeSD(mode,sd);
   real gammaModeSD_lpdf(real parm,real m,real sd) {
      real r=(m+sqrt(m^2+4*sd^2))/(2*sd^2);
      real s=1+m*r;
      return gamma_lpdf(parm|r,s);
   }

   // This function can be used via: p~betaModeConc(mode,concentration);
   real betaModeConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2)+1, (1-m)*(c-2)+1);
   }

   real betaMeanConc_lpdf(real parm,real m,real c) {
      return beta_lpdf(parm|m*(c-2), (1-m)*(c-2));
   }

   real dirichMultinom_lpmf(int[] x,int K,vector alpha) {
      real alpha0=sum(alpha);
      int n=sum(x);
      real logP=log(n)+lbeta(alpha0,n);
      for(k in 1:K)
         if(x[k]>0) logP-=(x[k]+lbeta(alpha[k],x[k]));
      return logP; }
}

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
   //real<lower=0,upper=1> p[N_VARIANTS]; // alt allele freq in DNA
   real<lower=0,upper=1> qi[N_VARIANTS,N_RNA]; // alt freqs in RNA
   real<lower=0> theta[N_VARIANTS];
   real<lower=0,upper=1> q[N_VARIANTS];
   real<lower=2> c1; // concentration parameter of beta prior for qi
   real<lower=2> c2; // concentration parameter of beta prior for p
   real<lower=0> s; // std.dev. parameter of lognormal prior for theta
}

transformed parameters { // ORDER MATTERS!
   //real<lower=0,upper=1> q[N_VARIANTS]; // alt allele freq in RNA
   //real<lower=0> theta[N_VARIANTS];
   real<lower=0,upper=1> p[N_VARIANTS];
   for(j in 1:N_VARIANTS)
      //theta[j] = (q[j]/p[j])/((1-q[j])/(1-p[j]));
      //q[j]=theta[j]*p[j]/(1.0-p[j]+theta[j]*p[j]);
      p[j]= q[j]/(theta[j]-theta[j]*q[j]+q[j]);
}

model {
   // Parameters
   c1 ~ gamma(1.1, 0.0005);
   c2 ~ gamma(1.1, 0.0005);
   s ~ gamma(1.1,3);
   for(j in 1:N_VARIANTS) {
      //p[j] ~ betaModeConc(v[j],c2);
      log(theta[j])/s ~ normal(0,1); // theta~lognormal(0,s)
      target+=-log(theta[j])-log(s); // Jacobian
      q[j] ~ beta(2,2);
      for(i in 1:N_RNA)
         qi[j,i] ~ betaModeConc(q[j],c1);
      for(i in 1:N_DNA)
         a[j,i] ~ binomial(a[j,i]+b[j,i],p[j]);
      for(i in 1:N_RNA)
         k[j,i] ~ binomial(k[j,i]+m[j,i],qi[j,i]);
   }
}
