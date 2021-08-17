
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
}


data {
   int<lower=1> N_RNA;                  // number of RNA replicates
   int<lower=1> N_VARIANTS;             // number of variants in the model
   real<lower=0,upper=1> qi[N_VARIANTS,N_RNA]; // alt allele freq in VCF file
}

parameters {
   real<lower=0,upper=1> q[N_VARIANTS];
   real<lower=2> c2; // concentration parameter of beta prior for p
}


model {
   // Parameters
   c2 ~ gamma(1.1, 0.005);
   for(j in 1:N_VARIANTS) {
      q[j] ~ uniform(0,1);
      for(i in 1:N_RNA)
         qi[j,i] ~ betaModeConc(q[j],c2);
   }
}
