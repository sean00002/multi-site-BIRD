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


