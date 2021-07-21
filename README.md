# multi-site-BIRD
## Introduction: Multi-site BIRD model and its other derived models.
1. `multiple2.stan` 
    - The originla multi-site BIRD model 
2. `functions.stan` 
    - The functions for `multiple2.stan`
3. `multiple2_merge_function.stan`
    - the model merged both `multiple2.stan` and `functions.stan`. All derived models are based on this one in order to prevent bugs of cmdStan reading function block file. 
4. `multiple2_merge_function_new.stan`
    - Modification: Use `theta` as deterministic parameter instead of `q`
    ```
    transformed parameters {
        real<lower=0> theta[N_VARIANTS];
        for(j in 1:N_VARIANTS)
            theta[j] = (q[j]/p[j])/((1-q[j])/(1-p[j]));
    }
    ```
5. `multiple2_merge_function_new2.stan`
    - Modification: Use `p` as deterministic parameter instead of `q`
    ```
    transformed parameters { 
        real<lower=0,upper=1> p[N_VARIANTS];
        for(j in 1:N_VARIANTS)
            p[j]= q[j]/(theta[j]-theta[j]*q[j]+q[j]);
    }
    ```
6. `multiple2_merge_function_new3.stan`
    - Modification: Give `theta` a fixed prior. 
    ```
     theta[j] ~ normal(1,1);
    ```
7. `multiple2_merge_function_new3.stan`
    - Modification: Give `theta` a fixed prior. 
    ```
     theta[j] ~ normal(1,1);
    ```
8. `multiple2_merge_function_normal.stan`
    - Modification: Turn q into normal distribution through logit transformation
    ```
    transformed parameters { // ORDER MATTERS!
        real q_logit[N_VARIANTS]; // alt allele freq in RNA
        for(j in 1:N_VARIANTS)
            q_logit[j]=logit((theta[j]*p[j])/(1.0-p[j]+theta[j]*p[j]));
    ...
    for(i in 1:N_RNA)
         qi[j,i] ~ betaModeConc(inv_logit(q_logit[j]),c1);
    ...
}
