# multi-site-BIRD

## Requirments:
1. Please download and extract STANINPUTS.tar.gz and STANOUTPUTS.tar.gz from https://prodduke-my.sharepoint.com/:u:/g/personal/sl548_duke_edu/ET0c_cDRkzVDnT7vFbTgwQ0BC4fq9cNGR4WSI-fKhbaE4w?e=P03Tv9 and https://prodduke-my.sharepoint.com/:u:/g/personal/sl548_duke_edu/EXhZBSp7HwBHiofY9Xd3KlEBUWcZ9f_7G95aRlZKMBIQTQ?e=wHByU3. 
2. Please download cmdstan-2.27.0 from https://github.com/stan-dev/cmdstan/releases 

## cmdStan: Multi-site BIRD model and its other derived models.
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
    model{
    ...
    theta[j] ~ normal(1,1);
    ...
    }
    ```
    
7. `multiple2_merge_function_normal3.stan`
    - Modification: Turn `qi`,`p` and `theta` into normal distribution through logit and log transformation
    ```
    parameters {
    ...
        real p_logit[N_VARIANTS];
        real qi_logit[N_VARIANTS,N_RNA];
        real theta_log[N_VARIANTS];
    ...
    }
    transformed parameters { // ORDER MATTERS!
        real<lower=0,upper=1> p[N_VARIANTS]; // untransformed p
        for(j in 1:N_VARIANTS)
            p[j]=inv_logit(p_logit[j]);
        real<lower=0,upper=1> qi[N_VARIANTS,N_RNA]; //untransformed qi
        for(j in 1:N_VARIANTS){
            for(i in 1:N_RNA){
                qi[j,i] = inv_logit(qi_logit[j,i]);
            }
        }
        real<lower=0> theta[N_VARIANTS];
        for(j in 1:N_VARIANTS){
            theta[j]=exp(theta_log[j]);
        ...
   }
   ```
8. `qi_q.stan`
    - BIRD model breakdown, to simply infer `q` from `qi`

## cmdStan: Input and Output Naming (__IMPORTANT!!__)
1. `STANINPUTS/0.5_data_0.01_1000_read100.txt`
    - This a human readable simulated input where `theta` is simulated as 0.5 for regulatory variants, `0.01` means there are 1 percent regulatory variants (when this number = 1.0, this means all of variants are regulatory variants), `1000` means there are 1000 variants in total, `read100` means 100 reads for DNA and 100/10=10 reads (10 is the number of RNA sites) for each RNA pair. If there is something like `read5050`, this means there are 50 reads for DNA and 50 reads for each of RNA pair. 
2. `STANINPUTS/0.5_data_0.01_1000_read100.txt.staninputs`
    - Stan readable inputs.
3. `STANOUTPUTS/multiple2_merge_function_new_VI_0.5_data_0.01_1000_read100.txt.staninputs.stanoutputs`
    - Stan outputs for this input. `multiple2_merge_function_new` is the version name of the model and `VI` is the algorithm (variational inference) or it could be `MCMC`(NUTS). 

## cmdStan: Other accessory files
1. `data_simulator_together.ipynb` 
    - Used to simulate data (__`Required_Python_Packages` needed__)
    - __Please see in the notebook for details of model version control and input control.__
        
        a. The data of all variants to have same simulated `theta`(usually `theta=2`). To investigate the bias of different algorithm. 
      
        b. The data of small percentage of regulatory variants(usually `theta=2`) and rest of null variants(`theta=1`). To investigate the accuracy of the model and algorithm. 
2. `sim-equal.R`
    - Used in `data_simulator_together.ipynb` to simulate data. 
3. `cmdStan_results_analysis.ipynb`
    - Used to run different models in different algorithms(MCMC or VI in Stan)
    - Used to plot distribution of all variants `theta` posterior medians. (Use data from __1a__)
    - Used to plot auc_roc (Need to use data from __1b__)
    - __Please see the details of model version control and inputs control in the notebook__
4. `compile.py`
    - Used to compile Stan model in current directory. 
    - `./compile.py xxx` (for xxx.stan)
    - __IMPORTANT! make sure to edit the directories in `compile.py` to be the directories where your models are__.

## cmdStan: Other accessory directories 
1. `STANINPUTS`
    - To store all simulated inputs including the visible inputs and inputs for Stan. 
2. `STANOUTPUTS`
    - To store all Stan's outputs 
3. `Required_Python_Packages`
    - All __required__ customized packages for this project. 
4. `cmdstan-2.27.0`
    - cmdStan for this project.
5. `P_simulation`
    - simulation to investigate `p`'s effect on logit transformation. 

## RStan (MCMC)
1. `rstan.Rmd`
    - Use RStan to run 4 chains independently in order to compare with RJAGS.
    - 6 different datasets for `multiple2_merge_function.stan` model. 
        a. `stan_model_30`:1000 variants with `theta=2`, DNA depth each pair = RNA depth each pair = 30 (To test bias)
        
        b. `stan_model_50`:1000 variants with `theta=2`, DNA depth each pair = RNA depth each pair = 50 (To test bias)
        
        c. `stan_model_001_30`: 1000 variants, 1 percent are regulatory with `theta=2`, 99 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 30 (To test accuracy)
        
        d. `stan_model_001_50`: 1000 variants, 1 percent are regulatory with `theta=2`, 99 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 50 (To test accuracy)
        
        e. `stan_model_005_30`: 1000 variants, 5 percent are regulatory with `theta=2`, 95 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 30 (To test accuracy)
        
        f. `stan_model_005_50`: 1000 variants, 5 percent are regulatory with `theta=2`, 95 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 50 (To test accuracy)

## RJAGS (MCMC)
1. `JAGS/rjags.Rmd`
    - Use RJAGS to run 4 chains independently in order to compare with RStan.
    - 6 different datasets for `multiple2_merge_function.stan` model. 
        a. `model_30`:1000 variants with `theta=2`, DNA depth each pair = RNA depth each pair = 30 (To test bias)
        
        b. `model_50`:1000 variants with `theta=2`, DNA depth each pair = RNA depth each pair = 50 (To test bias)
        
        c. `model_001_30`: 1000 variants, 1 percent are regulatory with `theta=2`, 99 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 30 (To test accuracy)
        
        d. `model_001_50`: 1000 variants, 1 percent are regulatory with `theta=2`, 99 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 50 (To test accuracy)
        
        e. `model_005_30`: 1000 variants, 5 percent are regulatory with `theta=2`, 95 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 30 (To test accuracy)
        
        f. `model_005_50`: 1000 variants, 5 percent are regulatory with `theta=2`, 95 percent are null variants with `theta=1`, DNA depth each pair = RNA depth each pair = 50 (To test accuracy)


## Tensorflow-probability (Variational Inference)
1. `simple_VI_tfp.ipynb`
2. The variational inference feature is only in `tfp-nightly` or the newest tensorflow-probability version (0.13.0+).
3. check https://www.tensorflow.org/probability/examples/TFP_Release_Notebook_0_13_0 or https://www.tensorflow.org/probability/examples/Variational_Inference_and_Joint_Distributions or http://hyperion.usc.edu/UQ-SummerSchool/pres/Dillon.pdf for further variational inference tutorial in TFP. 

## Pyro (Variational Inference)
1. `simple_VI_pyro.ipynb`
2. check http://pyro.ai/examples/svi_part_i.html for futher stochastic variational inference tutorial in Pyro.

