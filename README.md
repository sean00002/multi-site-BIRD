# multi-site-BIRD
## Introduction: Multi-site BIRD model and its other derived models.
1. `multiple2.stan` 
    - The originla multi-site BIRD model 
2. `functions.stan` 
    - The functions for `multiple2.stan`
3. `multiple2_merge_function.stan`
    - the model merged both `multiple2.stan` and `functions.stan`. All derived models are based on this one in order to prevent bugs of cmdStan reading function block file. 
4. `multiple2_merge_function_new.stan`
    - Modified model: Use `theta` as deterministic parameter instead of `q`
