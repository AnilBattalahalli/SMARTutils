### Install the library

``` r
library(devtools)
install_github("AnilBattalahalli/SMARTutils")
```


### Import the library

``` r
library(SMARTutils)
```

### Generate data from conditional parameters

``` r
recipe <- nocovairiates.from_conditional(cell_mu = c(10,8,2,4,-2,3), cell_var = c(100,100,100,100,100,100),
                                           cell_cor = c(0.1, 0.2, 0.15, 0.12, 0.11, 0.07), N=200, m=5, p_1=0.3, p_2=0.2)
data <- nocovariates.treat(recipe)
```

``` r
head(data)
```

    ##   i j A1 R A2           Y
    ## 1 1 1  1 0  1   6.6070590
    ## 2 1 2  1 0  1   7.3485029
    ## 3 1 3  1 0  1 -10.8454299
    ## 4 1 4  1 0  1   0.9431890
    ## 5 1 5  1 0  1  -1.7065543
    ## 6 2 1 -1 1  1  -0.9599292

### Get the treatment summary of the generated data

``` r
treatment_summary(recipe)
```

    ##         treat_mu treat_var  treat_cor betas
    ## (1, 1)       8.6    100.84 0.17691392  3.85
    ## (1, -1)      4.4    113.44 0.23748237  2.65
    ## (-1, 1)     -0.8    105.76 0.16036309  0.05
    ## (-1,-1)      3.2    100.16 0.08146965  2.05

### Get the true effect size of two DTRs

``` r
nocovariates.true_effect_size(recipe, '1,1', '1,-1')
```

    ## [1] 0.405764

### Estimate parameters from data

``` r
estimates <- nocovariates.estimate(data)
```

``` r
estimates
```

    ## $beta_hat
    ##            [,1]
    ## [1,]  3.9494134
    ## [2,]  3.2022105
    ## [3,] -0.0290067
    ## [4,]  1.7596609
    ## 
    ## $var_hat
    ## $var_hat$`1,1`
    ## [1] 102.8521
    ## 
    ## $var_hat$`1,-1`
    ## [1] 123.0341
    ## 
    ## $var_hat$`-1,1`
    ## [1] 109.0926
    ## 
    ## $var_hat$`-1,-1`
    ## [1] 108.178
    ## 
    ## 
    ## $rho_hat
    ## $rho_hat$`1,1`
    ## [1] 0.2394931
    ## 
    ## $rho_hat$`1,-1`
    ## [1] 0.0876658
    ## 
    ## $rho_hat$`-1,1`
    ## [1] 0.1158456
    ## 
    ## $rho_hat$`-1,-1`
    ## [1] 0.1155997
    ## 
    ## 
    ## $cov_hat_beta_hat
    ##             [,1]        [,2]        [,3]        [,4]
    ## [1,] 0.145137408 0.009633173 0.002219290 0.015796011
    ## [2,] 0.009633173 0.145137408 0.015796011 0.002219290
    ## [3,] 0.002219290 0.015796011 0.145137408 0.009633173
    ## [4,] 0.015796011 0.002219290 0.009633173 0.145137408
    ## 
    ## $mu_hat_dtr
    ##           [,1]
    ## [1,]  8.882278
    ## [2,]  5.420970
    ## [3,] -1.041465
    ## [4,]  2.535870

### Estimate the effect size of two DTRs

``` r
nocovariates.estimate_effectsize(estimates, '1,1', '1,-1')
```

    ## [1] 0.3256945

### Estimate the effect size of all DTR pairs

``` r
nocovariates.estimate_effectsizes(estimates)
```

    ## $`(1,1), (1,-1)`
    ## [1] 0.3256945
    ## 
    ## $`(1,1), (-1,1)`
    ## [1] 0.964005
    ## 
    ## $`(1,1), (-1,-1)`
    ## [1] 0.6178326
    ## 
    ## $`(1,-1), (-1,1)`
    ## [1] 0.5998583
    ## 
    ## $`(1,-1), (-1,-1)`
    ## [1] 0.2683308
    ## 
    ## $`(-1,1), (-1,-1)`
    ## [1] -0.3432213

### Estimate the variance of the effect size of two DTRs

``` r
nocovariates.estimate_var_effectsize(estimates, '1,1', '-1,1')
```

    ## [1] 0.01112416

### Draw bootstrap samples from the clusters

``` r
head(clusterBoot(data))
```

    ##   i j A1 R A2         Y
    ## 1 1 1 -1 0  1  9.665735
    ## 2 1 2 -1 0  1  5.380086
    ## 3 1 3 -1 0  1 -7.883186
    ## 4 1 4 -1 0  1 -1.432520
    ## 5 1 5 -1 0  1  8.417216
    ## 6 2 1  1 0  1  8.492830
