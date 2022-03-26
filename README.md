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

    ##   i j A1 R A2         Y
    ## 1 1 1  1 1  1 24.662398
    ## 2 1 2  1 1  1 22.581827
    ## 3 1 3  1 1  1  3.480175
    ## 4 1 4  1 1  1 18.550622
    ## 5 1 5  1 1  1  6.832073
    ## 6 2 1  1 0  1 -2.041407

### Get treatment summary of the generated data

``` r
treatment_summary(recipe)
```

    ##         treat_mu treat_var  treat_cor betas
    ## (1, 1)       8.6    100.84 0.17691392  3.85
    ## (1, -1)      4.4    113.44 0.23748237  2.65
    ## (-1, 1)     -0.8    105.76 0.16036309  0.05
    ## (-1,-1)      3.2    100.16 0.08146965  2.05

### Estimate parameters from data

``` r
estimates <- nocovariates.estimate(data)
```

``` r
estimates
```

    ## $beta_hat
    ##           [,1]
    ## [1,] 4.5457326
    ## [2,] 1.9717549
    ## [3,] 0.2690568
    ## [4,] 2.2854796
    ## 
    ## $var_hat
    ## $var_hat$`1,1`
    ## [1] 91.34229
    ## 
    ## $var_hat$`1,-1`
    ## [1] 115.5893
    ## 
    ## $var_hat$`-1,1`
    ## [1] 109.4682
    ## 
    ## $var_hat$`-1,-1`
    ## [1] 105.2085
    ## 
    ## 
    ## $rho_hat
    ## $rho_hat$`1,1`
    ## [1] 0.08410177
    ## 
    ## $rho_hat$`1,-1`
    ## [1] 0.15955
    ## 
    ## $rho_hat$`-1,1`
    ## [1] 0.13929
    ## 
    ## $rho_hat$`-1,-1`
    ## [1] 0.06764574
    ## 
    ## 
    ## $cov_hat_beta_hat
    ##              [,1]         [,2]         [,3]         [,4]
    ## [1,]  0.125390212 -0.001382212  0.001485363 -0.003719345
    ## [2,] -0.001382212  0.125390212 -0.003719345  0.001485363
    ## [3,]  0.001485363 -0.003719345  0.125390212 -0.001382212
    ## [4,] -0.003719345  0.001485363 -0.001382212  0.125390212
