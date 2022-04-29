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
recipe <- nocovariates.from_conditional()
data <- nocovariates.treat(recipe)
```

``` r
data <- nocovariates.treat(recipe)
head(data)
```

    ##   i j A1 R A2           Y
    ## 1 1 1  1 1  1   5.0282097
    ## 2 1 2  1 1  1   0.2076947
    ## 3 1 3  1 1  1   4.0099962
    ## 4 1 4  1 1  1 -12.9397841
    ## 5 1 5  1 1  1  13.3996666
    ## 6 2 1  1 0  1  20.3180844

### Get the treatment summary of the generated data

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
estimates <- clustered.estimate(Y~a1+a2+I(a1*a2), data)
```

Get the estimated betas

``` r
estimates$beta_hat
```

    ##           [,1]
    ## [1,] 4.2998086
    ## [2,] 1.7644759
    ## [3,] 0.2578912
    ## [4,] 2.4001336

Get variance-covariance matrix of betas

``` r
estimates$cov_hat_beta_hat
```

    ##             [,1]        [,2]        [,3]        [,4]
    ## [1,]  0.16103459  0.02497196  0.02029938 -0.01140415
    ## [2,]  0.02497196  0.16103459 -0.01140415  0.02029938
    ## [3,]  0.02029938 -0.01140415  0.16103459  0.02497196
    ## [4,] -0.01140415  0.02029938  0.02497196  0.16103459

The formula that was called

``` r
estimates$formula
```

    ## Y ~ a1 + a2 + I(a1 * a2)
