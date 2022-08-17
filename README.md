### Import the library

``` r
install.packages("devtools")
library(devtools)
```

``` r
install_github("AnilBattalahalli/SMARTutils")
library(SMARTutils)
```

### Generate data from conditional parameters

``` r
data <- cSMART.dgen()
```

    ##  DTR      Treatment Mean      Treatment Variance      Treatment Correlation 
    ## ______    ______________      __________________      _____________________ 
    ##  1, 1         8.6           100.84          0.1769139 
    ##  1,-1         4.4           113.44          0.2374824 
    ## -1, 1         -0.8          105.76          0.1603631 
    ## -1,-1         3.2           100.16          0.08146965 
    ## ______    ______________      __________________      _____________________ 
    ## 
    ## True Beta:
    ## 
    ## Beta[0] - 3.85 
    ## Beta[1] - 2.65 
    ## Beta[2] - 0.05 
    ## Beta[3] - 2.05

``` r
head(data)
```

    ##   i j A1 R A2         Y
    ## 1 1 1 -1 0  1 -3.452152
    ## 2 1 2 -1 0  1 -2.735390
    ## 3 1 3 -1 0  1 -1.758195
    ## 4 1 4 -1 0  1  8.243450
    ## 5 1 5 -1 0  1  3.232144
    ## 6 2 1 -1 0  1 13.987650

### Get the estimates for the regression parameters

``` r
report <- cSMART.mm(Y~a1+a2+I(a1*a2), data, verbose=T)
```

    ## Parameter     Estimate    Std.Err     Z Score     Pr(>|z|) 
    ## ___________   --------    -------     -------     ------- 
    ## (Intercept)   4.47717     0.36671     12.20888    2.786906e-34 
    ## a1            2.74937     0.36671     7.497308    6.514162e-14 
    ## a2            -0.22003    0.36671     -0.5999946      0.5485099 
    ## I(a1 * a2)    1.10251     0.36671     3.006454    0.00264314 
    ## ___________   --------    -------     -------     ------- 
    ## 
    ## Marginal Mean Model:  Y ~ a1 + a2 + I(a1 * a2) 
    ## 
    ## Working covariance structure: 'EXCH'  (Homogeneous-Exchangeable covariance structure)
    ## Variance      102.2465
    ## Correlation   0.1591008
    ## 
    ## Variance-Covariance matrix of the estimates
    ##              (Intercept)           a1           a2   I(a1 * a2)
    ## (Intercept)  0.134479319  0.003630063 -0.013093342 -0.004863585
    ## a1           0.003630063  0.134479319 -0.004863585 -0.013093342
    ## a2          -0.013093342 -0.004863585  0.134479319  0.003630063
    ## I(a1 * a2)  -0.004863585 -0.013093342  0.003630063  0.134479319
