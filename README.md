## SMARTutils
### Install the package using the 'install_github' method
```
library(devtools)
install_github("AnilBattalahalli/SMARTutils")
```
### Example:

```
recipe <- nocovairiates.from_conditional()
generated_data <- nocovariates.treat(recipe)
```
Generate new data from the same recipe

```
data_1 <- nocovariates.treat(recipe)
data_2 <- nocovariates.treat(recipe)
```
Get treatment summary

```
summary <- treatment_summary(recipe)
```

```
var <- list("1,1" = 100, "1,-1"=100, "-1,1"=100, "-1,-1"=100)
rho <- list("1,1" = 0.2, "1,-1"=0.1, "-1,1"=0.1, "-1,-1"=0.3)
nocovariates.estimate_betas(var, rho, generated_data)
```
