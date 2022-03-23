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
