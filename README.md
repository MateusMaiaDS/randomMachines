# randomMachines
## Installation

You can install the development version of randomMachines from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MateusMaiaDS/randomMachines")
```

## Example

This is a basic example which shows you how to solve a common binary classification problem:

``` r
library(randomMachines)
## Simple classification example
sim_train <- randomMachines::sim_class(n=100)
sim_test <- randomMachines::sim_class(n=100)
rm_mod <- randomMachines::randomMachines(y~.,train = sim_trainboots_size = 25,prob_model = F)
rm_mod_pred <- predict(rm_mod,sim_test)
```


For a regression task we would have similarly

``` r
library(randomMachines)
## Simple regression example
sim_train <- randomMachines::sim_reg(n=100)
sim_test <- randomMachines::sim_reg(n=100)
rm_mod <- randomMachines::randomMachines(y~.,train = sim_train,boots_size = 25)
rm_mod_pred <- predict(rm_mod,sim_test)
```
