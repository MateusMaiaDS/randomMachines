# rmachines
## Installation

You can install the development version of rmachines from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MateusMaiaDS/rmachinesproject")
```

## Example

This is a basic example which shows you how to solve a common binary classification problem:

``` r
library(rmachines)
## Simple classification example
sim_train <- rmachines::sim_class(n=100)
sim_val <- rmachines::sim_class(n=25)
sim_test <- rmachines::sim_class(n=100)
rm_mod <- rmachines::random_machines(y~.,train = sim_train,validation = sim_val,boots_size = 25,prob_model = F)
rm_mod_pred <- predict(rm_mod,sim_test)
```


For a regression task we would have similarly

``` r
library(rmachines)
## Simple regression example
sim_train <- rmachines::sim_reg(n=100)
sim_val <- rmachines::sim_reg(n=25)
sim_test <- rmachines::sim_reg(n=100)
rm_mod <- rmachines::random_machines(y~.,train = sim_train,validation = sim_val,boots_size = 25)
rm_mod_pred <- predict.rm_model_reg(rm_mod,sim_test)
```
