\name{predict.rm_reg}
\alias{predict.rm_reg}
\alias{predict,rm_reg-method}
\title{Prediction function for the rm_reg_model}
\usage{
\S4method{predict}{rm_reg}(object,newdata)
}
\arguments{
\item{object}{A fitted RM model object of class \code{rm_reg}.}

\item{newdata}{A data frame or matrix containing the new data to be predicted.}

}
\value{
Predicted values \code{newdata} object from the Random Machines model.
}
\description{
This function predicts the outcome for a RM object model using new data for continuous \eqn{y}
}
\examples{
# Generating a sample for the simulation
library(randomMachines)
sim_data <- sim_reg1(n = 75)
sim_new <- sim_reg1(n = 25)
rm_mod_reg <- randomMachines(y~., train = sim_data)
y_hat <- predict(rm_mod_reg, newdata = sim_new)
}
