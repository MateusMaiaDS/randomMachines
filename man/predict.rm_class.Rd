\name{predict.rm_class}
\alias{predict.rm_class}
\alias{predict,rm_class-method}
\title{Prediction function for the rm_class_model}
\usage{
\S4method{predict}{rm_class}(object,newdata)
}
\arguments{
\item{object}{A fitted RM model object of class \code{rm_class}.}

\item{newdata}{A data frame or matrix containing the new data to be predicted.}
}
\value{
A vector of predicted outcomes: probabilities in case of `prob_model = TRUE` and classes in case of `prob_model = FALSE`.
}
\description{
This function predicts the outcome for a RM object model using new data
}
\examples{
# Generating a sample for the simulation
library(randomMachines)
sim_data <- sim_class(n = 75)
sim_new <- sim_class(n = 25)
rm_mod <- randomMachines(y~., train = sim_data)
y_hat <- predict(rm_mod, newdata = sim_new)
}
