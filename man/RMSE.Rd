% Please edit documentation in R/main_functions.R
\name{RMSE}
\alias{RMSE}
\title{Root Mean Squared Error (RMSE) Function}
\usage{
RMSE(predicted, observed)
}
\arguments{
\item{predicted}{A vector of predicted values \eqn{\hat{\mathbf{y}}}.}

\item{observed}{A vector of observed values \eqn{\mathbf{y}}.}
}
\value{
a the Root Mean Squared error calculated by the formula in the description.
}
\description{
Computes the Root Mean Squared Error (RMSE), a widely used metric for evaluating the accuracy of predictions in regression tasks. The formula is given by
\deqn{RMSE = \sqrt{\frac{1}{n}\sum_{i=1}^{n}\left(y_{i}-\hat{y}_{i}\right)^{2}}}
}
