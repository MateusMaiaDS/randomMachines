\name{brier_score}
\alias{brier_score}
\title{Brier Score function}
\usage{
brier_score(prob, observed, levels)
}
\arguments{
\item{prob}{predicted probabilities}

\item{observed}{\eqn{y} observed values (it assumed that the positive class is coded is equal to one and the negative 0)}

\item{levels}{A string vector with the original levels from the target variable}
}
\value{
     Returns the Brier Score, a numeric value indicating the accuracy of the predictions.
}
\description{
 Calculate the Brier Score for a set of predicted probabilities and observed outcomes.
 The Brier Score is a measure of the accuracy of probabilistic predictions. It is commonly used in the evaluation of predictive models.
 }
