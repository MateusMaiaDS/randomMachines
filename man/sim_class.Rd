% Please edit documentation in R/simulation_examples.R
\name{sim_class}
\alias{sim_class}
\title{Generate a binary classification data set from normal distribution}
\usage{
sim_class(
  n,
  p = 2,
  ratio = 0.5,
  mu_a = 0,
  sigma_a = 1,
  mu_b = 1,
  sigma_b = 1
)
}
\arguments{
\item{n}{Sample size}

\item{p}{Number of predictors}

\item{ratio}{Ratio between class A and class B}

\item{mu_a}{Mean of \eqn{X_{1}}.}

\item{sigma_a}{Standard deviation of \eqn{X_{1}}.}

\item{mu_b}{Mean of \eqn{X_{2}}}

\item{sigma_b}{Standard devation of \eqn{X_{2}}}

}
\value{
A simulated data.frame with two predictors for a binary classification problem
}
\description{
Simulation used as example of a classification task based on a separation of two
normal multivariate distributions with different vector of means and differerent covariate matrices.
For the label \eqn{A} the \eqn{\mathbf{X}_{A}} are sampled from a normal distribution \eqn{{MVN}\left(\mu_{A}\mathbf{1}_{p},\sigma_{A}^{2}\mathbf{I}_{p}\right)}
while for label \eqn{B} the samples \eqn{\mathbf{X}_{B}} are from a normal distribution \eqn{{MVN} \left(\mu_{B}\mathbf{1}_{p},\sigma_{B}^{2}\mathbf{I}_{p}\right)}. For more details see Ara \emph{et. al} (2021), and Breiman L (1998).
}

\references{

Ara, Anderson, et al. "Random machines: A bagged-weighted support vector model with free kernel choice." Journal of Data Science 19.3 (2021): 409-428.

Breiman, L. (1998). Arcing classifier (with discussion and a rejoinder by the author). The annals of statistics, 26(3), 801-849.
}

\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}
\examples{
library(randomMachines)
sim_data <- sim_class(n = 100)
}
