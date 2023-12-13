\name{sim_reg5}
\alias{sim_reg5}
\title{Simulation for a regression toy examples from Random Machines Regression 3}
\usage{
sim_reg5(n, sigma)
}
\arguments{
\item{n}{Sample size}

\item{sigma}{Standard deviation of residual noise}

}
\value{
A simulated data.frame with two predictors and the target variable.
}
\description{
Simulation toy example initially found in Van der Laan, \emph{et.al} (2016), and used and escribed by Ara \emph{et. al} (2022).
Inputs are 6 independent variables sampled from \eqn{N(0,1)}. Outputs are generated following the equation
\deqn{Y=X_{1}+0.707 X^{2}_{2} + 2\mathcal{1}_{(X_{3}>0)}+0.873 \log (X_{1})|X_{3}|+0.894 X_{2} X_{4}+2\mathcal{1}_{(X_{5}>0)}+0.464e^{X_{6}}+ \varepsilon, \varepsilon \sim \mathcal{N}(0,\sigma^{2})}
}
\references{

Ara, Anderson, et al. "Regression random machines: An ensemble support vector regression model with free kernel choice." Expert Systems with Applications 202 (2022): 117107.

Roy, M. H., & Larocque, D. (2012). Robustness of random forests for regression. Journal of Nonparametric Statistics, 24(4), 993-1006.
}

\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}

\examples{
library(randomMachines)
sim_data <- sim_reg5(n=100)
}

