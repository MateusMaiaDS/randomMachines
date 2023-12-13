\name{sim_reg4}
\alias{sim_reg4}
\title{Simulation for a regression toy examples from Random Machines Regression 3}
\usage{
sim_reg4(n, sigma)
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
Inputs are 6 independent variables uniformly distributed on the interval \eqn{[-1,1]}. Outputs are generated following the equation
\deqn{Y={X^{2}_{1}}+{X}^{2}_{2}{X_{3}}e^{-|{X_{4}}|}+{X_{6}}-{X_{5}}+ \varepsilon, \varepsilon \sim \mathcal{N}(0,\sigma^{2})}
}
\references{

Ara, Anderson, et al. "Regression random machines: An ensemble support vector regression model with free kernel choice." Expert Systems with Applications 202 (2022): 117107.

Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. Statistical applications in genetics and molecular biology, 6(1).
}

\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}

\examples{
library(randomMachines)
sim_data <- sim_reg4(n=100)
}

