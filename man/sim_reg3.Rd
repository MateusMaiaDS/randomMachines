\name{sim_reg3}
\alias{sim_reg3}
\title{Simulation for a regression toy examples from Random Machines Regression 3}
\usage{
sim_reg3(n, sigma)
}
\arguments{
\item{n}{Sample size}

\item{sigma}{Standard deviation of residual noise}

}
\value{
A simulated data.frame with two predictors and the target variable.
}
\description{
Simulation toy example initially found in Scornet (2016), and used and escribed by Ara \emph{et. al} (2022).
Inputs are 4 independent variables uniformly distributed on the interval \eqn{[-1,1]}. Outputs are generated following the equation
\deqn{Y= -\sin({X}_{1})+{X}^{2}_{2}+{X}_{3}-e^{{-X^{2}_{4}}} + \varepsilon, \varepsilon \sim \mathcal{N}(0,0.5)}
}
\references{

Ara, Anderson, et al. "Regression random machines: An ensemble support vector regression model with free kernel choice." Expert Systems with Applications 202 (2022): 117107.

Scornet, E. (2016). Random forests and kernel methods. IEEE Transactions on Information Theory, 62(3), 1485-1500.

}

\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}

\examples{
library(randomMachines)
sim_data <- sim_reg3(n=100)
}

