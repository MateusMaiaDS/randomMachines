\name{random_machines}
\alias{random_machines}
\alias{predict.rm_class}
\alias{predict.rm_reg}
\title{Random Machines}
\description{Random Machines is an ensemble model which use the combination of different kernel functions to improve the diversity in the bagging approach improving the predictions.

\itemize{
     \item For a binary classification problem \eqn{\mathbin{{ G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right), }}{\mathbin{{ G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right). }}}}
     \item For a numeric response \eqn{G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right)}{G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right)}
}

}

\usage{
random_machines(
     y~.,
     train,validation,
     B = 25, cost = 1,
     seed.bootstrap = NULL,
     automatic_tuning = FALSE,
     gamma_rbf = 1,
     gamma_lap = 1,
     degree = 2,
     poly_scale = 1,
     offset = 0,
     gamma_cau = 1,
     d_t = 2,
     kernels = c("rbfdot", "polydot", "laplacedot", "vanilladot", "cauchydot"),
     prob_model = T,
     loss_function = RMSE,
     epsilon = 0.1,
     beta = 2
)
}

\arguments{
     \item{formula}{
     formula an object of class \link[=formula]{formula}: it should contain a symbolic description of the model to be fitted, indicating the dependent variable and all predictors that should be included.
     }
     \item{train}{
     the training data \eqn{\left\{\left( \mathbf{x}_{i},y_{i} \right)\right\}_{i=1}^{N}} used to train the model.
     }
     \item{validation}{
     the validation data \eqn{\left\{\left( \mathbf{x}_{i},y_{i}\right)  \right\}_{i=1}^{V}} used to calculate probabilities \eqn{\lambda_{r}}.
     }
     \item{B}{
     number of bootstrap samples.
     }
     \item{cost}{
     the \eqn{C}-constant term of the regularization on the Lagrange formulation.
     }
     \item{seed.boostrap}{
     setting a seed to replicate bootstrap sampling. The default value is \code{{NULL}}.
     }
     \item{automatic_tuning}{
     boolean to define if the kernel hyperparameters will be selected using the \code{sigest} from the \code{ksvm} function.
     }
     \item{gamma_rbf}{
     the hyperparameter \eqn{\gamma_{RBF}} used in the RBF kernel.
     }
     \item{gamma_lap}{
     the hyperparameter \eqn{\gamma_{LAP}} used in the Laplacian kernel.
     }
     \item{degree}{
     the degree used in the Polynomial kernel.
     }
     \item{poly_scale}{
     the scale parameter from Polynomial kernel.
     }
     \item{offset}{
     the offset parameter from the Polynomial kernel.
     }
     \item{gamma_cau}{
     the offset arameter from the Cauchy kernel.
     }
     \item{d_t}{
     ??
     }
     \item{kernels}{
     a vector with the name of kernel functions that will be used in the Random Machines model. The default includes all the functions: \code{c("rbfdot", "polydot", "laplacedot", "vanilladot", "cauchydot").}
     }
     \item{prob_model}{
     a boolean to define if the algorithm will be using a probabilistic approach to the define the predictions (default = \code{T}).
     }
     \item{loss_function}{
     Define which loss function is going to be used in the regression approach. The default is the \code{RMSE} function but other can be added.
     }
     \item{epsilon}{
     The epsilon in the loss function used from the SVR implementation.
     }
     \item{beta}{
     The correlation parameter \eqn{\beta} which calibrate the penalisation of each kernel performance.
     }
}
\details{
     bla bla bla bla
}
\references{

Ara, Anderson, et al. "Regression random machines: An ensemble support vector regression model with free kernel choice." Expert Systems with Applications 202 (2022): 117107.

Maia, Mateus, Arthur R. Azevedo, and Anderson Ara. "Predictive comparison between random machines and random forests." Journal of Data Science 19.4 (2021): 593-614.

Ara, Anderson, et al. "Random machines: A bagged-weighted support vector model with free kernel choice." Journal of Data Science 19.3 (2021): 409-428.

}
\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}
\examples{
#' library(rmachines)
#' sim_data <- rmachines::sim_class(n = 100)
#' rm_mod <- rmachines::random_machines(y~., train = sim_data, validation = sim_data)
}
