\name{random_machines}
\alias{random_machines}
\title{Random Machines}
\description{Random Machines is an ensemble model which use the combination of different kernel functions to improve the diversity in the bagging approach improving the predictions. Random Machines was developed for classification and regression problems by bagging multiple kernel functions in support vector models.

\itemize{
     \item For a binary classification problem \eqn{\mathbin{{ G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right), }}}
     \item For a numeric response \eqn{G(\boldsymbol{x_{i}})= \sum_{b}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})}
}

}

\usage{
random_machines(
     formula,
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
     kernels = c("rbfdot", "polydot", "laplacedot", "vanilladot"),
     prob_model = T,
     loss_function = RMSE,
     epsilon = 0.1,
     beta = 2
)
}

\arguments{
     \item{formula}{
     formula an object of class \code{\link{formula}}: it should contain a symbolic description of the model to be fitted, indicating the dependent variable and all predictors that should be included.
     }
     \item{train}{
     the training data \eqn{\left\{\left( \mathbf{x}_{i},y_{i} \right)\right\}_{i=1}^{N}} used to train the model.
     }
     \item{validation}{
     the validation data \eqn{\left\{\left( \mathbf{x}_{i},y_{i}\right)  \right\}_{i=1}^{V}} used to calculate probabilities \eqn{\lambda_{r}}. If \code{validation = NULL},the validation set is going be selected as 0.25 partition from the training data, and the remaining partition is selected as the new training sample.
     }
     \item{B}{
     number of bootstrap samples.
     }
     \item{cost}{
     the \eqn{C}-constant term of the regularization on the Lagrange formulation.
     }
     \item{seed.bootstrap}{
     setting a seed to replicate bootstrap sampling. The default value is \code{NULL}.
     }
     \item{automatic_tuning}{
     boolean to define if the kernel hyperparameters will be selected using the \code{sigest} from the \code{ksvm} function.
     }
     \item{gamma_rbf}{
     the hyperparameter \eqn{\gamma_{g}} used in the RBF kernel.
     }
     \item{gamma_lap}{
     the hyperparameter \eqn{\gamma_{l}} used in the Laplacian kernel.
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
     the hyperparameter \eqn{\gamma_{c}} used in the Cauchy kernel.
     }
     \item{d_t}{
     the \eqn{d_{t}}-norm from the Sudent's kernel.
     }
     \item{kernels}{
     a vector with the name of kernel functions that will be used in the Random Machines model. The default include the kernel functions: \code{c("rbfdot", "polydot", "laplacedot", "vanilladot").} The other kernel functions as \code{"cauchydot"} and \code{"tdot"} are exclusive to the binary classification setting, being prebuilt in the package.
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
\value{ \code{random_machines()} returns an object of \code{\link{class}} "rm_class" for classification tasks or "rm_reg" for if the target variable is a continuous numerical response. See \code{\link[=predict.rm_class]{predict.rm_class}} or \code{\link{predict.rm_reg}} for more details of how obtain predictions from each model respectively.}
\details{
     The Random Machines is an ensemble method which combines the bagging procedure proposed by Breiman (1996), using Support Vector Machine models as base learners jointly with a random selection of kernel functions that add diversity to the ensemble without harm its predictive performance. The kernel functions \eqn{k(x,y)} are described by the functions below,

\itemize{
\item {Linear Kernel}: \eqn{k(x,y) = (x\cdot y)}
\item {Polynomial Kernel}: \eqn{k(x,y) = \left(scale( x\cdot y) + offset\right)^{degree}}
\item {Gaussian Kernel}: \eqn{k(x,y) = e^{-\gamma_{g}||x-y||^2}}
\item {Laplacian Kernel}:  \eqn{k(x,y) = e^{-\gamma_{\ell}||x-y||}}
\item {Cauchy Kernel}: \eqn{k(x,y) = \frac{1}{1 + \frac{||x-y||^{2}}{\gamma_{c}}}}
\item {Student's t Kernel}: \eqn{k(x,y) = \frac{1}{1 + ||x-y||^{d_{t}}}}
}

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
#'
#' # Simulation from a binary output context
#' sim_data <- rmachines::sim_class(n = 100)
#' # Setting the training and validation set
#' sim_new <- rmachines::sim_class(n = 100)
#'
#' # Modelling Random Machines (probabilistic output)
#' rm_mod_prob <- rmachines::random_machines(y~., train = sim_data,)
#' # Modelling Random Machines (class output)
#' rm_mod_label <-rmachines::::random_machines(y~., train = sim_data,prob_model = FALSE)
#' # Predicting for new data
#' y_hat <- predict(rm_mod,sim_new)
}
