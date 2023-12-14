\name{randomMachines}
\alias{randomMachines}
\title{Random Machines}
\description{Random Machines is an ensemble model which uses the combination of different kernel functions to improve the diversity in the bagging approach, improving the predictions in general. Random Machines was developed for classification and regression problems by bagging multiple kernel functions in support vector models.

Random Machines uses SVMs (Cortes and Vapnik, 1995) as base learners in the bagging procedure with a random sample of kernel functions to build them.

Let a training sample given by \eqn{(\boldsymbol{x_{i}},y_i)} with \eqn{i=1,\dots, n} observations, where \eqn{\boldsymbol{x_{i}}} is the vector of independent variables and \eqn{y_{i}} the dependent one. The kernel bagging method initializes by training of the \eqn{r} single learner, where \eqn{r=1,\dots,R} and \eqn{R} is the total number of different kernel functions that could be used in support vector models. In this implementation the default value is \eqn{R=4} (gaussian, polynomial, laplacian and linear). See more details below.

Each single learner is internally validated and the weights \eqn{\lambda_{r}} are calculated proportionally to the strength from the single predictive performance.

Afterwards, \eqn{B} bootstrap samples are sampled from the training set. A support vector machine model \eqn{g_{b}} is trained for each bootstrap sample, \eqn{b=i,\dots,B} and the kernel function that will be used for \eqn{g_{b}} will be determined by a random choice with probability \eqn{\lambda_{r}}. The final weight \eqn{w_b} in the bagging procedure is calculated by out-of-bag samples.

The final model \eqn{G(\boldsymbol{x}_i)} for a new \eqn{\boldsymbol{x}_i} is given by,


The weights \eqn{\lambda_{r}} and \eqn{w_b} are different calculated for each task (classification, probabilistic classification and regression). See more details in the references.



\itemize{
     \item For a binary classification problem \eqn{\mathbin{{ G(\boldsymbol{x_{i}})= \text{sgn} \left( \sum_{b=1}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})\right)}}}, where \eqn{g_b} are single binary classification outputs;
     \item For a probabilistic binary classification problem \eqn{\mathbin{{ G(\boldsymbol{x_{i}})= \sum_{b=1}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})}}}, where \eqn{g_b} are single probabilistic classification outputs;
     \item For a regression problem \eqn{G(\boldsymbol{x_{i}})= \sum_{b=1}^{B}w_{b}g_{b}(\boldsymbol{x_{i}})}, , where \eqn{g_b} are single regression outputs.
}

}

\usage{
randomMachines(
     formula,
     train,validation,
     B = 25, cost = 1,
     automatic_tuning = FALSE,
     gamma_rbf = 1,
     gamma_lap = 1,
     degree = 2,
     poly_scale = 1,
     offset = 0,
     gamma_cau = 1,
     d_t = 2,
     kernels = c("rbfdot", "polydot", "laplacedot", "vanilladot"),
     prob_model = TRUE,
     loss_function = RMSE,
     epsilon = 0.1,
     beta = 2
)
}

\arguments{
     \item{formula}{
     an object of class \code{\link{formula}}: it should contain a symbolic description of the model to be fitted, indicating the dependent variable and all predictors that should be included.
     }
     \item{train}{
     the training data \eqn{\left\{\left( \mathbf{x}_{i},y_{i} \right)\right\}_{i=1}^{n}} used to train the model.
     }
     \item{validation}{
     the validation data \eqn{\left\{\left( \mathbf{x}_{i},y_{i}\right)  \right\}_{i=1}^{V}} used to calculate probabilities \eqn{\lambda_{r}}. If \code{validation = NULL},the validation set is going be selected as 0.25 partition from the training data, and the remaining partition is selected as the new training sample.
     }
     \item{B}{
     number of bootstrap samples. The default value is \code{B=25}.
     }
     \item{cost}{
     the \eqn{C}-constant term of the regularization on soft margins at support vector models. The default value is \code{cost=1}.
     }
     \item{automatic_tuning}{
     boolean to define if the kernel hyperparameters will be selected using the \code{sigest} from the \code{ksvm} function. The default value is \code{FALSE}.
     }
     \item{gamma_rbf}{
     the hyperparameter \eqn{\gamma_{g}} used in the RBF kernel. The default value is \code{gamma_rbf=1}.
     }
     \item{gamma_lap}{
     the hyperparameter \eqn{\gamma_{l}} used in the Laplacian kernel. The default value is \code{gamma_lap=1}.
     }
     \item{degree}{
     the degree used in the Polynomial kernel. The default value is \code{degree=2}.
     }
     \item{poly_scale}{
     the scale parameter from the Polynomial kernel. The default value is \code{poly_scale=1}.
     }
     \item{offset}{
     the offset parameter from the Polynomial kernel. The default value is \code{offset=0}.
     }
     \item{gamma_cau}{
     the hyperparameter \eqn{\gamma_{c}} used in the Cauchy kernel. The default value is \code{gamma_cau=1}.
     }
     \item{d_t}{
     the \eqn{d_{t}}-norm from the t-Student kernel. The default value is \code{d_t=2}.
     }
     \item{kernels}{
     a vector with the name of kernel functions that will be used in the Random Machines model. The default include the kernel functions: \code{c("rbfdot", "polydot", "laplacedot", "vanilladot").} The other kernel functions as \code{"cauchydot"} and \code{"tdot"} are exclusive to the binary classification setting.
     }
     \item{prob_model}{
     a boolean to define if the algorithm will be using a probabilistic approach to the define the predictions (default = \code{TRUE}).
     }
     \item{loss_function}{
     Define which loss function is going to be used in the regression approach. The default is the \code{RMSE} function but others can be added.
     }
     \item{epsilon}{
     The epsilon in the loss function used from the SVR implementation. The default value is \code{epsilon=0.1}.
     }
     \item{beta}{
     The correlation parameter \eqn{\beta} which calibrates the penalisation of each kernel performance in regression tasks.  The default value is \code{beta=2}.
     }
}
\value{ \code{randomMachines()} returns an object of \code{\link{class}} "rm_class" for classification tasks or "rm_reg" for if the target variable is a continuous numerical response. See \code{\link[=predict.rm_class]{predict.rm_class}} or \code{\link{predict.rm_reg}} for more details of how to obtain predictions from each model respectively.}
\details{
     The Random Machines is an ensemble method which combines the bagging procedure proposed by Breiman (1996), using Support Vector Machine models as base learners jointly with a random selection of kernel functions that add diversity to the ensemble without harming its predictive performance. The kernel functions \eqn{k(x,y)} are described by the functions below,

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

Ara, Anderson, et al. "Random machines: A bagged-weighted support vector model with free kernel choice." Journal of Data Science 19.3 (2021): 409-428.

Breiman, L. (1996). Bagging predictors. Machine learning, 24, 123-140.

Cortes, C., and Vapnik, V. (1995). Support-vector networks. Machine learning, 20, 273-297.

Maia, Mateus, Arthur R. Azevedo, and Anderson Ara. "Predictive comparison between random machines and random forests." Journal of Data Science 19.4 (2021): 593-614.





}
\author{
Mateus Maia: \email{mateusmaia11@gmail.com},
Gabriel Felipe Ribeiro: \email{brielribeiro08@gmail.com},
Anderson Ara: \email{ara@ufpr.br}
}
\examples{
library(randomMachines)

# Simulation from a binary output context
sim_data <- sim_class(n = 75)

## Setting the training and validation set
sim_new <- sim_class(n = 75)

# Modelling Random Machines (probabilistic output)
rm_mod_prob <- randomMachines(y~., train = sim_data)

## Modelling Random Machines (binary class output)
rm_mod_label <- randomMachines(y~., train = sim_data,prob_model = FALSE)

## Predicting for new data
y_hat <- predict(rm_mod_label,sim_new)
}
