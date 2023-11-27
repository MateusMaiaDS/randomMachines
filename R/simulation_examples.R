#' Function to generate a two classification data set from normal distribution
#'
#' This function predicts the outcome for a RM object model using new data
#'
#' @param n Sample size
#' @param p Number of predictors
#' @param ratio Ratio between class A and class B
#' @param mean_one Mean of \eqn{X_{1}}.
#' @param sd_one Standard deviation of \eqn{X_{1}}.
#' @param mean_two Mean of \eqn{X_{2}}
#' @param sd_two Standard devation of \eqn{X_{2}}
#' @param seed Setting a seed for reproducibility of results. The default is NULL
#' @return A simulated data.frame with two predictors for a binary classification problem
#'
#' @examples
#' library(rmachines)
#' sim_data <- sim_class(n = 100)
#' @export
sim_class <- function(n, p = 2 ,ratio = 0.5 , mean_one = 0,
                      sd_one = 1,mean_two = 1,sd_two = 1,
                      seed = NULL){
     # Setting the seed
     set.seed(seed)

     # Setting the number of observations from the first data set
     n_a <- round(n*abs(1-ratio))
     n_b <- round(n*ratio)

     # Generating values from the X observations
     x_a <- replicate(p,stats::rnorm(n_a,mean = mean_one,sd = sd_one))
     colnames(x_a) <- paste("x",1:p)

     x_b <- replicate(p,stats::rnorm(n_b,mean = mean_two,sd = sd_two))
     colnames(x_b) <- paste("x",1:p)

     # Formating the complete dataset
     x <- rbind(x_a,x_b)
     y <- as.factor(c(rep("A",n_a),rep("B",n_b)))

     simulated_data <- data.frame(x,y)

     return(simulated_data[sample(nrow(simulated_data)),])
}


#' Simple regression \eqn{y = x_{1}^{2} + e^{x_{2}^{2}} + \varepsilon} case based on Ara et. al 2022
#'
#' @param n Sample size
#' @param seed Define a seed to run the simulation. NULL is the default
#'
#' @return A simulated data.frame with two predictors and the target variable.
#' @export
#'
#' @examples
#' library(rmachines)
#' sim_data <- sim_reg(n=100)
sim_reg <- function(n,seed= NULL){

        # Setting the seed.
        set.seed(seed)

        # Generating the x
        x <- replicate(2,stats::runif(n,min = -1,max = 1))

        colnames(x) <- paste0("x.",1:2)

        # Generating the y
        y  <- x[,1]^2 + exp(-x[,2]^2) + stats::rnorm(n = n,mean = 0,sd = sqrt(0.1))

        return(data.frame(x,y=y))
}


