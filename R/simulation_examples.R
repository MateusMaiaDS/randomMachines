#' Function to generate a two classification dataset from normal distribution
#'
#' This function predicts the outcome for a RM object model using new data
#'
#' @param n Sample size
#' @param p Number of predictors
#' @param ratio Ratio between class A and class B
#' @param mean_one Mean of X1
#' @param sd_one Standard deviation of X1
#' @param mean_two Mean of X2
#' @param sd_two Standard devation of X2
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
     x_a <- replicate(p,rnorm(n_a,mean = mean_one,sd = sd_one))
     colnames(x_a) <- paste("x",1:p)

     x_b <- replicate(p,rnorm(n_b,mean = mean_two,sd = sd_two))
     colnames(x_b) <- paste("x",1:p)

     # Formating the complete dataset
     x <- rbind(x_a,x_b)
     y <- as.factor(c(rep("A",n_a),rep("B",n_b)))

     simulated_data <- data.frame(x,y)

     return(simulated_data[sample(nrow(simulated_data)),])
}
