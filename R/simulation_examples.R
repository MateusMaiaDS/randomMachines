#'
#' @export
sim_class <- function(n, p = 2 ,ratio = 0.5 , mu_a = 0,
                      sigma_a = 1,mu_b = 1,sigma_b = 1){

     # Setting the number of observations from the first data set
     n_a <- round(n*abs(1-ratio))
     n_b <- round(n*ratio)

     # Generating values from the X observations
     x_a <- replicate(p,stats::rnorm(n_a,mean = mu_a,sd = sigma_a))
     colnames(x_a) <- paste("x",1:p)

     x_b <- replicate(p,stats::rnorm(n_b,mean = mu_b,sd = sigma_b))
     colnames(x_b) <- paste("x",1:p)

     # Formating the complete dataset
     x <- rbind(x_a,x_b)
     y <- as.factor(c(rep("A",n_a),rep("B",n_b)))

     simulated_data <- data.frame(x,y)

     return(simulated_data[sample(nrow(simulated_data)),])
}


#' @export
sim_reg1 <- function(n,sigma = 0.5){


        # Generating the x
        x <- replicate(2,stats::runif(n,min = -1,max = 1))

        colnames(x) <- paste0("x.",1:2)

        # Generating the y
        y  <- x[,1]^2 + exp(-x[,2]^2) + stats::rnorm(n = n,mean = 0,sd = sigma)

        return(data.frame(x,y=y))
}

#' @export
sim_reg2 <- function(n,sigma = 0.5){


        # Generating the x
        x <- replicate(8,stats::runif(n,min = -1,max = 1))

        colnames(x) <- paste0("x.",1:8)

        # Generating the y
        y  <- x[,1]*x[,2] + x[,3]^2 -x[,4]*x[,7] +x[,5]*x[,8] -x[,6]^2 + stats::rnorm(n = n,mean = 0,sd = sigma)

        return(data.frame(x,y=y))
}

#' @export
sim_reg3 <- function(n,sigma = 0.5){

        # Generating the x
        x <- replicate(4,stats::runif(n,min = -1,max = 1))

        colnames(x) <- paste0("x.",1:4)

        # Generating the y
        y  <- -sin(x[,1]) + x[,4]^2 + x[,3] - exp(-x[,4]^2) + stats::rnorm(n = n,mean = 0,sd = sqrt(0.5))

        return(data.frame(x,y=y))
}

#' @export
sim_reg4 <- function(n,sigma = 0.5){

        # Generating the x
        x <- replicate(4,stats::runif(n,min = -1,max = 1))

        colnames(x) <- paste0("x.",1:4)

        # Generating the y
        y  <- -sin(x[,1]) + x[,4]^2 + x[,3] - exp(-x[,4]^2) + stats::rnorm(n = n,mean = 0,sd = sqrt(0.5))

        return(data.frame(x,y=y))
}


#' @export
sim_reg5 <- function(n,sigma = 0.5){


        # Generating the x
        x <- replicate(6,stats::rnorm(n = n))

        colnames(x) <- paste0("x.",1:6)

        # Generating the y
        y  <- x[,1] + 0.707*x[,2]^2 + 2*ifelse(x[,3]>0,1,0) + 0.873*log(abs(x[,1]))*abs(x[,3]) + 0.894*x[,2]*x[,4] +
                2*ifelse(x[,5]>0,1,0) + 0.464*exp(x[,6]) + stats::rnorm(n = n,mean = 0,sd = sigma)

        return(data.frame(x,y=y))
}
