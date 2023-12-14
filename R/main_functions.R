#' @importMethodsFrom kernlab predict
#' @export
randomMachines <- function(formula,
                                train,
                                validation = NULL,
                                B = 25,
                                cost = 1,
                                automatic_tuning = FALSE,
                                gamma_rbf = 1,
                                gamma_lap = 1,
                                degree = 2,
                                poly_scale = 1,
                                offset = 0,
                                gamma_cau = 1, d_t = 2,
                                kernels = c("rbfdot", "polydot", "laplacedot",
                                            "vanilladot"),
                                prob_model = TRUE,
                                # Default parameters for the regression case
                                loss_function = RMSE,
                                epsilon = 0.1,
                                beta = 2
                                ) {

   # Selecting a validation set
   if(is.null(validation)){
     train_index <- sample(1:nrow(train), round(0.75*nrow(train)))
     validation <- train[train_index,,drop = FALSE]
     train <- train[-train_index,,drop = FALSE]
   }

   # Checking the class of the training data
   if(is.numeric(train[[formula[[2]]]])){
     reg_rm <- TRUE
   } else {
     reg_rm <- FALSE
   }

   # Preventing
   if(reg_rm){
     valid_kernels<- c("rbfdot", "polydot", "laplacedot","vanilladot")
     if(any(!(kernels %in% valid_kernels))){
       stop("Insert only valid kernel functions.")
     }

   } else {
     valid_kernels<- c("rbfdot", "polydot", "laplacedot","vanilladot", "cauchydot","tdot")
     if(any(!(kernels %in% valid_kernels))){
       stop("Insert only valid kernel functions.")
     }
   }

   # Selecting
   if(prob_model & (!reg_rm)){
      rm_mod <- random_machines_prob(formula = formula,
                           train = train,
                           validation = validation,
                           B = B,
                           cost = cost,
                           automatic_tuning = automatic_tuning,
                           gamma_rbf = gamma_rbf,
                           gamma_lap = gamma_lap,
                           degree = degree,
                           poly_scale = poly_scale,
                           offset = offset,
                           gamma_cau = gamma_cau,
                           d_t = d_t,
                           kernels = kernels)
   } else if((!prob_model) & (!reg_rm)){
     rm_mod <- random_machines_acc(formula = formula,
                                    train = train,
                                    validation = validation,
                                    B = B,
                                    cost = cost,
                                    automatic_tuning = automatic_tuning,
                                    gamma_rbf = gamma_rbf,
                                    gamma_lap = gamma_lap,
                                    degree = degree,
                                    poly_scale = poly_scale,
                                    offset = offset,
                                    gamma_cau = gamma_cau,
                                    d_t = d_t,
                                    kernels = kernels)
   } else if(reg_rm) {
     if(length(kernels)!=4){
       warning("The regression version of Random Machines uses Linear, Polynomial, Gaussian and Laplacian kernels.")
     }
     rm_mod <- regression_random_machines(formula = formula,
                                          train = train,
                                          validation = validation,
                                          B = B,
                                          cost = cost,
                                          automatic_tuning = automatic_tuning,
                                          gamma_rbf = gamma_rbf,
                                          gamma_lap = gamma_lap,
                                          degree = degree,
                                          poly_scale = poly_scale,
                                          offset = offset,
                                          epsilon = epsilon,beta = beta,loss_function = loss_function)
   } else {
     stop("Insert a valid data.frame() with a valid response.")
   }

  return(rm_mod)
}

#' @importMethodsFrom kernlab predict
random_machines_prob <- function(formula,
                                 train,
                                 validation,
                                 B = 25,
                                 cost = 1,
                                 automatic_tuning = FALSE,
                                 gamma_rbf = 1,
                                 gamma_lap = 1,
                                 degree = 2,
                                 poly_scale = 1,
                                 offset = 0,
                                 gamma_cau = 1, d_t = 2,
                                 kernels = c("rbfdot", "polydot", "laplacedot",
                                             "vanilladot", "cauchydot")) {

  # New kernel functions used by the new article.
  cauchydot <- function(sigma = 1) {
    norma <- function(x, y) {
      return(norm(matrix(x - y),type = "2"))
    }

    rval <- function(x, y = NULL) {
      if (!is.vector(x)) {
        stop("x must be a vector")
      }
      if (!is.vector(y) && !is.null(y)) {
        stop("y must a vector")
      }
      if (is.vector(x) && is.null(y)) {
        return(1)
      }
      if (is.vector(x) && is.vector(y)) {
        if (!length(x) == length(y)) {
          stop("number of dimension must be the same on both data points")
        }
        return(1 / (1 + (norma(x, y) / sigma)))
      }
    }
    return(methods::new("kernel", .Data = rval, kpar = list(sigma = sigma)))
  }


  tdot <- function(d = 2) {
    norma <- function(x, y, d) {
      return(sqrt(norm(matrix(x - y),type = "2"))^d)
    }

    rval <- function(x, y = NULL) {
      if (!is.vector(x)) {
        stop("x must be a vector")
      }
      if (!is.vector(y) && !is.null(y)) {
        stop("y must a vector")
      }
      if (is.vector(x) && is.null(y)) {
        return(1)
      }
      if (is.vector(x) && is.vector(y)) {
        if (!length(x) == length(y)) {
          stop("number of dimension must be the same on both data points")
        }
        return(1 / (1 + norma(x, y, d)))
      }
    }
    return(methods::new("kernel", .Data = rval, kpar = list(d = d)))
  }

  # Calculating the probabilities
  test <- validation
  class_name <- as.character(formula[[2]])
  prob_weights <- list()
  kernel_type <- kernels

  # Getting the levels of the target variable
  y_levels <- levels(train[[class_name]])
  # Checking wether or not use automatic tuning
  if (automatic_tuning) {
    early_model <- lapply(kernel_type, function(k_type){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = train, type = "C-svc",
                    kernel = if (k_type == "vanilladot") {
                      "polydot"
                    } else if (k_type == "cauchydot") {
                      cauchydot()
                    } else if (k_type == "tdot") {
                      tdot()
                    } else {
                      k_type
                    }, C = cost, kpar = if (k_type == "laplacedot" || k_type == "rbfdot") {
                      "automatic"
                    } else if (k_type == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (k_type == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (k_type == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  } else {
    early_model <- lapply(kernel_type, function(k_type){

      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = train, type = "C-svc",
                    kernel = if (k_type == "vanilladot") {
                      "polydot"
                    } else if (k_type == "cauchydot") {
                      cauchydot()
                    } else if (k_type == "tdot") {
                      tdot()
                    } else {
                      k_type
                    }, C = cost, kpar = if (k_type == "laplacedot") {
                      list(sigma = gamma_lap)
                    } else if (k_type == "rbfdot") {
                      list(sigma = gamma_rbf)
                    } else if (k_type == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (k_type == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (k_type == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  }

  # Getting tprobabilities
  predict <- lapply(early_model,function(y){predict(y,newdata = test, type = "probabilities")[, 2]})


  # CONTINUE FROM HERE
  m_brier <- unlist(lapply(predict, function(p) {brier_score(prob = p, observed = (test[[class_name]]),levels = y_levels)}))

  log_brier <- log((1 - m_brier) / m_brier)
  log_brier[is.infinite(log_brier)] <- 1 # Sometimes the brier can be equal to 0, so this line certify to not produce any NA
  prob_weights <- log_brier / sum(log_brier)
  prob_weights <- ifelse(prob_weights < 0, 0, prob_weights) # Avoiding negative values of probabilities

  models <- rep(list(0), B)
  boots_sample <- list(rep(B))
  out_of_bag <- list(rep(B))
  boots_index_row <- rep(list(nrow(train)) ,B)
  at_least_one <- NULL

  # Sampling the bootstrap samples
  while (is.null(at_least_one)) {
    boots_index_row_new <- lapply(
      boots_index_row,function(x){
        sample(1:x, x, replace = TRUE)
      }
    )
    boots_sample <- lapply(boots_index_row_new, function(x) {train[x, ]})
    out_of_bag <- lapply(boots_index_row_new, function(x) {train[-unique(x), ]})

    ## The class 1 needs to be contained at least once so you have a valid bootstrap sample
    for (p in 1:length(boots_sample)) {
      while (table(boots_sample[[p]][class_name])[2] < 2) {
        boots_index_row_new_new <- lapply(
          boots_index_row, function(x){
            sample(1:x, x, replace = TRUE)
          }
        )
        boots_sample_new <- lapply(boots_index_row_new_new, function(x){train[x, ]})
        out_of_bag_new <- lapply(boots_index_row_new_new, function(x){train[-unique(x), ]})

        boots_sample[[p]] <- boots_sample_new[[1]]
        out_of_bag[[p]] <- out_of_bag_new[[1]]
      }
    }
    ## Checking if any has length 0
    if (any(unlist(lapply(boots_sample, function(x) { table(x[[class_name]]) == 0 })))) {
      at_least_one <- NULL
    } else {
      at_least_one <- 1
    }

  }


  # Sampling different kernel functions
  random_kernel <- sample(kernel_type, B, replace = TRUE, prob = prob_weights)

  if (automatic_tuning) {
    models <- mapply(boots_sample, random_kernel,FUN = function(boot_sample,rand_kern){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = boot_sample, type = "C-svc", kernel = if (rand_kern == "vanilladot") {
                      "polydot"
                    } else if (rand_kern == "cauchydot") {
                      cauchydot()
                    } else if (rand_kern == "tdot") {
                      tdot()
                    } else {
                      rand_kern
                    }, C = cost, kpar = if (rand_kern == "laplacedot" || rand_kern ==
                                            "rbfdot") {
                      "automatic"
                    } else if (rand_kern == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (rand_kern == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (rand_kern == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  } else {
    models <- mapply(boots_sample, random_kernel,FUN = function(boot_sample,rand_kern){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = boot_sample, type = "C-svc", kernel = if (rand_kern == "vanilladot") {
                      "polydot"
                    } else if (rand_kern == "cauchydot") {
                      cauchydot()
                    } else if (rand_kern == "tdot") {
                      tdot()
                    } else {
                      rand_kern
                    }, C = cost, kpar = if (rand_kern == "laplacedot") {
                      list(sigma = gamma_lap)
                    } else if (rand_kern == "rbfdot") {
                      list(sigma = gamma_rbf)
                    } else if (rand_kern == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (rand_kern == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (rand_kern == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  }

  # Predicting for each model and getting the weights
  predict <- lapply(models, function(mod){predict(mod, newdata = test, type = "probabilities")[, 2]})
  predict_oobg <- mapply(models, out_of_bag,FUN = function(mod,oob){
    predict(mod,
                     newdata = oob,
                     type = "probabilities")[, 2]})

  kernel_weight_raw <- unlist(mapply(predict_oobg, out_of_bag, FUN = function(pred_oob,oob){
    brier_score(pred_oob, oob[[class_name]], levels = levels(oob[[class_name]]))}))

  kernel_weight <- 1 / kernel_weight_raw^2

  # Re-creating the vector of names
  kern_names_final <- kernel_type
  for(i in 1:length(kernel_type)){
    kern_names_final[i] <- if(kernel_type[i] == "rbfdot" ) {
      "RBF_Kern"
    } else if (kernel_type[i] =="polydot"){
      "POL_Kern"
    } else if (kernel_type[i] == "laplacedot"){
      "LAP_Kern"
    } else if(kernel_type[i] == "vanilladot"){
      "LIN_Kern"
    } else if(kernel_type[i] == "cauchydot"){
      "CAU_Kern"
    } else if(kernel_type[i] == "tdot") {
      "T_Kern"
    }
  }

  if (length(kern_names_final) == 2) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = TRUE
    )
  } else if (length(kern_names_final) == 3) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = TRUE
    )
  } else if (length(kern_names_final) == 4) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = TRUE
    )
  } else if (length(kern_names_final) == 5) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4],
        prob_weights[5]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = TRUE

    )
  } else if (length(kern_names_final) == 6) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4],
        prob_weights[5],
        prob_weights[6]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = TRUE
    )
  } else {
    stop("The number of kernel isn't compatible")
  }

  attr(model_result, "class") <- "rm_class"
  return(model_result)
}


random_machines_acc <- function(formula,
                                 train,
                                 validation,
                                 B = 25,
                                 cost = 1,
                                 automatic_tuning = FALSE,
                                 gamma_rbf = 1,
                                 gamma_lap = 1,
                                 degree = 2,
                                 poly_scale = 1,
                                 offset = 0,
                                 gamma_cau = 1, d_t = 2,
                                 kernels = c("rbfdot", "polydot", "laplacedot",
                                             "vanilladot", "cauchydot")) {

  # New kernel functions used by the new article.
  cauchydot <- function(sigma = 1) {
    norma <- function(x, y) {
      return(norm(matrix(x - y),type = "2"))
    }

    rval <- function(x, y = NULL) {
      if (!is.vector(x)) {
        stop("x must be a vector")
      }
      if (!is.vector(y) && !is.null(y)) {
        stop("y must a vector")
      }
      if (is.vector(x) && is.null(y)) {
        return(1)
      }
      if (is.vector(x) && is.vector(y)) {
        if (!length(x) == length(y)) {
          stop("number of dimension must be the same on both data points")
        }
        return(1 / (1 + (norma(x, y) / sigma)))
      }
    }
    return(methods::new("kernel", .Data = rval, kpar = list(sigma = sigma)))
  }


  tdot <- function(d = 2) {
    norma <- function(x, y, d) {
      return(sqrt(norm(matrix(x - y),type = "2"))^d)
    }

    rval <- function(x, y = NULL) {
      if (!is.vector(x)) {
        stop("x must be a vector")
      }
      if (!is.vector(y) && !is.null(y)) {
        stop("y must a vector")
      }
      if (is.vector(x) && is.null(y)) {
        return(1)
      }
      if (is.vector(x) && is.vector(y)) {
        if (!length(x) == length(y)) {
          stop("number of dimension must be the same on both data points")
        }
        return(1 / (1 + norma(x, y, d)))
      }
    }
    return(methods::new("kernel", .Data = rval, kpar = list(d = d)))
  }

  # Calculating the probabilities
  test <- validation
  class_name <- as.character(formula[[2]])
  prob_weights <- list()
  kernel_type <- kernels

  # Getting the levels of the target variable
  y_levels <- levels(train[[class_name]])
  # Checking wether or not use automatic tuning
  if (automatic_tuning) {
    early_model <- lapply(kernel_type, function(k_type){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = train, type = "C-svc",
                    kernel = if (k_type == "vanilladot") {
                      "polydot"
                    } else if (k_type == "cauchydot") {
                      cauchydot()
                    } else if (k_type == "tdot") {
                      tdot()
                    } else {
                      k_type
                    }, C = cost, kpar = if (k_type == "laplacedot" || k_type == "rbfdot") {
                      "automatic"
                    } else if (k_type == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (k_type == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (k_type == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  } else {
    early_model <- lapply(kernel_type, function(k_type){

      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = train, type = "C-svc",
                    kernel = if (k_type == "vanilladot") {
                      "polydot"
                    } else if (k_type == "cauchydot") {
                      cauchydot()
                    } else if (k_type == "tdot") {
                      tdot()
                    } else {
                      k_type
                    }, C = cost, kpar = if (k_type == "laplacedot") {
                      list(sigma = gamma_lap)
                    } else if (k_type == "rbfdot") {
                      list(sigma = gamma_rbf)
                    } else if (k_type == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (k_type == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (k_type == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  }

  # Getting tprobabilities
  predict <- lapply(early_model,function(y){predict(y,newdata = test)})

  # Create a function to calculate acc
  acc <- function(observed,pred) {
    # Getting the table
    return(sum(diag(table(observed,pred)))/sum(table(observed,pred)))
  }

  # CONTINUE FROM HERE
  acc_kern <- unlist(lapply(predict, function(p) {acc(pred= p, observed = (test[[class_name]]))}))

  log_acc <- log(acc_kern / (1 - acc_kern))
  log_acc[is.infinite(log_acc)] <- 1
  prob_weights <- log_acc / sum(log_acc)
  prob_weights <- ifelse(prob_weights < 0, 0, prob_weights)

  models <- rep(list(0), B)
  boots_sample <- list(rep(B))
  out_of_bag <- list(rep(B))
  boots_index_row <- rep(list(nrow(train)) ,B)
  at_least_one <- NULL

  # Sampling the bootstrap samples
  while (is.null(at_least_one)) {
    boots_index_row_new <- lapply(
      boots_index_row,function(x){
        sample(1:x, x, replace = TRUE)
      }
    )
    boots_sample <- lapply(boots_index_row_new, function(x) {train[x, ]})
    out_of_bag <- lapply(boots_index_row_new, function(x) {train[-unique(x), ]})

    ## The class 1 needs to be contained at least once so you have a valid bootstrap sample
    for (p in 1:length(boots_sample)) {
      while (table(boots_sample[[p]][class_name])[2] < 2) {
        boots_index_row_new_new <- lapply(
          boots_index_row, function(x){
            sample(1:x, x, replace = TRUE)
          }
        )
        boots_sample_new <- lapply(boots_index_row_new_new, function(x){train[x, ]})
        out_of_bag_new <- lapply(boots_index_row_new_new, function(x){train[-unique(x), ]})

        boots_sample[[p]] <- boots_sample_new[[1]]
        out_of_bag[[p]] <- out_of_bag_new[[1]]
      }
    }
    ## Checking if any has length 0
    if (any(unlist(lapply(boots_sample, function(x) { table(x[[class_name]]) == 0 })))) {
      at_least_one <- NULL
    } else {
      at_least_one <- 1
    }

  }


  # Sampling different kernel functions
  random_kernel <- sample(kernel_type, B, replace = TRUE, prob = prob_weights)

  if (automatic_tuning) {
    models <- mapply(boots_sample, random_kernel,FUN = function(boot_sample,rand_kern){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = boot_sample, type = "C-svc", kernel = if (rand_kern == "vanilladot") {
                      "polydot"
                    } else if (rand_kern == "cauchydot") {
                      cauchydot()
                    } else if (rand_kern == "tdot") {
                      tdot()
                    } else {
                      rand_kern
                    }, C = cost, kpar = if (rand_kern == "laplacedot" || rand_kern ==
                                            "rbfdot") {
                      "automatic"
                    } else if (rand_kern == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (rand_kern == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (rand_kern == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  } else {
    models <- mapply(boots_sample, random_kernel,FUN = function(boot_sample,rand_kern){
      kernlab::ksvm(formula,
                    prob.model = TRUE,
                    data = boot_sample, type = "C-svc", kernel = if (rand_kern == "vanilladot") {
                      "polydot"
                    } else if (rand_kern == "cauchydot") {
                      cauchydot()
                    } else if (rand_kern == "tdot") {
                      tdot()
                    } else {
                      rand_kern
                    }, C = cost, kpar = if (rand_kern == "laplacedot") {
                      list(sigma = gamma_lap)
                    } else if (rand_kern == "rbfdot") {
                      list(sigma = gamma_rbf)
                    } else if (rand_kern == "polydot") {
                      list(degree = 2, scale = poly_scale, offset = 0)
                    } else if (rand_kern == "cauchydot") {
                      list(sigma = gamma_cau)
                    } else if (rand_kern == "tdot") {
                      list(d = d_t)
                    } else {
                      list(degree = 1, scale = poly_scale, offset = 0)
                    }
      )})
  }

  # Predicting for each model and getting the weights
  predict <- lapply(models, function(mod){predict(mod, newdata = test, type = "probabilities")[, 2]})
  predict_oobg <- mapply(models, out_of_bag,FUN = function(mod,oob){
    predict(mod,
                     newdata = oob)})

  kernel_weight_raw <- unlist(mapply(predict_oobg, out_of_bag, FUN = function(pred_oob,oob){
    acc(pred = pred_oob, oob[[class_name]])}))

  kernel_weight <- 1 / kernel_weight_raw^2

  # Re-creating the vector of names
  kern_names_final <- kernel_type
  for(i in 1:length(kernel_type)){
    kern_names_final[i] <- if(kernel_type[i] == "rbfdot" ) {
      "RBF_Kern"
    } else if (kernel_type[i] =="polydot"){
      "POL_Kern"
    } else if (kernel_type[i] == "laplacedot"){
      "LAP_Kern"
    } else if(kernel_type[i] == "vanilladot"){
      "LIN_Kern"
    } else if(kernel_type[i] == "cauchydot"){
      "CAU_Kern"
    } else if(kernel_type[i] == "tdot") {
      "T_Kern"
    }
  }

  if (length(kern_names_final) == 2) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = FALSE

    )
  } else if (length(kern_names_final) == 3) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = FALSE

    )
  } else if (length(kern_names_final) == 4) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = FALSE
    )
  } else if (length(kern_names_final) == 5) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4],
        prob_weights[5]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = FALSE
    )
  } else if (length(kern_names_final) == 6) {
    model_result <- list(
      train = train, class_name = class_name,
      kernel_weight = kernel_weight, lambda_values = stats::setNames(list(
        prob_weights[1],
        prob_weights[2],
        prob_weights[3],
        prob_weights[4],
        prob_weights[5],
        prob_weights[6]
      ), kern_names_final), model_params = list(
        class_name = class_name,
        B = B, cost = cost, gamma_rbf = gamma_rbf,
        gamma_lap = gamma_lap, degree = degree
      ), bootstrap_models = models,
      bootstrap_samples = boots_sample,
      prob_model = FALSE

    )
  } else {
    stop("The number of kernel isn't compatible")
  }

  attr(model_result, "class") <- "rm_class"
  return(model_result)
}


regression_random_machines<-function(formula,#Formula that will be used
                                     train,#The Training set
                                     validation,#The validation set
                                     B=25, #B correspoding to the number of bootstrap samples
                                     cost=1,#Cost parameter of SVM
                                     gamma_rbf=1,#Gamma used in Table 1.
                                     gamma_lap=1,
                                     poly_scale = 1, # Scale factor from polynomial kernel
                                     offset = 0, # offset from the polynomial kernel
                                     degree=2,#Degree used in Table 1.
                                     epsilon=0.1,beta=2,
                                     loss_function,automatic_tuning=FALSE #Choose a loss-fucntion

){

  # Creating the class name variable
  class_name <- as.character(formula[[2]])

  #Probability associated with each kernel function

  prob_weights<-list()

  #The Kernel types used in the algorithm
  kernel_type<-c('rbfdot','polydot','laplacedot','vanilladot')

  #TUNING AUTOMÁTICO
  if(automatic_tuning){

    early_model<- lapply(kernel_type,function(kern_type){kernlab::ksvm(formula,data=train,type="eps-svr",
                                        kernel=if(kern_type=="vanilladot"){
                                          "polydot"
                                        }else{
                                          kern_type
                                        },
                                        C=cost,
                                        kpar=if(kern_type=='laplacedot' ||kern_type=='rbfdot')
                                        {
                                          "automatic"
                                        }else if(kern_type=='polydot'){
                                          list(degree=2,scale=poly_scale,offset=offset)
                                        }else{
                                          list(degree=1,scale=poly_scale,offset=offset)
                                        },
                                        epsilon=epsilon)})
  }else{
    #The early model that will calculate the probabilities that will be used during the sort process
    early_model<- lapply(kernel_type,function(kern_type){kernlab::ksvm(formula,data=train,type="eps-svr",
                                       kernel=if(kern_type=="vanilladot"){
                                         "polydot"
                                       }else{
                                         kern_type
                                       },
                                       C=cost,
                                       kpar=if(kern_type=='laplacedot')
                                       {
                                         list(sigma=gamma_lap)
                                       }else if(kern_type=='rbfdot'){

                                         list(sigma=gamma_rbf)

                                       }else if(kern_type=='polydot'){
                                         list(degree=2,scale=poly_scale,offset=offset)
                                       }else{
                                         list(degree=1,scale=poly_scale,offset=offset)
                                       },
                                       epsilon=epsilon)})
  }
  #Calculando o predict para cada modelo
  predict <- lapply(early_model,function(x)predict(x,newdata=validation))

  #Calculating the weights (Equation 9)
  rmse <- unlist(lapply(predict,function(x){loss_function(predicted=x,observed=validation[,class_name])}))
  rmse<-rmse/stats::sd(rmse)
  # std_rmse<-rmse/(range(validation[,class_name])[2]-range(validation[,class_name])[1])
  inv_rmse<-(exp(-rmse*beta))

  prob_weights<-inv_rmse/sum(inv_rmse)
  prob_weights<-ifelse(prob_weights<0,0,prob_weights)#To not heve negative values of probabilities


  #----Defining the variables----
  models<-rep(list(0),B)#Creating the list of models
  boots_sample<-list(rep(B)) #Argument that will be passed in the map function
  out_of_bag<-list(rep(B)) #OOB samples object
  boots_index_row<- rep(list(nrow(train)),B)

  #====================================================

  #======Selecting the Bootstraping samples============
  #Defining which rows will be sampled
  boots_index_row<- lapply(boots_index_row,function(x){sample(1:x,x,replace=TRUE)})#Generating the boots_sample index


  #Defining out_of the bags_sample
  #Defining the Boots samples
  boots_sample<- lapply(boots_index_row,function(x){train[x,]}) #Without feature subsction
  out_of_bag<-lapply(boots_index_row,function(x){train[-unique(x),]})

  #=====================================================

  #=================Generating the models===============
  #Calculating the models

  #Here is defined which kernel will be used to heach model
  random_kernel<-sample(c('rbfdot','polydot','laplacedot','vanilladot'),
                        B,replace = TRUE,prob = prob_weights)

  if(automatic_tuning){
    models <- mapply(boots_sample,random_kernel,FUN = function(boot_sample,rand_kern){kernlab::ksvm(formula, data=boot_sample,type="eps-svr",
                                                  kernel=if(rand_kern=="vanilladot"){
                                                    "polydot"
                                                  }else{
                                                    rand_kern
                                                  },
                                                  C=cost,
                                                  kpar=if(rand_kern=='laplacedot' ||rand_kern=='rbfdot')
                                                  {
                                                    "automatic"
                                                  }else if(rand_kern=='polydot'){
                                                    list(degree=2,scale=poly_scale,offset=offset)
                                                  }else{
                                                    list(degree=1,scale=poly_scale,offset=offset)
                                                  },
                                                  epsilon=epsilon)})

  }else{
    models <- mapply(boots_sample,random_kernel, FUN  = function(boot_sample, rand_kern){kernlab::ksvm(formula, data=boot_sample,type="eps-svr",
                                                  kernel=if(rand_kern=="vanilladot"){
                                                    "polydot"
                                                  }else{
                                                    rand_kern
                                                  },
                                                  C=cost,
                                                  kpar=if(rand_kern=='laplacedot')
                                                  {
                                                    list(sigma=gamma_lap)
                                                  }else if(rand_kern=='rbfdot'){
                                                    list(sigma=gamma_rbf)
                                                  }else if(rand_kern=='polydot'){
                                                    list(degree=2,scale=poly_scale,offset=offset)
                                                  }else{
                                                    list(degree=1,scale=poly_scale,offset=offset)
                                                  },epsilon=epsilon)})

  }


  #Prediction of each mode ( for validation purpose)
  predict <- lapply(models,function(mod){predict(mod,newdata=validation)})
  predict_train <- lapply(models,function(mod){mod@fitted})

  #Prediction of OOB samples
  predict_oobg <- mapply(models,out_of_bag,FUN = function(mod,oob){predict(mod,newdata=oob)})

  #Calculating weights from equation 10
  kernel_weight<- mapply(predict_oobg,out_of_bag, FUN = function(pred_oob, oob){
    loss_function(predicted = pred_oob,observed = oob[,class_name])})

  boots_error<-mapply(predict_oobg,out_of_bag, FUN  = function(pred_oob,oob){loss_function(predicted = pred_oob,observed = oob[,class_name])})

  # Penalising by beta
  kernel_weight<- sapply(kernel_weight/stats::sd(kernel_weight),function(kern_weight){exp(-kern_weight*beta)})

  # Normalizing it
  kernel_weight_norm <-kernel_weight/sum((kernel_weight))

  #Predictions finals
  predict_df<- matrix(unlist(predict_train),ncol=nrow(train),byrow = TRUE) # Generating a matrix with where the the rows are each bootstrap sample
                                                                    #and the columns are each observation from test set

  # Accessing training error
  pred_df_train<-apply(mapply(models,kernel_weight_norm,FUN = function(mod, k_w_n){(predict(mod,train)*k_w_n)}),1,sum)#Multiplying the weights


  model_result <- list(y_train_hat=pred_df_train,lambda_values=list(Lin_Kern=prob_weights[4],
                                                                Pol_Kern=prob_weights[2],
                                                                RBF_Kern=prob_weights[1],
                                                                LAP_Kern=prob_weights[3]),
                       model_params=list(class_name=class_name,
                                         B=B,
                                         cost=cost,
                                         gamma=gamma,
                                         degree=degree),bootstrap_models=models,bootstrap_samples=boots_sample,
                       kernel_weight_norm=kernel_weight_norm,
                       predict_oob=predict_oobg)
  attr(model_result,"class")<-"rm_reg"
  #=============================
  return(model_result)
}




#' S4 class for RM classification
#'
#' @slot train a \code{data.frame} corresponding to the training data used into the model
#' @slot class_name a string with target variable used in the model
#' @slot kernel_weight a numeric vector corresponding to the weights for each bootstrap model contribution
#' @slot lambda_values a named list with value of the vector of \eqn{\boldsymbol{\lambda}} sampling probabilities associated with each each kernel function
#' @slot model_params a list with all used model specifications
#' @slot bootstrap_models a list with all \code{ksvm} objects for each bootstrap sample
#' @slot bootstrap_samples a list with all bootstrap samples used to train each base model of the ensemble
#' @slot prob a boolean indicating if a probabilitistic approch was used in the classification Random Machines
#' @details For more details see Ara, Anderson, et al. "Random machines: A bagged-weighted support vector model with free kernel choice." Journal of Data Science 19.3 (2021): 409-428.
#' @importFrom methods new
rm_class <- setClass("rm_class",
  slots = list(train = "data.frame",
               class_name = "character",
               kernel_weight = "numeric",
               lambda_values = "list",
               model_params = "list",
               bootstrap_models = "list",
               bootstrap_samples = "list")
)

#' @method predict rm_class
#' @export
predict.rm_class <- function(object, newdata,...) {
  # UseMethod(predict,rm_model)
  if (object$prob_model) {
    predict_new <- lapply(object$bootstrap_models, function(x) {predict(x, newdata = newdata, type = "probabilities")[, 2]})
    predict_df <- matrix(unlist(predict_new), ncol = nrow(newdata), byrow = TRUE)
    predict_df_new <- lapply(seq(1:nrow(newdata)), function(x) {predict_df[,x]})
    pred_df_fct <- lapply(predict_df_new, function(x) {stats::weighted.mean(x, object$kernel_weight)})
    return(unlist(pred_df_fct))

  } else {
    models <- object$bootstrap_models
    train <- object$train
    class_name <- object$class_name
    kernel_weight <- object$kernel_weight
    predict_new <- lapply(object$bootstrap_models, function(x){predict(x,newdata = newdata)})
    predict_df <- matrix(unlist(predict_new),ncol = nrow(newdata), byrow = TRUE)
    predict_df_new <- lapply(seq(1:nrow(newdata)), function(x){predict_df[,x]})
    pred_df_fct <- lapply(predict_df_new, function(x) {ifelse(x==unlist(levels(object$train[[object$class_name]]))[1], 1, -1)})
    pred_df_fct_final <- as.factor(unlist(lapply(pred_df_fct, function(x) {ifelse(sign(sum(x/((1+1e-10)-object$kernel_weight)^2))==1,levels(object$train[[object$class_name]])[1],levels(object$train[[object$class_name]])[2]) })))
    return(pred_df_fct_final)
  }
}

#' S4 class for RM regression
#' @slot y_train_hat a numeric  corresponding to the predictions \eqn{\hat{y}_{i}} for the training set
#' @slot lambda_values a named list with value of the vector of \eqn{\boldsymbol{\lambda}} sampling probabilities associated with each each kernel function
#' @slot model_params a list with all used model specifications
#' @slot bootstrap_models a list with all \code{ksvm} objects for each bootstrap sample
#' @slot bootstrap_samples a list with all bootstrap samples used to train each base model of the ensemble
#' @slot kernel_weight_norm a numeric vector corresponding to the normalised weights for each bootstrap model contribution
#' @importFrom methods new
#' @details For more details see Ara, Anderson, et al. "Regression random machines: An ensemble support vector regression model with free kernel choice." Expert Systems with Applications 202 (2022): 117107.
rm_reg <- setClass("rm_reg",
                     slots = list(y_train_hat = "numeric",
                                  lambda_values = "list",
                                  model_params = "list",
                                  bootstrap_models = "list",
                                  bootstrap_samples = "list")
)


#' @method predict rm_reg
#' @export
predict.rm_reg <- function(object, newdata,...) {
  # UseMethod(predict,rm_reg)
  # Accessing training error
  pred_df_test <- apply(mapply(object$bootstrap_models, object$kernel_weight_norm, FUN = function(mod, k_w_n) {(predict(mod, newdata) * k_w_n)}), 1, sum) # Multiplying the weights
  return(pred_df_test)
}

#' @export
#'
brier_score <- function(prob, observed, levels){
  y <- ifelse(observed==levels[1],1,0)
  b_score <- mean((y-prob)^2)
  return(b_score)
}

#' @export
RMSE<-function(predicted,observed){
  min<-min(observed)
  max<-max(observed)
  sqrt(mean(unlist((predicted-observed)^2)))
}
