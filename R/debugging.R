# # Importing some data
# library(kernlab)
# library(mlbench)
# library(magrittr)
#
# set.seed(42)
# # Taking a subset of it
# subset_spam_train <- mlbench.circle(n = 50,d = 2) %>% as.data.frame()
# subset_spam_validation <- mlbench.circle(n = 50,d = 2) %>% as.data.frame()
#
# formula <- classes ~ .
# train <- subset_spam_train
# validation <- subset_spam_validation
# boots_size = 25
# cost = 10
# seed.bootstrap = NULL
# automatic_tuning = TRUE
# gamma_rbf = 1
# gamma_lap = 1
# degree = 2
# poly_scale = 1
# offset = 0
# gamma_cau = 1
# d_t = 2
# kernels = c("rbfdot", "polydot", "laplacedot", "vanilladot")
# prob_model = F
#
#  # mod <- random_machines(formula = formula,train = train,validation = validation,kernels = kernels)
# # #
# sim_data <- rmachines::sim_class(n = 50)
# rm_mod <- rmachines::random_machines(y~., train = sim_data, validation =  sim_data,prob_model = F)
# rm_pred <- predict(rm_mod, newdata = sim_data)
# table(rm_pred,sim_data$y)
