# Best Parameters for Accuracy (as R data frames)
# params_acc_rfmut <- data.frame(
#   num_trees = 100,
#   mtry = 4,
#   min_node_size = 10,
#   max_depth = 100,
#   sample_fraction = 1
# )

params_acc_rfmut_seats <- data.frame(
  num_trees = 500,
  mtry = 4,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 1,
  num_seats = 19,
  gamma = 1.6
)

# params_acc_rfexp <- data.frame(
#   num_trees = 100,
#   mtry = 40,
#   min_node_size = 5,
#   max_depth = 100,
#   sample_fraction = 0.8
# )

params_acc_rfexp_seats <- data.frame(
  num_trees = 500,
  mtry = 2,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 1,
  num_seats = 25,
  gamma = 10
)

# Best Parameters for Error (as R data frames)
# params_err_rfmut <- data.frame(
#   num_trees = 200,
#   mtry = 2,
#   min_node_size = 3,
#   max_depth = 100,
#   sample_fraction = 0.8
# )

# params_err_rfmut_seats <- data.frame(
#   num_trees = 200,
#   mtry = 69,
#   min_node_size = 10,
#   max_depth = 100,
#   sample_fraction = 0.7,
#   num_seats = 19,
#   gamma = 1.6
# )

# params_err_rfexp <- data.frame(
#   num_trees = 100,
#   mtry = 2,
#   min_node_size = 5,
#   max_depth = 100,
#   sample_fraction = 0.7
# )

# params_err_rfexp_seats <- data.frame(
#   num_trees = 100,
#   mtry = 8,
#   min_node_size = 5,
#   max_depth = 100,
#   sample_fraction = 1,
#   num_seats = 22,
#   gamma = 1.6
# )

# Best Parameters for XGBoost Accuracy (as R data frames)
# params_acc_xgboostmut <- data.frame(
#   max_depth = 5,
#   eta = 0.01,
#   gamma = 1,
#   min_child_weight = 1,
#   subsample = 1,
#   colsample_bytree = 0.8,
#   nrounds = 200
# )
params_acc_xgboostmut <- data.frame( # use estos sin querer pero funcionan bien
  max_depth = 5,
  eta = 0.01,
  gamma = 1,
  min_child_weight = 1,
  subsample = 1,
  colsample_bytree = 0.8,
  nrounds = 200
)

params_acc_xgboostmut_seats <- data.frame(
  max_depth = 3,
  eta = 0.01,
  gamma = 1,
  min_child_weight = 5,
  subsample = 0.7,
  colsample_bytree = 0.8,
  nrounds = 100,
  num_seats = 9,
  gamma_seats = 10
)

# params_acc_xgboostexp <- data.frame(
#   max_depth = 7,
#   eta = 0.01,
#   gamma = 1,
#   min_child_weight = 1,
#   subsample = 0.8,
#   colsample_bytree = 1,
#   nrounds = 200
# )

# params_acc_xgboostexp_seats <- data.frame(
#   max_depth = 7,
#   eta = 0.01,
#   gamma = 0,
#   min_child_weight = 1,
#   subsample = 0.8,
#   colsample_bytree = 1,
#   nrounds = 2,
#   num_seats = 6,
#   gamma_seats = 20
# )

# Best Parameters for XGBoost Error (as R data frames)
# params_err_xgboostmut <- data.frame(
#   max_depth = 7,
#   eta = 0.1,
#   gamma = 0,
#   min_child_weight = 5,
#   subsample = 0.7,
#   colsample_bytree = 0.8,
#   nrounds = 50
# )

# params_err_xgboostmut_seats <- data.frame(
#   max_depth = 5,
#   eta = 0.01,
#   gamma = 5,
#   min_child_weight = 1,
#   subsample = 0.7,
#   colsample_bytree = 0.8,
#   nrounds = 100,
#   num_seats = 16,
#   gamma_seats = 1.6
# )

# params_err_xgboostexp <- data.frame(
#   max_depth = 7,
#   eta = 0.1,
#   gamma = 0,
#   min_child_weight = 5,
#   subsample = 0.7,
#   colsample_bytree = 0.8,
#   nrounds = 50
# )

params_err_xgboostexp_seats <- data.frame(
  max_depth = 7,
  eta = 0.01,
  gamma = 0,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 1,
  nrounds = 2,
  num_seats = 6,
  gamma_seats = 1.6
)
