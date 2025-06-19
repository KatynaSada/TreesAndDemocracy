# Best Parameters for Accuracy (as R data frames)
params_acc_rfmut <- data.frame(
  num_trees = 100,
  mtry = 4,
  min_node_size = 10,
  max_depth = 100,
  sample_fraction = 1
)

params_acc_rfmut_seats <- data.frame(
  num_trees = 500,
  mtry = 4,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 1,
  num_seats = 19,
  gamma = 1.6
)

params_acc_rfexp <- data.frame(
  num_trees = 100,
  mtry = 40,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 0.8
)

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
params_err_rfmut <- data.frame(
  num_trees = 200,
  mtry = 2,
  min_node_size = 3,
  max_depth = 100,
  sample_fraction = 0.8
)

params_err_rfmut_seats <- data.frame(
  num_trees = 200,
  mtry = 69,
  min_node_size = 10,
  max_depth = 100,
  sample_fraction = 0.7,
  num_seats = 19,
  gamma = 1.6
)

params_err_rfexp <- data.frame(
  num_trees = 100,
  mtry = 2,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 0.7
)

params_err_rfexp_seats <- data.frame(
  num_trees = 100,
  mtry = 8,
  min_node_size = 5,
  max_depth = 100,
  sample_fraction = 1,
  num_seats = 22,
  gamma = 1.6
)

# Best Parameters for XGBoost Accuracy (as R data frames)
params_acc_xgboostmut <- data.frame(
  max_depth = 5,
  