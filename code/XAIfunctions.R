library(randomForest)
library(xgboost)

#############################################################
# SEATS voting function
#############################################################

SEATS_assignment <- function(proportion, seats) {
  # Assign seats according to a vector of probabilities.
  
  # normalize the proportion vector to ensure it sums to one
  proportion <- proportion / sum(proportion)
  # calculate the initial seats assigned by rounding
  seats_assigned <- round(seats * proportion)
  
  # adjust the seats to ensure the total matches the specified number of seats
  # "Discrepancy" refers to the difference between the current proportion of seats assigned to a class and the desired proportion.
  while (sum(seats_assigned) != seats) {
    if (sum(seats_assigned) > seats) {
      # Largest Discrepancy: if there are too many seats, reduce one vote from the class with the largest discrepancy
      remove_from <- which.max(seats_assigned / seats - proportion)
      seats_assigned[remove_from] <- seats_assigned[remove_from] - 1
    } else {
      # Smallest Discrepancy:if there are too few seats, add one vote to the class with the smallest discrepancy
      add_to <- which.min(seats_assigned / seats - proportion)
      seats_assigned[add_to] <- seats_assigned[add_to] + 1
    }
  }
  
  # return the adjusted seats
  return(seats_assigned)
}

#############################################################
# Random Forest
#############################################################
require(randomForest)

train_rf_hard <- function(genomics_matrix, response_matrix, num_trees = 500, mtry = sqrt(ncol(genomics_matrix)), min_node_size = 1, max_depth = NULL, sample_fraction = 1) {
  # Train a Random Forest model using hard class assignments.
  # Args:
  #   genomics_matrix: A matrix of features (samples x features).
  #   response_matrix: A matrix of responses (samples x class_labels)
  #   num_trees: Number of trees in the Random Forest.
  #   mtry: Number of variables randomly sampled at each split
  #   min_node_size: Minimum size of terminal nodes
  #   max_depth: Maximum depth of trees (NULL means unlimited)
  #   sample_fraction: Fraction of cases to sample for each tree
  
  # Convert the response matrix to a factor of class labels (hard assignment)
  class_labels <- as.factor(colnames(response_matrix)[max.col(response_matrix)])
  
  # Train the Random Forest model
  model <- randomForest(x = genomics_matrix, 
             y = class_labels, 
             ntree = num_trees,
             mtry = mtry,
             nodesize = min_node_size,
             maxnodes = if(is.null(max_depth)) NULL else 2^max_depth,
             sampsize = ceiling(nrow(genomics_matrix) * sample_fraction))
  return(model)
}

train_rf_seats <- function(genomics_matrix, pseudoprobabilities, num_trees = 500, mtry = sqrt(ncol(genomics_matrix)), min_node_size = 1, max_depth = NULL, sample_fraction = 1, num_seats = 10) {
  # Expand the expression matrix
  expanded_expression <- genomics_matrix[rep(1:nrow(genomics_matrix), each = num_seats),]
  
  # Assign seats to each class
  seats_matrix <- apply(pseudoprobabilities, 1, SEATS_assignment, num_seats)
  
  # Create class labels vector
  class_labels <- as.factor(rep(rep(colnames(pseudoprobabilities), nrow(pseudoprobabilities)), as.numeric(seats_matrix)))
  
  # Train the Random Forest model
  model <- randomForest(x = expanded_expression, 
             y = class_labels, 
             ntree = num_trees,
             mtry = mtry,
             nodesize = min_node_size,
             maxnodes = if(is.null(max_depth)) NULL else 2^max_depth,
             sampsize = ceiling(nrow(expanded_expression) * sample_fraction))
  
  return(model)
}

predict_rf <- function(model, genomics_matrix) {
  # Predict class labels using a trained Random Forest model.
  # Args:
  #   model: A trained Random Forest model.
  #   genomics_matrix: A matrix of features (samples x features) for prediction.
  # Returns:
  #   A vector of predicted class labels as character strings.
  
  predictions <- as.character(predict(model, genomics_matrix))
  return(predictions)
}

# XGBoost
#############################################################
require(xgboost)

train_xgb_hard <- function(genomics_matrix, response_matrix, max_depth = 3, eta = 0.3, gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 8, early_stopping_rounds = 20) {
  # Train an XGBoost model using hard class assignments.
  # Args:
  #   genomics_matrix: A matrix of features (samples x features).
  #   response_matrix: A matrix of responses (samples x class_labels).
  #   max_depth: Maximum depth of trees
  #   eta: Learning rate
  #   gamma: Minimum loss reduction for partition
  #   min_child_weight: Minimum sum of instance weight in a child
  #   subsample: Subsample ratio of training instances
  #   colsample_bytree: Subsample ratio of columns for each tree
  #   nrounds: Number of boosting rounds
  #   nthread: Number of parallel threads
  #   early_stopping_rounds: Stop if performance doesn't improve
  
  class_names <- colnames(response_matrix)
  class_labels <- factor(class_names[max.col(response_matrix)], levels = class_names)
  numeric_labels <- as.numeric(class_labels) - 1
  
  xgb_params <- list("objective" = "multi:softmax",
           "eval_metric" = "mlogloss",
           "num_class" = ncol(response_matrix),
           "booster" = "gbtree")
  
  model <- xgboost(params = xgb_params, 
           data = genomics_matrix, 
           label = numeric_labels,
           max_depth = max_depth,
           eta = eta,
           gamma = gamma,
           min_child_weight = min_child_weight,
           subsample = subsample,
           colsample_bytree = colsample_bytree,
           nrounds = nrounds,
           nthread = nthread,
           early_stopping_rounds = early_stopping_rounds,
           verbose = 1)
  
  model$class_names <- class_names
  return(model)
}

train_xgb_seats <- function(genomics_matrix, pseudoprobabilities, max_depth = 3, eta = 0.3, gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 8, early_stopping_rounds = NULL, num_seats = 10) {
  
  expanded_expression <- genomics_matrix[rep(1:nrow(pseudoprobabilities), each = num_seats),]
  seats_matrix <- apply(pseudoprobabilities, 1, SEATS_assignment, num_seats)
  class_names <- colnames(pseudoprobabilities)
  class_labels <- factor(rep(rep(1:ncol(pseudoprobabilities), nrow(pseudoprobabilities)), 
               as.numeric(seats_matrix)), 
            levels = 1:ncol(pseudoprobabilities))
  numeric_labels <- as.numeric(class_labels) - 1

  xgb_params <- list("objective" = "multi:softmax",
           "eval_metric" = "mlogloss",
           "num_class" = ncol(pseudoprobabilities),
           "booster" = "gbtree")
           
  model <- xgboost(params = xgb_params, 
           data = expanded_expression, 
           label = numeric_labels,
           max_depth = max_depth,
           eta = eta,
           gamma = gamma,
           min_child_weight = min_child_weight,
           subsample = subsample,
           colsample_bytree = colsample_bytree,
           nrounds = nrounds,
           nthread = nthread,
           early_stopping_rounds = early_stopping_rounds,
           verbose = 1)
  
  model$class_names <- class_names
  return(model)
}

predict_xgb <- function(model, genomics_matrix) {
  # For multi:softmax, predict() directly returns class indices (0-based)
  predictions <- predict(model, genomics_matrix)
  # Add 1 to convert from 0-based to 1-based indexing and map to class names
  return(model$class_names[predictions + 1])
}
