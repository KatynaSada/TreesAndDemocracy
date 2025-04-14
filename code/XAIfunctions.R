library(randomForest)
library(xgboost)

votes_assignment <- function(proportion, votes=10) {
  # Assign votes according to a vector of probabilities.
  
  # normalize the proportion vector to ensure it sums to one
  proportion <- proportion / sum(proportion)
  # calculate the initial votes assigned by rounding
  votes_assigned <- round(votes * proportion)
  
  # adjust the votes to ensure the total matches the specified number of votes
  while (sum(votes_assigned) != votes) {
    if (sum(votes_assigned) > votes) {
      # if there are too many votes, reduce one vote from the class with the largest discrepancy
      remove_from <- which.max(votes_assigned / votes - proportion)
      votes_assigned[remove_from] <- votes_assigned[remove_from] - 1
    } else {
      # if there are too few votes, add one vote to the class with the smallest discrepancy
      add_to <- which.min(votes_assigned / votes - proportion)
      votes_assigned[add_to] <- votes_assigned[add_to] + 1
    }
  }
  
  # return the adjusted votes
  return(votes_assigned)
}

#############################################################
# Random Forest
#############################################################
require(randomForest)

train_rf_hard <- function(expression_matrix, response_matrix, num_trees = 500) {
  # Train a Random Forest model using hard class assignments.
  # Args:
  #   expression_matrix: A matrix of features (samples x features).
  #   response_matrix: A matrix of responses (samples x classes), where each row has a single 1 indicating the class.
  #   ntree: Number of trees in the Random Forest.
  # Returns:
  #   A trained Random Forest model.
  
  # Convert the response matrix to a factor of class labels (hard assignment).
  class_labels <- as.factor(max.col(response_matrix))
  
  # Train the Random Forest model.
  model <- randomForest(x = expression_matrix, y = class_labels, ntree = num_trees)
  return(model)
}

predict_rf_hard <- function(model, expression_matrix) {
  # Predict class labels using a trained Random Forest model.
  # Args:
  #   model: A trained Random Forest model.
  #   expression_matrix: A matrix of features (samples x features) for prediction.
  # Returns:
  #   A vector of predicted class labels.
  
  predictions <- predict(model, expression_matrix)
  return(predictions)
}

train_rf_votes <- function(expression_matrix, response_matrix, num_trees = 500, num_votes = 10) {
  # Train a Random Forest model where each drug is repeated according to the number of assigned votes.
  
  # Expand the expression matrix by repeating each row num_votes times.
  expanded_expression <- expression_matrix[rep(1:nrow(response_matrix), each = num_votes),]
  
  # Assign votes to each class (drug) for every sample.
  votes_matrix <- apply(response_matrix, 1, votes_assignment, num_votes)
  
  # Create a vector of class labels, repeated according to the assigned votes.
  class_labels <- as.factor(rep(rep(1:ncol(response_matrix), nrow(response_matrix)), as.numeric(votes_matrix)))
  
  # Train the Random Forest model.
  model <- randomForest(x = expanded_expression, y = class_labels, ntree = num_trees)
  
  return(model)
}

predict_rf <- function(model, expression_matrix) {
  predictions <- predict(model, expression_matrix)
  return(predictions)
}

# XGBoost
#############################################################
require(xgboost)

train_xgb_hard <- function(expression_matrix, response_matrix, num_rounds = 100, params = list()) {
  # Train an XGBoost model using hard class assignments.
  # Args:
  #   expression_matrix: A matrix of features (samples x features).
  #   response_matrix: A matrix of responses (samples x classes), where each row has a single 1 indicating the class.
  #   num_rounds: Number of boosting rounds.
  #   params: A list of parameters for the XGBoost model.
  # Returns:
  #   A trained XGBoost model.
  
  # Convert the response matrix to a factor of class labels (hard assignment).
  class_labels <- max.col(response_matrix) - 1  # XGBoost uses 0-based indexing for labels. 
  
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = ncol(response_matrix),
                    #"tree_method" = "gpu_hist",
                     "booster" = "gbtree")
  
  model <- xgboost(params = xgb_params, data = expression_matrix, 
                                   label = as.numeric(class_labels),
                                   max_depth = 3, eta = .3, nthread = 8, 
                                   nrounds = num_rounds, subsample = 1,
                                   print_every_n = 1)

  
  # Train the XGBoost model.
  #model <- xgboost(data = dtrain, params = params, nrounds = num_rounds, verbose = 0)
  return(model)
}

train_xgb_votes <- function(expression_matrix, response_matrix, num_rounds = 100, num_votes = 10, max_depth=3) {
  # Train an XGBoost model where each drug is repeated according to the number of assigned votes.

  # Expand the expression matrix by repeating each row num_votes times.
  expanded_expression <- expression_matrix[rep(1:nrow(response_matrix), each = num_votes),]
  votes_matrix <- apply(response_matrix,1,votes_assignment,num_votes)
  classes <- factor(rep(rep(1:ncol(response_matrix), nrow(response_matrix)), as.numeric(votes_matrix)), levels = 1:ncol(response_matrix))
  
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = ncol(response_matrix),
                     #"tree_method" = "gpu_hist",
                     "booster" = "gbtree")
  model <- xgboost(params = xgb_params, data = expanded_expression, 
                    label = as.numeric(classes)-1,
                    max_depth = max_depth, eta = .3, nthread = 8, 
                    nrounds = num_rounds, subsample = 1, verbose =1,
                    print_every_n = 1)
  
  
  return(model)
}

predict_xgb <- function(model, expression_matrix) {
  Treatment <- matrix(predict(model, expression_matrix), ncol = model$params$num_class, byrow = T)
  Treatment <- max.col(Treatment)
  return(Treatment)
}
