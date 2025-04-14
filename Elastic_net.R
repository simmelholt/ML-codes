rm(list = ls())

library(stargazer)
library(elasticnet)
library(quantmod)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(xts)
library(caret)
library(ggplot2)
library(pracma)
library(glmnet)




#Data
{
# Get factor returns and SPX returns
getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = Sys.Date())
spx_close <- Cl(GSPC)
spx_returns <- diff(spx_close) / lag(spx_close, 1)
spx_returns <- spx_returns[-1]
spx_returns <- round(spx_returns, 4)
spx_returns <- data.frame(date = index(spx_returns), returns = coredata(spx_returns))


file_path <- "C:/Users/ulrik/Downloads/[usa]_[all_factors]_[daily]_[vw_cap].csv"
factor_returns <- read.csv(file_path)
factor_returns$date <- as.Date(factor_returns$date, format="%Y-%m-%d")
factor_returns_filtered <- factor_returns %>%
  select(name, date, ret) %>%
  filter(date >= as.Date("1990-01-01"))
factor_returns <- factor_returns_filtered %>%
  pivot_wider(names_from = name, values_from = ret)
merged_data <- merge(spx_returns, factor_returns, by = "date", all = TRUE)
merged_data <- merged_data[order(merged_data$date), ]
colnames(merged_data)[2] <- "spx_ret"
merged_data <- merged_data[merged_data$date >= as.Date("1990-01-03"), ]


merged_data <- na.omit(merged_data)

#Demarket 
demarket <- function(r, mkt, b = NULL) {
  # Ensure inputs are matrices for matrix operations
  r <- as.matrix(r)
  mkt <- as.numeric(mkt)
  
  # Compute market beta if not provided
  if (is.null(b)) {
    b <- colMeans(r * mkt) / mean(mkt^2)  # Equivalent to regression slope
  }
  
  # De-market the returns
  rme <- r - outer(mkt, b)
  
  return(list(rme = rme, b = b))
}

# Apply demarket function to data from the 3rd column onward
demarket_returns <- demarket(merged_data[, 3:ncol(merged_data)], merged_data[, 2])
demarket_returns <- demarket_returns$rme

#merged_data <- demarket(merged_data[, 3:ncol(merged_data)], merged_data[, 2])
#merged_data <- merged_data$rme

#merged_data[, -2] <- merged_data[, -2] * 252
n <- nrow(demarket_returns)

# Compute the split index (80% train, 20% test)
split_index <- floor(0.6 * n)  

# Split into train (first 80%) and test (last 20%)
merged_data_train <- demarket_returns[1:split_index, ]
merged_data_test <- demarket_returns[(split_index + 1):n, ]



#Training data
#x_train <- as.matrix(merged_data_train[, 1:ncol(merged_data_train)])
#x_train_mu <- round(colMeans(x_test),5) 
#x_train_mu <- matrix(x_test_mu, ncol = 1)
#x_train_sigma <- cov(x_test)

#Test data
#x_test <- as.matrix(merged_data_test[, 1:ncol(merged_data_test)])
#x_test_mu <- round(colMeans(x_test),5) 
#x_test_mu <- matrix(x_test_mu, ncol = 1)
#x_test_sigma <- cov(x_test)

}

install.packages("HDRFA")
library(HDRFA)
# Step 1: Compute sample covariance matrix (no centering)
cov_matrix <- cov(merged_data_train*252)

# Step 2: Extract eigenvalues in decreasing order
eigenvalues <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values

# Step 3: Compute eigenvalue ratios
ratios <- eigenvalues[1:(length(eigenvalues) - 1)] / eigenvalues[2:length(eigenvalues)]

# Step 4: Estimate number of factors as index of max ratio
rhat <- which.max(ratios)

# Improved eigenvalue ratio plot without highlight
plot(1:length(ratios), ratios, type = "b", pch = 19, col = "blue",
     xlab = "Number of Factors (r)", ylab = "Eigenvalue Ratio λ₍ᵣ₎ / λ₍ᵣ₊₁₎",
     main = "Eigenvalue Ratios (Ahn & Horenstein)",
     cex.axis = 0.9, cex.lab = 1, cex.main = 1.1)

# Add grid for clarity
grid()




lambda_seq <- exp(seq(log(10^(-3)), log(10^(1)), length.out = 200))

{
  #Plot Enet path
  mu_plot <- colMeans(merged_data_train) * 252
  cov_plot <- cov(merged_data_train) * 252
  object <- enet(cov_plot, mu_plot, lambda = 0.1, max.steps = 153, normalize = FALSE, intercept = FALSE)
  par(mfrow=c(2,2))
  
  plot(object, xvar="step", use.color = TRUE)
  
  svd_result_plot <- svd(cov_plot)
  Q_plot <- svd_result_plot$v
  PC_rotation_plot <- merged_data_train %*% Q_plot
  mu_plot_PC <- colMeans(PC_rotation_plot) * 252
  cov_plot_PC <- cov(PC_rotation_plot) * 252
  object_PC <- enet(cov_plot_PC, mu_plot_PC, lambda = 0.1, max.steps = 153, normalize = FALSE, intercept = FALSE)
  plot(object_PC, xvar="step", use.color = TRUE)
}


library(purrr)
library(viridis)
library(akima)

demarket_df <- as.data.frame(demarket_returns)



#Enet function with raw data
{


log_min = -2
log_max = 1
  
merged_data_enet <- merged_data[, -c(1,2)]
  
lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = 400))



split_and_compute_stats <- function(data, train_frac = 0.8) {
  # Ensure only numeric columns are used (keeps all columns but filters non-numeric ones for calculations)
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Get number of rows
  n <- nrow(numeric_data)
  
  # Define the split index (70% for training, 30% for testing)
  split_index <- floor(train_frac * n)
  
  # Split the data chronologically
  train_data <- numeric_data[1:split_index, ]  # First 70% rows
  test_data <- numeric_data[(split_index + 1):n, ]  # Last 30% rows
  
  # Compute mean (mu) and covariance matrix (sigma) for training and test sets
  mu_train <- colMeans(train_data, na.rm = TRUE) * 252   # Multiply by 252
  sigma_train <- cov(train_data, use = "pairwise.complete.obs") * 252  # Multiply by 252
  mu_test <- colMeans(test_data, na.rm = TRUE) * 252   # Multiply by 252
  sigma_test <- cov(test_data, use = "pairwise.complete.obs") * 252  # Multiply by 252
  
  # Return only the mu and sigma values for both sets
  return(list(
    mu_train = mu_train,
    sigma_train = sigma_train,
    mu_test = mu_test,
    sigma_test = sigma_test
  ))
}
stats_result <- split_and_compute_stats(merged_data_enet)

mu_train <- stats_result$mu_train
sigma_train <- stats_result$sigma_train
mu_test <- stats_result$mu_test
sigma_test <- stats_result$sigma_test
elastic_net_analysis <- function(x_train_sigma, x_train_mu, x_test_sigma, x_test_mu) {
  # Define lambda sequence
  lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = 200))
  
  # Initialize results list
  L1L2 <- list(coef = list(), lambda = numeric(length(lambda_seq)))
  
  # Fit enet model for each lambda
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    
    enet_model <- enet(as.matrix(x_train_sigma), x_train_mu, 
                       lambda = lambda,
                       
                       normalize = FALSE, 
                       intercept = FALSE, 
                       trace = FALSE)
        
    
    # Store coefficients
    L1L2$coef[[as.character(lambda)]] <- enet_model$beta.pure
    L1L2$lambda[i] <- lambda
  }
  
  # Extract coefficients where the number of columns is 153
  all_coefs <- Filter(function(x) {
    dim_x <- dim(x)  # Get dimensions
    !is.null(dim_x) && dim_x[2] == 153
  }, L1L2$coef)
  
  # Ensure test data is in the correct format
  x_test_mu <- as.numeric(x_test_mu )  # Should be a numeric vector of length 151
  x_test_sigma <- as.matrix(x_test_sigma )  # Should be a 151 x 151 matrix
  
  # Compute denominator for R² (constant term)
  mu_norm_sq <- sum(x_test_mu^2)  # Equivalent to μ_2' * μ_2
  
  # Function to compute R² for each row of coefficients
  compute_R2 <- function(row_b) {
    row_b <- as.numeric(row_b)  # Convert to numeric vector (length 151)
    
    sigma_b <- x_test_sigma %*% row_b  # Compute Σ_2 * b̂ (151 x 1 result)
    diff_vec <- x_test_mu - sigma_b  # Compute (μ_2 - Σ_2 * b̂) (151 x 1)
    
    # Compute squared norm (numerator of R² equation)
    diff_norm_sq <- sum(diff_vec^2) 
    
    # Compute R²
    R2 <- 1 - (diff_norm_sq / mu_norm_sq)
    return(R2)
  }
  
  # Apply R² computation to each matrix in all_coefs
  all_coefs_R2 <- lapply(all_coefs, function(coef_matrix) {
    R2_values <- apply(coef_matrix, 1, compute_R2)  # Compute R² for each row
    coef_matrix <- cbind(coef_matrix, R2_values)  # Append R² as a new column
    return(coef_matrix)
  })
  
  # Find the highest R² for each lambda
  max_R2_values <- lapply(all_coefs_R2, function(coef_matrix) {
    max(coef_matrix[, ncol(coef_matrix)], na.rm = TRUE)  # Extract max R² (last column)
  })
  
  # Return all results
  return(list(
    coefficients = all_coefs,
    all_coefs_R2 = all_coefs_R2,
    max_R2_values = max_R2_values,
    lambda_values = L1L2$lambda
  ))
}
test <- elastic_net_analysis(sigma_train, mu_train, sigma_test, mu_test)


# Extract the all_coefs_R2 list from test
all_coefs_R2 <- test$all_coefs_R2  

# Define the lambda sequence (log-scaled X-axis)
lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = length(all_coefs_R2)))

# Find the minimum lambda value dynamically
min_lambda <- min(lambda_seq, na.rm = TRUE)

# Convert the list into a dataframe for plotting
heatmap_data <- do.call(rbind, lapply(seq_along(lambda_seq), function(i) {
  lambda <- lambda_seq[i]
  
  # Safe lookup in all_coefs_R2
  coef_matrix <- all_coefs_R2[[as.character(lambda)]]
  
  # If lambda key is missing, find the closest match
  if (is.null(coef_matrix)) {
    lambda_str <- names(all_coefs_R2)[which.min(abs(as.numeric(names(all_coefs_R2)) - lambda))]
    coef_matrix <- all_coefs_R2[[lambda_str]]
  }
  
  if (!is.null(coef_matrix) && !all(is.na(coef_matrix))) {  
    R2_values <- coef_matrix[, ncol(coef_matrix)]  # Extract R² values (last column)
    model_step <- 1:length(R2_values)  # Y-axis: 154 steps
    
    data.frame(Lambda = rep(lambda, length(R2_values)), Model_Step = model_step, R2 = R2_values)
  } else {
    NULL  # Skip if the matrix is NULL or contains only NAs
  }
}))

# Ensure there is valid data for plotting
if (nrow(heatmap_data) == 0) stop("ERROR: No valid data found for plotting.")

heatmap_filter <- heatmap_data[heatmap_data$Model_Step <= 153, ]
heatmap_filter$R2 <- ifelse(heatmap_filter$R2 < -0.1, -0.1, heatmap_filter$R2)

tau <- sum(diag(sigma_train))
T <- 153
gamma <- heatmap_filter[1]
heatmap_filter$kappa <- sqrt(tau / (gamma * T))
heatmap_filter$kappa <- as.numeric(heatmap_filter$kappa[[1]])




ggplot(heatmap_filter, aes(x = Lambda, y = Model_Step, fill = R2)) +
  geom_tile() +  
  scale_x_log10(limits = c(min_lambda, 1)) +  
  scale_y_continuous(limits = c(min(heatmap_filter$Model_Step, na.rm = TRUE), 153)) + 
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       limits = range(heatmap_filter$R2, na.rm = TRUE),
                       breaks = round(seq(min(heatmap_filter$R2, na.rm = TRUE), 
                                          max(heatmap_filter$R2, na.rm = TRUE), length.out = 5), 2),
                       labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
  labs(
    x = expression("L2 Shrinkage (" * lambda * ")"),
    y = "Number of non-zero coefficients",
    fill = expression(R^2)
  ) +
  theme_minimal(base_size = 14) + 
  theme(axis.title = element_text(face = "bold"))

ggplot(heatmap_filter, aes(x = kappa, y = Model_Step, fill = R2)) +
  geom_tile() +  
  scale_x_log10(limits = c(min(heatmap_filter$kappa, na.rm = TRUE), max(heatmap_filter$kappa, na.rm = TRUE))) +  
  scale_y_continuous(limits = c(min(heatmap_filter$Model_Step, na.rm = TRUE), 153)) + 
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       limits = range(heatmap_filter$R2, na.rm = TRUE),
                       breaks = round(seq(min(heatmap_filter$R2, na.rm = TRUE), 
                                          max(heatmap_filter$R2, na.rm = TRUE), length.out = 5), 2),
                       labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
  labs(
    x = expression("Root Expected SR"^2*" (Prior), " * kappa),
    y = "Number of non-zero coefficients",
    fill = expression(R^2)
  ) +
  theme_minimal(base_size = 14) + 
  theme(axis.title = element_text(face = "bold"))



}

y_vec <- as.numeric(sigma_train[, 1])
dim(x_mat)         # should be 153 x 1
length(y_vec)      # should be 153

x_mat <- as.matrix(sigma_train)
y_vec <- as.numeric(mu_train)

# Check dimensions
stopifnot(nrow(x_mat) == length(y_vec))

# Fit elastic net with cross-validation
library(elasticnet)

cv_fit <- cv.enet(
  x = x_mat,
  y = y_vec,
  K = 3,
  lambda = 0.1,                          # adjust if needed
  s = 1:153,
  mode = "step",
  plot.it = TRUE,
  se = TRUE
)





par(lwd = 2)

# Fit the elastic net model with specified options
fit_enet <- enet(
  x = sigma_train,
  y = mu_train,
  lambda = 1,
  max.steps = 153,
  intercept = FALSE,
  normalize = FALSE
)

# Plot the model with custom y-axis label
plot.enet(fit_enet, use.color = TRUE, xvar = "step")





#Enet with PC rotation 

{
  PC_merged  <- merged_data_enet 
  
  split_and_compute_stats <- function(data, train_frac = 0.8) {
    # Ensure only numeric columns are used (keeps all columns but filters non-numeric ones for calculations)
    numeric_data <- data[, sapply(data, is.numeric)]
    
    # Get number of rows
    n <- nrow(numeric_data)
    
    # Define the split index (train_frac% for training, rest for testing)
    split_index <- floor(train_frac * n)
    
    # Split the data chronologically
    train_data <- numeric_data[1:split_index, ]  # First train_frac% rows
    test_data <- numeric_data[(split_index + 1):n, ]   # Last (1 - train_frac)% rows
    
    train_data <- train_data 
    test_data <- test_data 
    # Compute mean (mu) and covariance matrix (sigma) for training and test sets
    mu_train <- colMeans(train_data, na.rm = TRUE) * 252  # Multiply by 252
    sigma_train <- cov(train_data, use = "pairwise.complete.obs") * 252   # Multiply by 252
    mu_test <- colMeans(test_data, na.rm = TRUE) * 252   # Multiply by 252
    sigma_test <- cov(test_data, use = "pairwise.complete.obs") * 252   # Multiply by 252
    
    # Return train and test data along with computed statistics
    return(list(
      train_data = train_data,   # Add train data
      test_data = test_data,     # Add test data
      mu_train = mu_train,
      sigma_train = sigma_train,
      mu_test = mu_test,
      sigma_test = sigma_test
    ))
  }
  
  stats_result_PC <- split_and_compute_stats(PC_merged)
  
  PC_merged_annual <- PC_merged  
  
  r_train_PC <- stats_result_PC$train_data
  r_test_PC <- stats_result_PC$test_data
  
  
  cov_matrix <- cov(r_train_PC)
  svd_result <- svd(cov_matrix)
  Q <- svd_result$v
  
  r_train_PC <- as.matrix(data.frame(lapply(r_train_PC, as.numeric)))
  r_test_PC <- as.matrix(data.frame(lapply(r_test_PC, as.numeric)))
  
  r_train_PC <- r_train_PC %*% Q
  r_test_PC <- r_test_PC %*% Q
  
  sigma_train_PC <- cov(r_train_PC)
  sigma_test_PC <- cov(r_test_PC)
  mu_train_PC <- colMeans(r_train_PC)
  mu_test_PC <- colMeans(r_test_PC)
  
  elastic_net_analysis <- function(x_train_sigma, x_train_mu, x_test_sigma, x_test_mu) {
    # Define lambda sequence
    lambda_seq <- exp(seq(log(10^(-8)), log(10^(0)), length.out = 200))
    
    # Initialize results list
    L1L2 <- list(coef = list(), lambda = numeric(length(lambda_seq)))
    
    # Fit enet model for each lambda
    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]
      
      enet_model <- enet(as.matrix(x_train_sigma), x_train_mu, 
                         lambda = lambda,
                         
                         normalize = FALSE, 
                         intercept = FALSE, 
                         trace = FALSE)
      
      
      # Store coefficients
      L1L2$coef[[as.character(lambda)]] <- enet_model$beta.pure
      L1L2$lambda[i] <- lambda
    }
    
    # Extract coefficients where the number of columns is 153
    all_coefs <- Filter(function(x) {
      dim_x <- dim(x)  # Get dimensions
      !is.null(dim_x) && dim_x[2] == 153
    }, L1L2$coef)
    
    # Ensure test data is in the correct format
    x_test_mu <- as.numeric(x_test_mu )  # Should be a numeric vector of length 151
    x_test_sigma <- as.matrix(x_test_sigma )  # Should be a 151 x 151 matrix
    
    # Compute denominator for R² (constant term)
    mu_norm_sq <- sum(x_test_mu^2)  # Equivalent to μ_2' * μ_2
    
    # Function to compute R² for each row of coefficients
    compute_R2 <- function(row_b) {
      row_b <- as.numeric(row_b)  # Convert to numeric vector (length 151)
      
      sigma_b <- x_test_sigma %*% row_b  # Compute Σ_2 * b̂ (151 x 1 result)
      diff_vec <- x_test_mu - sigma_b  # Compute (μ_2 - Σ_2 * b̂) (151 x 1)
      
      # Compute squared norm (numerator of R² equation)
      diff_norm_sq <- sum(diff_vec^2) 
      
      # Compute R²
      R2 <- 1 - (diff_norm_sq / mu_norm_sq)
      return(R2)
    }
    
    # Apply R² computation to each matrix in all_coefs
    all_coefs_R2 <- lapply(all_coefs, function(coef_matrix) {
      R2_values <- apply(coef_matrix, 1, compute_R2)  # Compute R² for each row
      coef_matrix <- cbind(coef_matrix, R2_values)  # Append R² as a new column
      return(coef_matrix)
    })
    
    # Find the highest R² for each lambda
    max_R2_values <- lapply(all_coefs_R2, function(coef_matrix) {
      max(coef_matrix[, ncol(coef_matrix)], na.rm = TRUE)  # Extract max R² (last column)
    })
    
    # Return all results
    return(list(
      coefficients = all_coefs,
      all_coefs_R2 = all_coefs_R2,
      max_R2_values = max_R2_values,
      lambda_values = L1L2$lambda
    ))
  }
  test_PC <- elastic_net_analysis(sigma_train_PC, mu_train_PC, sigma_test_PC, mu_test_PC)
  
  
  # Extract the all_coefs_R2 list from test
  all_coefs_R2_PC <- test_PC$all_coefs_R2  
  
  # Define the lambda sequence (log-scaled X-axis)
  lambda_seq_PC <- exp(seq(log(10^(-8)), log(10^(0)), length.out = length(all_coefs_R2_PC)))
  
  # Find the minimum lambda value dynamically
  min_lambda_PC <- min(lambda_seq_PC, na.rm = TRUE)
  
  # Convert the list into a dataframe for plotting
  heatmap_data_PC <- do.call(rbind, lapply(seq_along(lambda_seq_PC), function(i) {
    lambda_PC <- lambda_seq_PC[i]
    
    # Safe lookup in all_coefs_R2
    coef_matrix_PC <- all_coefs_R2_PC[[as.character(lambda_PC)]]
    
    # If lambda key is missing, find the closest match
    if (is.null(coef_matrix_PC)) {
      lambda_str <- names(all_coefs_R2_PC)[which.min(abs(as.numeric(names(all_coefs_R2_PC)) - lambda_PC))]
      coef_matrix_PC <- all_coefs_R2_PC[[lambda_str]]
    }
    
    if (!is.null(coef_matrix_PC) && !all(is.na(coef_matrix_PC))) {  
      R2_values_PC <- coef_matrix_PC[, ncol(coef_matrix_PC)]  # Extract R² values (last column)
      model_step_PC <- 1:length(R2_values_PC)  # Y-axis: 154 steps
      
      data.frame(Lambda_PC = rep(lambda_PC, length(R2_values_PC)), Model_Step_PC = model_step_PC, R2_PC = R2_values_PC)
    } else {
      NULL  
    }
  }))
}  
{
heatmap_data_PC <- heatmap_data_PC[heatmap_data_PC$Model_Step_PC <= 153, ]
heatmap_data_PC$R2_PC <- ifelse(heatmap_data_PC$R2_PC < -0.1, -0.1, heatmap_data_PC$R2_PC)

tau_PC <- sum(diag(sigma_train_PC))
T <- 153
gamma_PC <- heatmap_data_PC[1]
heatmap_data_PC$kappa_PC <- sqrt(tau_PC / (gamma_PC * T))
heatmap_data_PC$kappa_PC <- as.numeric(heatmap_data_PC$kappa_PC[[1]])


  # Ensure there is valid data for plotting
  if (nrow(heatmap_data_PC) == 0) stop("ERROR: No valid data found for plotting.")
  
  ggplot(heatmap_data_PC, aes(x = Lambda_PC, y = Model_Step_PC, fill = R2_PC)) +
    geom_tile() +  
    scale_x_log10(limits = c(min_lambda_PC, 10)) +  
    scale_y_continuous(limits = c(min(heatmap_data_PC$Model_Step_PC, na.rm = TRUE), 153)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_data_PC$R2_PC, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_data_PC$R2_PC, na.rm = TRUE), 
                                            max(heatmap_data_PC$R2_PC, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("L2 Shrinkage (" * lambda * ")"),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))
  
  ggplot(heatmap_data_PC, aes(x = kappa_PC, y = Model_Step_PC, fill = R2_PC)) +
    geom_tile() +  
    scale_x_log10(limits = c(min(heatmap_data_PC$kappa_PC, na.rm = TRUE), max(heatmap_data_PC$kappa_PC, na.rm = TRUE))) +  
    scale_y_continuous(limits = c(min(heatmap_data_PC$Model_Step_PC, na.rm = TRUE), 153)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_data_PC$R2_PC, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_data_PC$R2_PC, na.rm = TRUE), 
                                            max(heatmap_data_PC$R2_PC, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("Root Expected SR"^2*" (Prior), " * kappa),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))

  
  
  

}

{
ggplot(heatmap_data_PC, aes(x = kappa_PC, y = Model_Step_PC, fill = R2_PC)) +
  geom_raster(interpolate = TRUE) +  
  scale_x_log10(limits = c(min(heatmap_data_PC$kappa_PC, na.rm = TRUE), max(heatmap_data_PC$kappa_PC, na.rm = TRUE))) +  
  scale_y_log10(
    limits = c(max(1, min(heatmap_data_PC$Model_Step_PC, na.rm = TRUE)), 153), 
    breaks = c(1, 5, 10, 20, 50, 100, 153)
  ) + 
  scale_fill_viridis_c(option = "viridis", direction = 1, 
                       limits = range(heatmap_data_PC$R2_PC, na.rm = TRUE),
                       breaks = round(seq(min(heatmap_data_PC$R2_PC, na.rm = TRUE), 
                                          max(heatmap_data_PC$R2_PC, na.rm = TRUE), length.out = 5), 2),
                       labels = function(x) format(x, nsmall = 2)) +  
  labs(
    x = expression("Root Expected SR"^2*" (Prior), " * kappa),  
    aspect.ratio = 0.4  
  )}

#Sparcity chart
{
# Group by Model_Step_PC and find the max value of R2_PC
raw_vs_sparcity <- heatmap_data_PC %>%
  group_by(Model_Step_PC) %>%
  summarise(Max_R2_PC = max(R2_PC, na.rm = TRUE)) %>%
  ungroup()

raw_vs_sparcity_raw <- heatmap_data %>%
  group_by(Model_Step) %>%
  summarise(Max_R2 = max(R2, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Model_Step <= 153) 

# Rename Model_Step to Model_Step_PC to match raw_vs_sparcity
raw_vs_sparcity_raw <- raw_vs_sparcity_raw %>%
  rename(Model_Step_PC = Model_Step)

# Merge the two datasets
raw_vs_sparcity <- raw_vs_sparcity %>%
  full_join(raw_vs_sparcity_raw, by = "Model_Step_PC")

raw_vs_sparcity_long <- raw_vs_sparcity %>%
  pivot_longer(cols = c(Max_R2_PC, Max_R2), names_to = "Type", values_to = "R2_Value") %>%
  mutate(Type = recode(Type, 
                       Max_R2_PC = "Principal Components", 
                       Max_R2 = "Raw Factors"))

# Create the plot
ggplot(raw_vs_sparcity_long, aes(x = Model_Step_PC, y = R2_Value, color = Type)) +
  geom_line(size = 1) +  # Line plot for both series
  labs(
    title = "Sparsity",
    x = "Number of Coefficients",
    y = expression("OOS " * R^2)  # Properly format R squared
  ) +
  scale_color_manual(
    values = c("Principal Components" = "orange", "Raw Factors" = "lightblue")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top", # Move legend to top
    legend.title = element_blank(), # Remove "Legend" title
    text = element_text(size = 14) # Increase font size for readability
  )


}


#EXport data 
{
# Convert list of lists into a data frame
all_coefs_df <- do.call(rbind, lapply(test$all_coefs_R2, as.data.frame))

# Ensure R2 is numeric
all_coefs_df$R2 <- as.numeric(as.character(all_coefs_df$R2))

# Find the row with the highest R² value
max_r2_index <- which.max(all_coefs_df$R2)

# Extract the corresponding row
best_model <- all_coefs_df[max_r2_index, ]

# Display the result
print(best_model)

rownames(best_model) <- "Raw Factors"



# Convert list of lists into a data frame
all_coefs_df_PC <- do.call(rbind, lapply(test_PC$all_coefs_R2, as.data.frame))

# Ensure R2 is numeric
all_coefs_df_PC$R2 <- as.numeric(as.character(all_coefs_df_PC$R2))

# Find the row with the highest R² value
max_r2_index_PC <- which.max(all_coefs_df_PC$R2)

# Extract the corresponding row
best_model_PC <- all_coefs_df_PC[max_r2_index_PC, ]

# Display the result
print(best_model_PC)

rownames(best_model_PC) <- "PC Factors"

# Download The Coefs of the best Raw Model And best PC Model
# Save best_model as CSV
write.csv(best_model, file = "best_model.csv", row.names = TRUE)

# Save best_PC_model as CSV
write.csv(best_PC_model, file = "best_PC_model.csv", row.names = TRUE)


}

#FF data 
{
ff_data <- read.csv("C:/Users/ulrik/Downloads/12_Industry_Portfolios_Daily.csv", 
                    skip = 4, 
                    header = FALSE)

colnames(ff_data) <- c("Date", "NoDur", "Durbl", "Manuf", "Enrgy", "Chems", 
                       "BusEq", "Telcm", "Utils", "Shops", "Hlth", "Money", "Other")


ff_data$Date <- as.Date(as.character(ff_data$Date), format = "%Y%m%d")

ff_data <- ff_data[ff_data$Date >= as.Date("1990-01-01"), ]

ff_data <- ff_data[!is.na(ff_data$Date), ]

getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = Sys.Date())
spx_close <- Cl(GSPC)
spx_returns <- diff(spx_close) / lag(spx_close, 1)
spx_returns <- spx_returns[-1]
spx_returns <- round(spx_returns, 4)
spx_returns <- data.frame(date = index(spx_returns), returns = coredata(spx_returns))

colnames(spx_returns)[colnames(spx_returns) == "date"] <- "Date"

merged_data_FF <- merge(ff_data, spx_returns, by = "Date")
#merged_data_FF[, 14] <- merged_data[, 14] * 252
merged_data_enet_ff <- merged_data_FF[, -c(1, ncol(merged_data_FF))]

sapply(merged_data_enet_ff, is.numeric)
# Convert all columns to numeric (if possible)
merged_data_enet_ff[] <- lapply(merged_data_enet_ff, function(col) as.numeric(as.character(col)))

# Now test again
sapply(merged_data_enet_ff, is.numeric)

merged_data_enet_ff <- merged_data_enet_ff /100
}
#Enet function with FF data 
{
  
  
  log_min = -7
  log_max = 0
  
  
  
  lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = 1000))
  
  
  
  split_and_compute_stats <- function(data, train_frac = 0.9) {
    # Ensure only numeric columns are used (keeps all columns but filters non-numeric ones for calculations)
    numeric_data <- data[, sapply(data, is.numeric)]
    
    # Get number of rows
    n <- nrow(numeric_data)
    
    # Define the split index (70% for training, 30% for testing)
    split_index <- floor(train_frac * n)
    
    # Split the data chronologically
    train_data <- numeric_data[1:split_index, ]  # First 70% rows
    test_data <- numeric_data[(split_index + 1):n, ]  # Last 30% rows
    
    # Compute mean (mu) and covariance matrix (sigma) for training and test sets
    mu_train <- colMeans(train_data, na.rm = TRUE)   # Multiply by 252
    sigma_train <- cov(train_data, use = "pairwise.complete.obs")   # Multiply by 252
    mu_test <- colMeans(test_data, na.rm = TRUE)    # Multiply by 252
    sigma_test <- cov(test_data, use = "pairwise.complete.obs")   # Multiply by 252
    
    # Return only the mu and sigma values for both sets
    return(list(
      mu_train = mu_train,
      sigma_train = sigma_train,
      mu_test = mu_test,
      sigma_test = sigma_test
    ))
  }
  stats_result <- split_and_compute_stats(merged_data_enet_ff)
  
  mu_train <- stats_result$mu_train
  sigma_train <- stats_result$sigma_train
  mu_test <- stats_result$mu_test
  sigma_test <- stats_result$sigma_test
  elastic_net_analysis <- function(x_train_sigma, x_train_mu, x_test_sigma, x_test_mu) {
    # Define lambda sequence
    lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = 1000))
    
    # Initialize results list
    L1L2 <- list(coef = list(), lambda = numeric(length(lambda_seq)))
    
    # Fit enet model for each lambda
    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]
      
      enet_model <- enet(as.matrix(x_train_sigma), x_train_mu, 
                         lambda = lambda,
                        
                         normalize = FALSE, 
                         intercept = FALSE, 
                         trace = FALSE)
      
      
      # Store coefficients
      L1L2$coef[[as.character(lambda)]] <- enet_model$beta.pure
      L1L2$lambda[i] <- lambda
    }
    
    # Extract coefficients where the number of columns is 153
    all_coefs <- Filter(function(x) {
      dim_x <- dim(x)  # Get dimensions
      !is.null(dim_x) && dim_x[2] == 12
    }, L1L2$coef)
    
    # Ensure test data is in the correct format
    x_test_mu <- as.numeric(x_test_mu )  # Should be a numeric vector of length 151
    x_test_sigma <- as.matrix(x_test_sigma )  # Should be a 151 x 151 matrix
    
    # Compute denominator for R² (constant term)
    mu_norm_sq <- sum(x_test_mu^2)  # Equivalent to μ_2' * μ_2
    
    # Function to compute R² for each row of coefficients
    compute_R2 <- function(row_b) {
      row_b <- as.numeric(row_b)  # Convert to numeric vector (length 151)
      
      sigma_b <- x_test_sigma %*% row_b  # Compute Σ_2 * b̂ (151 x 1 result)
      diff_vec <- x_test_mu - sigma_b  # Compute (μ_2 - Σ_2 * b̂) (151 x 1)
      
      # Compute squared norm (numerator of R² equation)
      diff_norm_sq <- sum(diff_vec^2) 
      
      # Compute R²
      R2 <- 1 - (diff_norm_sq / mu_norm_sq)
      return(R2)
    }
    
    # Apply R² computation to each matrix in all_coefs
    all_coefs_R2 <- lapply(all_coefs, function(coef_matrix) {
      R2_values <- apply(coef_matrix, 1, compute_R2)  # Compute R² for each row
      coef_matrix <- cbind(coef_matrix, R2_values)  # Append R² as a new column
      return(coef_matrix)
    })
    
    # Find the highest R² for each lambda
    max_R2_values <- lapply(all_coefs_R2, function(coef_matrix) {
      max(coef_matrix[, ncol(coef_matrix)], na.rm = TRUE)  # Extract max R² (last column)
    })
    
    # Return all results
    return(list(
      coefficients = all_coefs,
      all_coefs_R2 = all_coefs_R2,
      max_R2_values = max_R2_values,
      lambda_values = L1L2$lambda
    ))
  }
  test <- elastic_net_analysis(sigma_train, mu_train, sigma_test, mu_test)
  
  
  # Extract the all_coefs_R2 list from test
  all_coefs_R2 <- test$all_coefs_R2  
  
  # Define the lambda sequence (log-scaled X-axis)
  lambda_seq <- exp(seq(log(10^(log_min)), log(10^(log_max)), length.out = length(all_coefs_R2)))
  
  # Find the minimum lambda value dynamically
  min_lambda <- min(lambda_seq, na.rm = TRUE)
  
  # Convert the list into a dataframe for plotting
  heatmap_data <- do.call(rbind, lapply(seq_along(lambda_seq), function(i) {
    lambda <- lambda_seq[i]
    
    # Safe lookup in all_coefs_R2
    coef_matrix <- all_coefs_R2[[as.character(lambda)]]
    
    # If lambda key is missing, find the closest match
    if (is.null(coef_matrix)) {
      lambda_str <- names(all_coefs_R2)[which.min(abs(as.numeric(names(all_coefs_R2)) - lambda))]
      coef_matrix <- all_coefs_R2[[lambda_str]]
    }
    
    if (!is.null(coef_matrix) && !all(is.na(coef_matrix))) {  
      R2_values <- coef_matrix[, ncol(coef_matrix)]  # Extract R² values (last column)
      model_step <- 1:length(R2_values)  # Y-axis: 154 steps
      
      data.frame(Lambda = rep(lambda, length(R2_values)), Model_Step = model_step, R2 = R2_values)
    } else {
      NULL  # Skip if the matrix is NULL or contains only NAs
    }
  }))
  
  # Ensure there is valid data for plotting
  if (nrow(heatmap_data) == 0) stop("ERROR: No valid data found for plotting.")
  
  heatmap_filter <- heatmap_data[heatmap_data$Model_Step <= 12, ]
  heatmap_filter$R2 <- ifelse(heatmap_filter$R2 < -0.1, -0.1, heatmap_filter$R2)
  
  tau <- sum(diag(sigma_train))
  T <- 12
  gamma <- heatmap_filter[1]
  heatmap_filter$kappa <- sqrt(tau / (gamma * T))
  heatmap_filter$kappa <- as.numeric(heatmap_filter$kappa[[1]])
  
  
  
  
  ggplot(heatmap_filter, aes(x = Lambda, y = Model_Step, fill = R2)) +
    geom_tile() +  
    scale_x_log10(limits = c(min_lambda, 1)) +  
    scale_y_continuous(limits = c(min(heatmap_filter$Model_Step, na.rm = TRUE), 12)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_filter$R2, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_filter$R2, na.rm = TRUE), 
                                            max(heatmap_filter$R2, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("L2 Shrinkage (" * lambda * ")"),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))
  
  ggplot(heatmap_filter, aes(x = kappa, y = Model_Step, fill = R2)) +
    geom_tile() +  
    scale_x_log10(limits = c(min(heatmap_filter$kappa, na.rm = TRUE), max(heatmap_filter$kappa, na.rm = TRUE))) +  
    scale_y_continuous(limits = c(min(heatmap_filter$Model_Step, na.rm = TRUE), 12)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_filter$R2, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_filter$R2, na.rm = TRUE), 
                                            max(heatmap_filter$R2, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("Root Expected SR"^2*" (Prior), " * kappa),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))
  
  
  
}

#FF PC rotation

{
  PC_merged  <- merged_data_enet_ff 
  
  split_and_compute_stats <- function(data, train_frac = 0.9) {
    # Ensure only numeric columns are used (keeps all columns but filters non-numeric ones for calculations)
    numeric_data <- data[, sapply(data, is.numeric)]
    
    # Get number of rows
    n <- nrow(numeric_data)
    
    # Define the split index (train_frac% for training, rest for testing)
    split_index <- floor(train_frac * n)
    
    # Split the data chronologically
    train_data <- numeric_data[1:split_index, ]  # First train_frac% rows
    test_data <- numeric_data[(split_index + 1):n, ]   # Last (1 - train_frac)% rows
    
    train_data <- train_data 
    test_data <- test_data 
    # Compute mean (mu) and covariance matrix (sigma) for training and test sets
    mu_train <- colMeans(train_data, na.rm = TRUE) * 252  # Multiply by 252
    sigma_train <- cov(train_data, use = "pairwise.complete.obs") * 252   # Multiply by 252
    mu_test <- colMeans(test_data, na.rm = TRUE) * 252   # Multiply by 252
    sigma_test <- cov(test_data, use = "pairwise.complete.obs") * 252   # Multiply by 252
    
    # Return train and test data along with computed statistics
    return(list(
      train_data = train_data,   # Add train data
      test_data = test_data,     # Add test data
      mu_train = mu_train,
      sigma_train = sigma_train,
      mu_test = mu_test,
      sigma_test = sigma_test
    ))
  }
  
  stats_result_PC <- split_and_compute_stats(PC_merged)
  
  PC_merged_annual <- PC_merged  
  
  r_train_PC <- stats_result_PC$train_data
  r_test_PC <- stats_result_PC$test_data
  
  
  cov_matrix <- cov(r_train_PC)
  svd_result <- svd(cov_matrix)
  Q <- svd_result$v
  
  r_train_PC <- as.matrix(data.frame(lapply(r_train_PC, as.numeric)))
  r_test_PC <- as.matrix(data.frame(lapply(r_test_PC, as.numeric)))
  
  r_train_PC <- r_train_PC %*% Q
  r_test_PC <- r_test_PC %*% Q
  
  sigma_train_PC <- cov(r_train_PC)
  sigma_test_PC <- cov(r_test_PC)
  mu_train_PC <- colMeans(r_train_PC)
  mu_test_PC <- colMeans(r_test_PC)
  
  elastic_net_analysis <- function(x_train_sigma, x_train_mu, x_test_sigma, x_test_mu) {
    # Define lambda sequence
    lambda_seq <- exp(seq(log(10^(-8)), log(10^(0)), length.out = 200))
    
    # Initialize results list
    L1L2 <- list(coef = list(), lambda = numeric(length(lambda_seq)))
    
    # Fit enet model for each lambda
    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]
      
      enet_model <- enet(as.matrix(x_train_sigma), x_train_mu, 
                         lambda = lambda,
                         
                         normalize = FALSE, 
                         intercept = FALSE, 
                         trace = FALSE)
      
      
      # Store coefficients
      L1L2$coef[[as.character(lambda)]] <- enet_model$beta.pure
      L1L2$lambda[i] <- lambda
    }
    
    # Extract coefficients where the number of columns is 153
    all_coefs <- Filter(function(x) {
      dim_x <- dim(x)  # Get dimensions
      !is.null(dim_x) && dim_x[2] == 12
    }, L1L2$coef)
    
    # Ensure test data is in the correct format
    x_test_mu <- as.numeric(x_test_mu )  # Should be a numeric vector of length 151
    x_test_sigma <- as.matrix(x_test_sigma )  # Should be a 151 x 151 matrix
    
    # Compute denominator for R² (constant term)
    mu_norm_sq <- sum(x_test_mu^2)  # Equivalent to μ_2' * μ_2
    
    # Function to compute R² for each row of coefficients
    compute_R2 <- function(row_b) {
      row_b <- as.numeric(row_b)  # Convert to numeric vector (length 151)
      
      sigma_b <- x_test_sigma %*% row_b  # Compute Σ_2 * b̂ (151 x 1 result)
      diff_vec <- x_test_mu - sigma_b  # Compute (μ_2 - Σ_2 * b̂) (151 x 1)
      
      # Compute squared norm (numerator of R² equation)
      diff_norm_sq <- sum(diff_vec^2) 
      
      # Compute R²
      R2 <- 1 - (diff_norm_sq / mu_norm_sq)
      return(R2)
    }
    
    # Apply R² computation to each matrix in all_coefs
    all_coefs_R2 <- lapply(all_coefs, function(coef_matrix) {
      R2_values <- apply(coef_matrix, 1, compute_R2)  # Compute R² for each row
      coef_matrix <- cbind(coef_matrix, R2_values)  # Append R² as a new column
      return(coef_matrix)
    })
    
    # Find the highest R² for each lambda
    max_R2_values <- lapply(all_coefs_R2, function(coef_matrix) {
      max(coef_matrix[, ncol(coef_matrix)], na.rm = TRUE)  # Extract max R² (last column)
    })
    
    # Return all results
    return(list(
      coefficients = all_coefs,
      all_coefs_R2 = all_coefs_R2,
      max_R2_values = max_R2_values,
      lambda_values = L1L2$lambda
    ))
  }
  test_PC <- elastic_net_analysis(sigma_train_PC, mu_train_PC, sigma_test_PC, mu_test_PC)
  
  
  # Extract the all_coefs_R2 list from test
  all_coefs_R2_PC <- test_PC$all_coefs_R2  
  
  # Define the lambda sequence (log-scaled X-axis)
  lambda_seq_PC <- exp(seq(log(10^(-8)), log(10^(0)), length.out = length(all_coefs_R2_PC)))
  
  # Find the minimum lambda value dynamically
  min_lambda_PC <- min(lambda_seq_PC, na.rm = TRUE)
  
  # Convert the list into a dataframe for plotting
  heatmap_data_PC <- do.call(rbind, lapply(seq_along(lambda_seq_PC), function(i) {
    lambda_PC <- lambda_seq_PC[i]
    
    # Safe lookup in all_coefs_R2
    coef_matrix_PC <- all_coefs_R2_PC[[as.character(lambda_PC)]]
    
    # If lambda key is missing, find the closest match
    if (is.null(coef_matrix_PC)) {
      lambda_str <- names(all_coefs_R2_PC)[which.min(abs(as.numeric(names(all_coefs_R2_PC)) - lambda_PC))]
      coef_matrix_PC <- all_coefs_R2_PC[[lambda_str]]
    }
    
    if (!is.null(coef_matrix_PC) && !all(is.na(coef_matrix_PC))) {  
      R2_values_PC <- coef_matrix_PC[, ncol(coef_matrix_PC)]  # Extract R² values (last column)
      model_step_PC <- 1:length(R2_values_PC)  # Y-axis: 154 steps
      
      data.frame(Lambda_PC = rep(lambda_PC, length(R2_values_PC)), Model_Step_PC = model_step_PC, R2_PC = R2_values_PC)
    } else {
      NULL  
    }
  }))
}  
{
  heatmap_data_PC <- heatmap_data_PC[heatmap_data_PC$Model_Step_PC <= 12, ]
  heatmap_data_PC$R2_PC <- ifelse(heatmap_data_PC$R2_PC < -0.1, -0.1, heatmap_data_PC$R2_PC)
  
  tau_PC <- sum(diag(sigma_train_PC))
  T <- 12
  gamma_PC <- heatmap_data_PC[1]
  heatmap_data_PC$kappa_PC <- sqrt(tau_PC / (gamma_PC * T))
  heatmap_data_PC$kappa_PC <- as.numeric(heatmap_data_PC$kappa_PC[[1]])
  
  
  # Ensure there is valid data for plotting
  if (nrow(heatmap_data_PC) == 0) stop("ERROR: No valid data found for plotting.")
  
  ggplot(heatmap_data_PC, aes(x = Lambda_PC, y = Model_Step_PC, fill = R2_PC)) +
    geom_tile() +  
    scale_x_log10(limits = c(min_lambda_PC, 10)) +  
    scale_y_continuous(limits = c(min(heatmap_data_PC$Model_Step_PC, na.rm = TRUE), 12)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_data_PC$R2_PC, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_data_PC$R2_PC, na.rm = TRUE), 
                                            max(heatmap_data_PC$R2_PC, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("L2 Shrinkage (" * lambda * ")"),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))
  
  ggplot(heatmap_data_PC, aes(x = kappa_PC, y = Model_Step_PC, fill = R2_PC)) +
    geom_tile() +  
    scale_x_log10(limits = c(min(heatmap_data_PC$kappa_PC, na.rm = TRUE), max(heatmap_data_PC$kappa_PC, na.rm = TRUE))) +  
    scale_y_continuous(limits = c(min(heatmap_data_PC$Model_Step_PC, na.rm = TRUE), 12)) + 
    scale_fill_viridis_c(option = "viridis", direction = 1, 
                         limits = range(heatmap_data_PC$R2_PC, na.rm = TRUE),
                         breaks = round(seq(min(heatmap_data_PC$R2_PC, na.rm = TRUE), 
                                            max(heatmap_data_PC$R2_PC, na.rm = TRUE), length.out = 5), 2),
                         labels = function(x) format(x, nsmall = 2)) +  # Keep labels at 2 decimals
    labs(
      x = expression("Root Expected SR"^2*" (Prior), " * kappa),
      y = "Number of non-zero coefficients",
      fill = expression(R^2)
    ) +
    theme_minimal(base_size = 14) + 
    theme(axis.title = element_text(face = "bold"))
  
  
  
  
  
}





# Group by Model_Step_PC and find the max value of R2_PC
raw_vs_sparcity <- heatmap_data_PC %>%
  group_by(Model_Step_PC) %>%
  summarise(Max_R2_PC = max(R2_PC, na.rm = TRUE)) %>%
  ungroup()

raw_vs_sparcity_raw <- heatmap_data %>%
  group_by(Model_Step) %>%
  summarise(Max_R2 = max(R2, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Model_Step <= 153) 

# Rename Model_Step to Model_Step_PC to match raw_vs_sparcity
raw_vs_sparcity_raw <- raw_vs_sparcity_raw %>%
  rename(Model_Step_PC = Model_Step)

# Merge the two datasets
raw_vs_sparcity <- raw_vs_sparcity %>%
  full_join(raw_vs_sparcity_raw, by = "Model_Step_PC")

raw_vs_sparcity_long <- raw_vs_sparcity %>%
  pivot_longer(cols = c(Max_R2_PC, Max_R2), names_to = "Type", values_to = "R2_Value") %>%
  mutate(Type = recode(Type, 
                       Max_R2_PC = "Principal Components", 
                       Max_R2 = "Raw Factors"))

# Create the plot
ggplot(raw_vs_sparcity_long, aes(x = Model_Step_PC, y = R2_Value, color = Type)) +
  geom_line(size = 1) +  # Line plot for both series
  labs(
    title = "Sparsity",
    x = "Number of Coefficients",
    y = expression("OOS " * R^2)  # Properly format R squared
  ) +
  scale_color_manual(
    values = c("Principal Components" = "orange", "Raw Factors" = "lightblue")
  ) +
  theme_minimal() +
  theme(
    legend.position = "top", # Move legend to top
    legend.title = element_blank(), # Remove "Legend" title
    text = element_text(size = 14) # Increase font size for readability
  )




