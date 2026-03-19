#' Basic MR Family Estimation
#'
#' @param beta_hat_exp Exposure beta estimates (J x 3 matrix)
#' @param beta_hat_out Outcome beta estimates (J x 3 matrix)
#' @param beta_sigma_exp Exposure covariance matrices (list of 3x3)
#' @param beta_sigma_out Outcome covariance matrices (list of 3x3)
#' @param alpha_init Initial value for alpha (default 0)
#' @param max_iter Maximum iterations (default 100)
#' @param tol Convergence tolerance (default 1e-6)
#' @return List with alpha estimate and gamma values
run_mr_family <- function(
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    alpha_init = 0.0,
    max_iter = 100,
    tol = 1e-6
) {
    J <- nrow(beta_hat_exp)
    inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
    inv_beta_sigma_out <- lapply(beta_sigma_out, solve)
    alpha_k <- alpha_init

    for (k in 1:max_iter) {
        gamma_k_plus_1 <- matrix(0, nrow = J, ncol = 3)
        for (j in 1:J) {
            inv_beta_sigma_exp_j <- inv_beta_sigma_exp[[j]]
            inv_beta_sigma_out_j <- inv_beta_sigma_out[[j]]
            term1_inv <- inv_beta_sigma_exp_j + (alpha_k^2) * inv_beta_sigma_out_j
            term1 <- solve(term1_inv)
            term2 <- (inv_beta_sigma_exp_j %*% beta_hat_exp[j, ]) +
                (alpha_k * (inv_beta_sigma_out_j %*% beta_hat_out[j, ]))
            gamma_k_plus_1[j, ] <- term1 %*% term2
        }

        numerator <- 0.0
        denominator <- 0.0
        for (j in 1:J) {
            inv_beta_sigma_out_j <- inv_beta_sigma_out[[j]]
            gamma_j <- gamma_k_plus_1[j, ]
            numerator <- numerator + as.numeric(t(gamma_j) %*% inv_beta_sigma_out_j %*% beta_hat_out[j, ])
            denominator <- denominator + as.numeric(t(gamma_j) %*% inv_beta_sigma_out_j %*% gamma_j)
        }

        if (denominator == 0) break
        alpha_k_plus_1 <- numerator / denominator

        if (abs(alpha_k_plus_1 - alpha_k) < tol) {
            alpha_k <- alpha_k_plus_1
            gamma_k <- gamma_k_plus_1
            break
        }
        alpha_k <- alpha_k_plus_1
        gamma_k <- gamma_k_plus_1
    }
    list(alpha = alpha_k, gamma = gamma_k)
}

#' Calculate Variance for MR Alpha Estimate
#'
#' @param alpha_hat Estimated alpha value
#' @param beta_hat_out Outcome beta estimates
#' @param beta_hat_exp Exposure beta estimates
#' @param beta_sigma_out Outcome covariance matrices
#' @param beta_sigma_exp Exposure covariance matrices
#' @return List with variance and standard error
calculate_alpha_variance_internal_mr <- function(
    alpha_hat,
    beta_hat_out,
    beta_hat_exp,
    beta_sigma_out,
    beta_sigma_exp
) {
    J <- nrow(beta_hat_out)
    l_1 <- 0
    for (j in 1:J) {
        beta_sigma_out_j <- beta_sigma_out[[j]]
        beta_sigma_exp_j <- beta_sigma_exp[[j]]
        beta_hat_exp_j <- beta_hat_exp[j, ]
        u_j <- beta_hat_out[j, ] - alpha_hat * beta_hat_exp_j
        Sigma_j_alpha <- beta_sigma_out_j + alpha_hat^2 * beta_sigma_exp_j
        Sigma_j_alpha_inv <- solve(Sigma_j_alpha)
        A_j <- Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv
        term1 <- beta_hat_exp_j %*% Sigma_j_alpha_inv %*% beta_hat_exp_j - u_j %*% A_j %*% u_j
        term2 <- 4 * alpha_hat * beta_hat_exp_j %*% A_j %*% u_j
        term3 <- 2 * alpha_hat^2 * u_j %*% (A_j %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv +
                                                 Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% A_j) %*% u_j
        l_1 <- term1 + term2 + term3 + l_1
    }
    l <- l_1
    se <- 1 / sqrt(l)
    list(variance = se^2, std_error = se)
}

#' Mendelian Randomization Family Analysis
#'
#' Basic MR family analysis returning alpha estimate, standard error, and p-value
#'
#' @param beta_hat_exp Exposure beta estimates (J x 3 matrix)
#' @param beta_hat_out Outcome beta estimates (J x 3 matrix)
#' @param beta_sigma_exp Exposure covariance matrices (list of 3x3 matrices)
#' @param beta_sigma_out Outcome covariance matrices (list of 3x3 matrices)
#' @param alpha_init Initial alpha value (default 0)
#' @param max_iter Maximum iterations (default 100)
#' @param tol Convergence tolerance (default 1e-6)
#' @return Data frame with alpha_point, alpha_se, and p_value columns
#' @export
mr_family <- function(
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    alpha_init = 0.0,
    max_iter = 100,
    tol = 1e-6
) {
    a <- run_mr_family(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        beta_sigma_exp = beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        alpha_init = alpha_init,
        max_iter = max_iter,
        tol = tol
    )
    b <- calculate_alpha_variance_internal_mr(
        alpha_hat = a$alpha,
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        beta_sigma_exp = beta_sigma_exp,
        beta_sigma_out = beta_sigma_out
    )
    p_value <- 2 * pnorm(-abs(a$alpha / b$std_error))
    data.frame(alpha_point = a$alpha, alpha_se = b$std_error, p_value = p_value)
}
