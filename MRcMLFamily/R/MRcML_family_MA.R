#' Helper function to invert beta_sigma with masking
inv_M_j <- function(M_j, beta_sigma) {
    beta_sigma_clean <- M_j %*% beta_sigma %*% M_j
    if (all(beta_sigma_clean == 0)) {
        return(matrix(0, ncol = 3, nrow = 3))
    }
    non_zero_rows <- rowSums(abs(beta_sigma_clean)) > 0
    non_zero_cols <- colSums(abs(beta_sigma_clean)) > 0
    indices_to_keep <- non_zero_rows & non_zero_cols
    beta_sigma_clean <- beta_sigma_clean[indices_to_keep, indices_to_keep]
    inv_results <- matrix(0, ncol = 3, nrow = 3)
    inv_results[indices_to_keep, indices_to_keep] <- solve(beta_sigma_clean)
    inv_results
}

#' Internal cML Family Iterative Algorithm
#'
#' @keywords internal
run_cml_family_internal <- function(
    beta_hat_exp,
    beta_hat_out,
    inv_beta_sigma_exp,
    beta_sigma_out,
    a,
    b,
    alpha_init,
    max_iter,
    tol
) {
    J <- nrow(beta_hat_exp)
    alpha_k <- alpha_init
    gamma_k <- beta_hat_exp

    for (k in 1:max_iter) {
        f_values <- numeric(J)
        t_values <- numeric(J)
        for (j in 1:J) {
            sigma_beta_o_sq <- beta_sigma_out[[j]][1, 1]
            f_values[j] <- (beta_hat_out[j, 1] - alpha_k * gamma_k[j, 1])^2 / sigma_beta_o_sq
            d_fm <- beta_hat_out[j, 2:3] - alpha_k * gamma_k[j, 2:3]
            sigma_beta_fm <- beta_sigma_out[[j]][2:3, 2:3]
            inv_sigma_beta_fm <- solve(sigma_beta_fm)
            t_values[j] <- as.numeric(t(d_fm) %*% inv_sigma_beta_fm %*% d_fm)
        }
        sf_indices <- order(f_values)
        st_indices <- order(t_values)
        S_f <- sf_indices[1:a]
        S_t <- st_indices[1:b]
        if (a == 0) S_f <- NULL
        if (b == 0) S_t <- NULL

        gamma_k_plus_1 <- matrix(0, nrow = J, ncol = 3)
        for (j in 1:J) {
            M_j <- diag(c(
                ifelse(j %in% S_f, 1, 0),
                ifelse(j %in% S_t, 1, 0),
                ifelse(j %in% S_t, 1, 0)
            ))
            inv_beta_sigma_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
            term1_inv <- inv_beta_sigma_exp[[j]] + (alpha_k^2) * inv_beta_sigma_out_j
            term1 <- solve(term1_inv)
            term2 <- (inv_beta_sigma_exp[[j]] %*% beta_hat_exp[j, ]) +
                (alpha_k * (inv_beta_sigma_out_j %*% beta_hat_out[j, ]))
            gamma_k_plus_1[j, ] <- term1 %*% term2
        }

        numerator <- 0.0
        denominator <- 0.0
        for (j in 1:J) {
            M_j <- diag(c(
                ifelse(j %in% S_f, 1, 0),
                ifelse(j %in% S_t, 1, 0),
                ifelse(j %in% S_t, 1, 0)
            ))
            gamma_j <- gamma_k_plus_1[j, ]
            inv_beta_simga_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
            numerator <- numerator + as.numeric(t(gamma_j) %*% inv_beta_simga_out_j %*% beta_hat_out[j, ])
            denominator <- denominator + as.numeric(t(gamma_j) %*% inv_beta_simga_out_j %*% gamma_j)
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

    list(alpha = alpha_k, gamma = gamma_k, S_f = S_f, S_t = S_t)
}

#' Run cML Family with Pre-computed Inverses
#'
#' @param beta_hat_exp Exposure beta estimates (J x 3)
#' @param beta_hat_out Outcome beta estimates (J x 3)
#' @param beta_sigma_exp Exposure covariance matrices (list)
#' @param beta_sigma_out Outcome covariance matrices (list)
#' @param a Number of top exposure SNPs to select
#' @param b Number of top outcome SNPs to select
#' @param alpha_init Initial alpha value
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @export
run_cml_family <- function(
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    a = 10,
    b = 10,
    alpha_init = 0.0,
    max_iter = 100,
    tol = 1e-6
) {
    if (abs(a) + abs(b) == 0) stop("a and b cannot both be zero")

    inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
    inv_beta_sigma_out <- lapply(beta_sigma_out, solve)

    run_cml_family_internal(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        inv_beta_sigma_exp = inv_beta_sigma_exp,
        inv_beta_sigma_out = inv_beta_sigma_out,
        beta_sigma_out = beta_sigma_out,
        a = a,
        b = b,
        alpha_init = alpha_init,
        max_iter = max_iter,
        tol = tol
    )
}

#' Variance Calculation for cML Family
#'
#' @keywords internal
calculate_alpha_variance_internal <- function(
    alpha_hat,
    gamma_hat,
    beta_hat_out,
    beta_sigma_out,
    inv_beta_sigma_exp,
    S_f,
    S_t
) {
    J <- nrow(gamma_hat)
    H_matrix <- matrix(0, ncol = 3 * J + 1, nrow = 3 * J + 1)
    d2l_d_alpha2 <- 0
    for (j in 1:J) {
        M_j <- diag(c(
            ifelse(j %in% S_f, 1, 0),
            ifelse(j %in% S_t, 1, 0),
            ifelse(j %in% S_t, 1, 0)
        ))
        gamma_j <- gamma_hat[j, ]
        term2 <- inv_M_j(M_j, beta_sigma_out[[j]])
        d2l_d_alpha2 <- d2l_d_alpha2 + as.numeric(t(gamma_j) %*% term2 %*% gamma_j)
    }
    H_matrix[1, 1] <- d2l_d_alpha2

    for (j in 1:J) {
        gamma_j <- gamma_hat[j, ]
        hat_beta_j <- beta_hat_out[j, ]
        inv_s_gamma <- inv_beta_sigma_exp[[j]]
        M_j <- diag(c(
            ifelse(j %in% S_f, 1, 0),
            ifelse(j %in% S_t, 1, 0),
            ifelse(j %in% S_t, 1, 0)
        ))
        inv_beta_sigma_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
        H_gg <- inv_s_gamma + (alpha_hat^2) * inv_beta_sigma_out_j
        term_common <- inv_beta_sigma_out_j
        H_ga <- -term_common %*% hat_beta_j + 2 * alpha_hat * term_common %*% gamma_j
        idx <- (3 * j - 1):(3 * j + 1)
        H_matrix[idx, 1] <- H_ga
        H_matrix[1, idx] <- t(H_ga)
        H_matrix[idx, idx] <- H_gg
    }

    information <- tryCatch(solve(H_matrix), error = function(e) NULL)

    if (is.null(information)) return(list(variance = NA_real_, std_error = NA_real_))
    variance <- information[1, 1]
    if (!is.finite(variance) || variance < 0) return(list(variance = NA_real_, std_error = NA_real_))
    std_error <- sqrt(variance)
    if (!is.finite(std_error)) return(list(variance = NA_real_, std_error = NA_real_))

    list(variance = variance, std_error = std_error)
}

#' Variance Calculation Version 2
#'
#' @keywords internal
calculate_alpha_variance_internal_two <- function(
    alpha_hat,
    beta_hat_out,
    beta_hat_exp,
    beta_sigma_out,
    beta_sigma_exp,
    S_f,
    S_t
) {
    J <- nrow(beta_hat_out)
    S_combin <- intersect(S_f, S_t)
    S_fdt <- setdiff(S_f, S_t)
    S_tdf <- setdiff(S_t, S_f)
    l_1 <- 0
    l_2 <- 0
    l_3 <- 0

    for (j in 1:J) {
        if (j %in% S_combin) {
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
        if (j %in% S_fdt) {
            beta_hat_out_o <- beta_hat_out[j, 1]
            beta_hat_exp_o <- beta_hat_exp[j, 1]
            u_j <- beta_hat_out_o - alpha_hat * beta_hat_exp_o
            beta_sigma_out_o <- beta_sigma_out[[j]][1, 1]
            beta_sigma_exp_o <- beta_sigma_exp[[j]][1, 1]
            v_j <- beta_sigma_out_o + alpha_hat^2 * beta_sigma_exp_o
            term1 <- beta_hat_exp_o^2 / v_j
            term2 <- (4 * alpha_hat * beta_hat_exp_o * beta_sigma_exp_o * u_j) / v_j^2
            term3 <- (beta_sigma_exp_o * u_j^2) / v_j^2
            term4 <- (4 * alpha_hat^2 * beta_sigma_exp_o^2 * u_j^2) / v_j^3
            l_2 <- term1 + term2 - term3 + term4 + l_2
        }
        if (j %in% S_tdf) {
            beta_sigma_out_j <- beta_sigma_out[[j]][2:3, 2:3]
            beta_sigma_exp_j <- beta_sigma_exp[[j]][2:3, 2:3]
            beta_hat_exp_j <- beta_hat_exp[j, 2:3]
            u_j <- beta_hat_out[j, 2:3] - alpha_hat * beta_hat_exp_j
            Sigma_j_alpha <- beta_sigma_out_j + alpha_hat^2 * beta_sigma_exp_j
            Sigma_j_alpha_inv <- solve(Sigma_j_alpha)
            A_j <- Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv
            term1 <- beta_hat_exp_j %*% Sigma_j_alpha_inv %*% beta_hat_exp_j - u_j %*% A_j %*% u_j
            term2 <- 4 * alpha_hat * beta_hat_exp_j %*% A_j %*% u_j
            term3 <- 2 * alpha_hat^2 * u_j %*% (A_j %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv +
                                                     Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% A_j) %*% u_j
            l_3 <- term1 + term2 + term3 + l_3
        }
    }
    l <- l_1 + l_2 + l_3
    se <- 1 / sqrt(l)
    list(variance = se^2, std_error = se)
}

#' Calculate BIC for cML Family
#'
#' @keywords internal
calculate_bic <- function(
    alpha_hat,
    gamma_hat,
    beta_hat_out,
    beta_hat_exp,
    inv_beta_sigma_exp,
    beta_sigma_out,
    S_f,
    S_t,
    n = 1000,
    p = 1,
    q = 1 / 3 * p
) {
    J <- nrow(beta_hat_out)
    term1 <- 0
    term2 <- 0
    a <- length(S_f)
    b <- length(S_t)
    for (j in 1:J) {
        M_j <- diag(c(
            ifelse(j %in% S_f, 1, 0),
            ifelse(j %in% S_t, 1, 0),
            ifelse(j %in% S_t, 1, 0)
        ))
        inv_beta_sigma_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
        inv_beta_sigma_exp_j <- inv_beta_sigma_exp[[j]]
        term1 <- 1 / 2 * ((beta_hat_exp[j, ] - gamma_hat[j, ]) %*%
                              inv_beta_sigma_exp_j %*%
                              (beta_hat_exp[j, ] - gamma_hat[j, ])) + term1
        term21 <- (beta_hat_out[j, ] - alpha_hat * gamma_hat[j, ])
        term2 <- 1 / 2 * term21 %*% inv_beta_sigma_out_j %*% term21 + term2
    }
    penalty <- log(n) * 2 * J + log(n) * q * (2 * (J - b)) + log(n) * p * ((J - a))
    bic_weight <- penalty + 2 * (term1 + term2)
    bic_weight
}

#' Majority Voting Controller for Multi-Start
#'
#' @keywords internal
cml_majority_controller_ag <- function(
    beta_hat_exp,
    beta_hat_out,
    inv_beta_sigma_exp,
    beta_sigma_out,
    a,
    b,
    n_runs = 5,
    alpha_inits = NULL,
    max_iter = 100,
    tol = 1e-6,
    alpha_tol = NULL,
    gamma_tol = NULL,
    majority_frac = 0.6,
    verbose = FALSE
) {
    if (is.null(alpha_tol)) alpha_tol <- 5 * tol
    if (is.null(gamma_tol)) gamma_tol <- 5 * tol

    if (is.null(alpha_inits)) {
        alpha_inits <- seq(-1, 1, length.out = n_runs)
    } else {
        n_runs <- length(alpha_inits)
    }

    res_list <- vector("list", n_runs)
    AG <- matrix(NA_real_, n_runs, 2L, dimnames = list(NULL, c("alpha", "gamma")))
    ok <- rep(FALSE, n_runs)

    for (i in seq_len(n_runs)) {
        res <- tryCatch(
            run_cml_family_internal(
                beta_hat_exp = beta_hat_exp,
                beta_hat_out = beta_hat_out,
                inv_beta_sigma_exp = inv_beta_sigma_exp,
                beta_sigma_out = beta_sigma_out,
                a = a,
                b = b,
                alpha_init = alpha_inits[i],
                max_iter = max_iter,
                tol = tol
            ),
            error = function(e) NULL
        )
        res_list[[i]] <- res
        a_i <- suppressWarnings(as.numeric(if (!is.null(res)) res$alpha[1] else NA_real_))
        g_i <- suppressWarnings(as.numeric(if (!is.null(res)) res$gamma[1] else NA_real_))
        if (is.finite(a_i) && is.finite(g_i)) {
            AG[i, ] <- c(a_i, g_i)
            ok[i] <- TRUE
        }
    }

    if (!any(ok)) return(NA)

    AG_ok <- AG[ok, , drop = FALSE]
    idx_ok <- which(ok)
    m <- nrow(AG_ok)

    if (m == 1L) {
        if (1 >= ceiling(majority_frac * 1)) {
            return(res_list[[idx_ok[1L]]])
        }
        return(NA)
    }

    sim <- matrix(FALSE, m, m)
    for (i in seq_len(m)) {
        for (j in i:m) {
            same_cluster <- (abs(AG_ok[i, 1] - AG_ok[j, 1]) <= alpha_tol) &&
                (abs(AG_ok[i, 2] - AG_ok[j, 2]) <= gamma_tol)
            sim[i, j] <- same_cluster
            sim[j, i] <- same_cluster
        }
    }

    labels <- rep(NA_integer_, m)
    lab <- 0L
    for (i in seq_len(m)) {
        if (is.na(labels[i])) {
            lab <- lab + 1L
            q <- i
            labels[i] <- lab
            while (length(q)) {
                k <- q[1]
                q <- q[-1]
                nbrs <- which(sim[k, ] & is.na(labels))
                labels[nbrs] <- lab
                q <- c(q, nbrs)
            }
        }
    }

    tab <- sort(table(labels), decreasing = TRUE)
    best_lab <- as.integer(names(tab)[1L])
    size_best <- as.integer(tab[1L])
    prop_best <- size_best / m

    if (prop_best < majority_frac) {
        if (verbose) message(sprintf("no majority cluster (%.1f%% < %.1f%%)", 100 * prop_best, 100 * majority_frac))
        return(NA)
    }

    idx_in <- which(labels == best_lab)
    idx_runs_in <- idx_ok[idx_in]
    center <- colMeans(AG_ok[idx_in, , drop = FALSE])
    dists <- apply(AG_ok[idx_in, , drop = FALSE], 1, function(r) sum(abs(r - center)))
    rep_local <- which.min(dists)
    rep_run <- idx_runs_in[rep_local]

    res_list[[rep_run]]
}

#' BIC-Weighted Controller for Multi-Start
#'
#' @keywords internal
cml_bic_weighted_controller_ag <- function(
    beta_hat_exp,
    beta_hat_out,
    inv_beta_sigma_exp,
    beta_sigma_out,
    a,
    b,
    n_runs = 5,
    alpha_inits = NULL,
    max_iter = 100,
    tol = 1e-6,
    alpha_tol = NULL,
    gamma_tol = NULL,
    q = 1 / 3,
    n = 1000,
    p = 1
) {
    J <- nrow(beta_hat_out)
    cand_list <- vector("list", n_runs)
    bic_vec <- rep(NA_real_, n_runs)

    for (r in seq_len(n_runs)) {
        if (is.null(alpha_inits)) {
            alpha_inits <- seq(-1, 1, length.out = n_runs)
        } else {
            n_runs <- length(alpha_inits)
        }
        fit_r <- tryCatch(
            run_cml_family_internal(
                beta_hat_exp = beta_hat_exp,
                beta_hat_out = beta_hat_out,
                inv_beta_sigma_exp = inv_beta_sigma_exp,
                beta_sigma_out = beta_sigma_out,
                a = a,
                b = b,
                alpha_init = alpha_inits[r],
                max_iter = max_iter,
                tol = tol
            ),
            error = function(e) NULL
        )

        if (is.null(fit_r) || !is.list(fit_r) || is.null(fit_r$alpha) || is.null(fit_r$gamma)) {
            cand_list[[r]] <- NULL
            bic_vec[r] <- NA_real_
            next
        }

        bic_r <- tryCatch(
            calculate_bic(
                alpha_hat = fit_r$alpha,
                gamma_hat = fit_r$gamma,
                beta_hat_out = beta_hat_out,
                beta_hat_exp = beta_hat_exp,
                inv_beta_sigma_exp = inv_beta_sigma_exp,
                beta_sigma_out = beta_sigma_out,
                S_f = fit_r$S_f,
                S_t = fit_r$S_t,
                n = n,
                p = p,
                q = q
            ),
            error = function(e) NA_real_
        )

        cand_list[[r]] <- fit_r
        bic_vec[r] <- bic_r
    }

    ok <- which(is.finite(bic_vec))
    if (length(ok) == 0L) return(NULL)

    bic_ok <- bic_vec[ok]
    min_bic <- min(bic_ok)
    delta_bic <- bic_ok - min_bic
    exp_term <- exp(-0.5 * delta_bic)
    denom <- sum(exp_term)

    if (!is.finite(denom) || denom <= 0) {
        warning("BIC weight normalization issue.")
        return(NULL)
    }

    weights <- as.numeric(exp_term / denom)

    alpha_w <- 0
    gamma_w <- matrix(0, nrow = J, ncol = 3)
    for (i in seq_along(ok)) {
        r <- ok[i]
        wi <- weights[i]
        alpha_w <- alpha_w + wi * as.numeric(cand_list[[r]]$alpha[1])
        gamma_w <- gamma_w + wi * as.matrix(cand_list[[r]]$gamma)
    }

    best_idx <- ok[which.min(bic_ok)]
    Sf_best <- cand_list[[best_idx]]$S_f
    St_best <- cand_list[[best_idx]]$S_t

    out <- list(
        alpha = alpha_w,
        gamma = gamma_w,
        S_f = Sf_best,
        S_t = St_best,
        bic_vec = bic_vec,
        weights = setNames(rep(NA_real_, length(bic_vec)), paste0("run", seq_along(bic_vec)))
    )
    out$weights[ok] <- weights
    out
}

#' Fast Grid Search over (a, b) Parameter Space
#'
#' @keywords internal
tune_cml_family_grid_fast <- function(
    a_values,
    b_values,
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    n = 1000,
    p = 1,
    q = 1 / 3,
    inv_beta_sigma_exp_precomp = NULL,
    controller = cml_bic_weighted_controller_ag
) {
    inv_beta_sigma_exp <- if (is.null(inv_beta_sigma_exp_precomp)) {
        lapply(beta_sigma_exp, solve)
    } else {
        inv_beta_sigma_exp_precomp
    }

    na <- length(a_values)
    nb <- length(b_values)
    n_pairs <- na * nb
    results_list <- vector("list", n_pairs)
    idx <- 1L

    for (ia in seq_len(na)) {
        a_cur <- a_values[ia]
        for (ib in seq_len(nb)) {
            b_cur <- b_values[ib]

            est_results <- tryCatch(
                cml_bic_weighted_controller_ag(
                    beta_hat_exp = beta_hat_exp,
                    beta_hat_out = beta_hat_out,
                    inv_beta_sigma_exp = inv_beta_sigma_exp,
                    beta_sigma_out = beta_sigma_out,
                    a = a_cur,
                    b = b_cur,
                    n_runs = 5,
                    alpha_inits = NULL,
                    max_iter = 100,
                    tol = 1e-06,
                    alpha_tol = NULL,
                    gamma_tol = NULL,
                    q = q,
                    n = n,
                    p = p
                ),
                error = function(e) NULL
            )

            if (is.null(est_results) || !is.list(est_results)) {
                results_list[[idx]] <- data.frame(
                    a = a_cur, b = b_cur,
                    alpha_hat = NA_real_,
                    variance = NA_real_,
                    std_error = NA_real_,
                    p_value = NA_real_,
                    bic = NA_real_,
                    S_f = I(list(NA_integer_)),
                    S_t = I(list(NA_integer_))
                )
                idx <- idx + 1L
                next
            }

            var_results <- tryCatch(
                calculate_alpha_variance_internal_two(
                    alpha_hat = est_results$alpha,
                    beta_hat_out = beta_hat_out,
                    beta_hat_exp = beta_hat_exp,
                    beta_sigma_out = beta_sigma_out,
                    beta_sigma_exp = beta_sigma_exp,
                    S_f = est_results$S_f,
                    S_t = est_results$S_t
                ),
                error = function(e) list(variance = NA_real_, std_error = NA_real_)
            )

            alpha_hat_val <- suppressWarnings(as.numeric(est_results$alpha[1]))
            se_hat_val <- suppressWarnings(as.numeric(var_results$std_error[1]))
            if (!is.finite(alpha_hat_val)) alpha_hat_val <- NA_real_
            if (!is.finite(se_hat_val)) se_hat_val <- NA_real_

            p_val <- if (is.finite(alpha_hat_val) && is.finite(se_hat_val) && se_hat_val > 0) {
                2 * pnorm(abs(alpha_hat_val / se_hat_val), lower.tail = FALSE)
            } else {
                NA_real_
            }

            bic_val <- tryCatch(
                calculate_bic(
                    alpha_hat = est_results$alpha,
                    gamma_hat = est_results$gamma,
                    beta_hat_out = beta_hat_out,
                    beta_hat_exp = beta_hat_exp,
                    inv_beta_sigma_exp = inv_beta_sigma_exp,
                    beta_sigma_out = beta_sigma_out,
                    S_f = est_results$S_f,
                    S_t = est_results$S_t,
                    n = n,
                    p = p,
                    q = q
                ),
                error = function(e) NA_real_
            )

            results_list[[idx]] <- data.frame(
                a = a_cur,
                b = b_cur,
                alpha_hat = alpha_hat_val,
                variance = suppressWarnings(as.numeric(var_results$variance[1])),
                std_error = se_hat_val,
                p_value = p_val,
                bic = bic_val,
                S_f = I(list(est_results$S_f)),
                S_t = I(list(est_results$S_t))
            )

            idx <- idx + 1L
        }
    }

    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    results_df
}

#' Calculate BIC-Weighted Estimates (BMA)
#'
#' @keywords internal
calculate_bic_weighted_estimates_fast <- function(results_df) {
    stopifnot(is.data.frame(results_df))

    valid_mask <- with(
        results_df,
        !is.na(variance) & is.finite(variance) &
            !is.na(bic) & is.finite(bic) &
            !is.na(alpha_hat) & is.finite(alpha_hat)
    )

    results_with_weights <- results_df
    results_with_weights$weight <- NA_real_

    if (!any(valid_mask)) {
        return(list(
            alpha_bma = NA_real_,
            se_bma = NA_real_,
            results_with_weights = results_with_weights
        ))
    }

    vr <- results_df[valid_mask, , drop = FALSE]
    min_bic <- min(vr$bic)
    delta_bic <- vr$bic - min_bic
    exp_term <- exp(-0.5 * delta_bic)
    denom <- sum(exp_term)

    if (!is.finite(denom) || denom <= 0) {
        warning("Numerical issue computing BIC weights.")
        return(list(
            alpha_bma = NA_real_,
            se_bma = NA_real_,
            results_with_weights = results_with_weights
        ))
    }

    weights <- exp_term / denom
    results_with_weights$weight[valid_mask] <- weights

    alpha_bma <- sum(weights * vr$alpha_hat)
    within_model_var <- sum(weights * vr$variance)
    between_model_var <- sum(weights * (vr$alpha_hat - alpha_bma)^2)
    total_var <- within_model_var + between_model_var
    if (!is.finite(total_var) || total_var < 0) total_var <- NA_real_
    se_bma <- if (is.na(total_var)) NA_real_ else sqrt(total_var)

    list(
        alpha_bma = alpha_bma,
        se_bma = se_bma,
        results_with_weights = results_with_weights
    )
}

#' Resolve (a, b) Grid Parameters
#'
#' @keywords internal
.resolve_family_ab_grid <- function(a_values, b_values, n_snps) {
    if (missing(a_values) || is.null(a_values)) {
        a_values <- seq_len(n_snps)
    }
    if (missing(b_values) || is.null(b_values)) {
        b_values <- 0:n_snps
    }
    list(
        a_values = as.integer(a_values),
        b_values = as.integer(b_values)
    )
}

#' MRcML-Family Main Analysis Function
#'
#' Performs Mendelian Randomization analysis using cML family methods with
#' BIC-based Bayesian Model Averaging across (a, b) parameter grid.
#'
#' @param a_values Vector of a parameters (horizontal pleiotropy). Default: 1:n_snps
#' @param b_values Vector of b parameters (vertical pleiotropy). Default: 0:n_snps
#' @param beta_hat_exp Exposure beta estimates (J x 3 matrix)
#' @param beta_hat_out Outcome beta estimates (J x 3 matrix)
#' @param beta_sigma_exp Exposure covariance matrices (list of 3x3 matrices)
#' @param beta_sigma_out Outcome covariance matrices (list of 3x3 matrices)
#' @param n Sample size (default 1000)
#' @param p Number of exposures (default 1)
#' @param q Ratio parameter (default 1/3)
#' @return List with alpha_bma (BMA estimate), se_bma (BMA standard error),
#'   and results_with_weights (full grid results)
#' @export
MRcML_family <- function(
    a_values,
    b_values,
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    n = 1000,
    p = 1,
    q = 1 / 3
) {
    n_snps <- nrow(beta_hat_exp)
    grid <- .resolve_family_ab_grid(a_values, b_values, n_snps)
    inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)

    grid_fit <- tune_cml_family_grid_fast(
        a_values = grid$a_values,
        b_values = grid$b_values,
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        beta_sigma_exp = beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        n = n,
        p = p,
        q = q,
        inv_beta_sigma_exp_precomp = inv_beta_sigma_exp
    )

    bic_weighted <- calculate_bic_weighted_estimates_fast(grid_fit)
    bic_weighted
}

#' Alias for MRcML_family
#'
#' @rdname MRcML_family
#' @export
cml_family_new <- function(
    a_values,
    b_values,
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    n = 1000,
    p = 1,
    q = 1 / 3
) {
    MRcML_family(
        a_values = a_values,
        b_values = b_values,
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        beta_sigma_exp = beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        n = n,
        p = p,
        q = q
    )
}
