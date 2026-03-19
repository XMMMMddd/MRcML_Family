# %% 载入要用到的包
# Rcpp::sourceCpp("cML-family/cML-family-new/run_cml_family_internal_rcpp.cpp")
Rcpp::sourceCpp("cML-family/cML-family-new/alpha_variance.cpp")
# Rcpp::sourceCpp("cML-family/cML-family-new/run_cml_family_internal_rcpp.cpp")
Rcpp::sourceCpp("cML-family/cML-family-new/MRcML_family_MA.cpp")
# %% 定义函数
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
  return(inv_beta_sigma = inv_results)
}

run_cml_family_internal <- function(
  beta_hat_exp, beta_hat_out,
  inv_beta_sigma_exp,
  beta_sigma_out, # 仍然需要原始矩阵来提取对角元和子矩阵
  a, b, alpha_init, max_iter, tol
) {
  J <- nrow(beta_hat_exp)
  alpha_k <- alpha_init
  gamma_k <- beta_hat_exp

  for (k in 1:max_iter) {
    # --- Step 1 ---
    f_values <- numeric(J)
    t_values <- numeric(J)
    for (j in 1:J) {
      sigma_beta_o_sq <- beta_sigma_out[[j]][1, 1]
      f_values[j] <- (beta_hat_out[j, 1] -
        alpha_k * gamma_k[j, 1])^2 / sigma_beta_o_sq
      d_fm <- beta_hat_out[j, 2:3] - alpha_k * gamma_k[j, 2:3]
      sigma_beta_fm <- beta_sigma_out[[j]][2:3, 2:3]
      inv_sigma_beta_fm <- solve(sigma_beta_fm)
      t_values[j] <- as.numeric(t(d_fm) %*% inv_sigma_beta_fm %*% d_fm)
    }
    sf_indices <- order(f_values)
    st_indices <- order(t_values)
    S_f <- sf_indices[1:a]
    S_t <- st_indices[1:b]
    if (a == 0) {
      S_f <- NULL
    }
    if (b == 0) {
      S_t <- NULL
    }

    # --- Step 2 ---
    gamma_k_plus_1 <- matrix(0, nrow = J, ncol = 3)
    for (j in 1:J) {
      M_j <- diag(c(
        ifelse(j %in% S_f, 1, 0),
        ifelse(j %in% S_t, 1, 0),
        ifelse(j %in% S_t, 1, 0)
      ))
      inv_beta_sigma_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
      term1_inv <- inv_beta_sigma_exp[[j]] +
        (alpha_k^2) * inv_beta_sigma_out_j
      term1 <- solve(term1_inv)
      term2 <- (inv_beta_sigma_exp[[j]] %*% beta_hat_exp[j,]) +
        (alpha_k *
          (inv_beta_sigma_out_j
            %*% beta_hat_out[j,]))
      gamma_k_plus_1[j,] <- term1 %*% term2
    }

    # --- Step 3 ---
    numerator <- 0.0
    denominator <- 0.0
    for (j in 1:J) {
      M_j <- diag(c(
        ifelse(j %in% S_f, 1, 0),
        ifelse(j %in% S_t, 1, 0),
        ifelse(j %in% S_t, 1, 0)
      ))
      gamma_j <- gamma_k_plus_1[j,]
      inv_beta_simga_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
      term_middle <- inv_beta_simga_out_j
      numerator <- numerator + as.numeric(t(gamma_j) %*%
                                            term_middle %*%
                                            beta_hat_out[j,])
      denominator <- denominator + as.numeric(t(gamma_j) %*%
                                                term_middle %*%
                                                gamma_j)
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

  return(list(alpha = alpha_k, gamma = gamma_k, S_f = S_f, S_t = S_t))
}

run_cml_family <- function(
  beta_hat_exp,
  beta_hat_out, beta_sigma_exp, beta_sigma_out,
  a = 10, b = 10, alpha_init = 0.0,
  max_iter = 100, tol = 1e-6
) {
  if (abs(a) + abs(b) == 0) {
    stop()
  }

  # 在这里进行一次性求逆
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
  inv_beta_sigma_out <- lapply(beta_sigma_out, solve)

  # 调用内部函数进行计算
  run_cml_family_internal(
    beta_hat_exp = beta_hat_exp, beta_hat_out = beta_hat_out,
    inv_beta_sigma_exp = inv_beta_sigma_exp,
    inv_beta_sigma_out = inv_beta_sigma_out,
    beta_sigma_out = beta_sigma_out,
    a = a, b = b, alpha_init = alpha_init, max_iter = max_iter, tol = tol
  )
}

run_cml_family_rcpp <- function(beta_hat_exp,
                                beta_hat_out,
                                beta_sigma_exp, beta_sigma_out,
                                a = 10, b = 10, alpha_init = 0.0, max_iter = 100, tol = 1e-6) {
  if (a < 0 || b < 0 || (a == 0 && b == 0)) {
    stop("a 和 b 必须为非负整数, 且至少一个不为零。")
  }

  # 一次性求逆 (仅对 inv_beta_sigma_exp)
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)

  # 调用 Rcpp 内部函数
  # 注意：Rcpp 函数不需要 inv_beta_sigma_out，所以不传递它
  run_cml_family_internal_rcpp(
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    inv_beta_sigma_exp = inv_beta_sigma_exp,
    beta_sigma_out = beta_sigma_out,
    a = a,
    b = b,
    alpha_init = alpha_init,
    max_iter = max_iter,
    tol = tol
  )
}

calculate_alpha_variance_internal <- function(
  alpha_hat, gamma_hat,
  beta_hat_out,
  beta_sigma_out, inv_beta_sigma_exp,
  S_f, S_t
) {
  J <- nrow(gamma_hat)
  H_matrix <- matrix(0, ncol = 3 * J + 1, nrow = 3 * J + 1)

  # d2l/dalpha^2
  d2l_d_alpha2 <- 0
  for (j in 1:J) {
    M_j <- diag(c(
      ifelse(j %in% S_f, 1, 0),
      ifelse(j %in% S_t, 1, 0),
      ifelse(j %in% S_t, 1, 0)
    ))
    gamma_j <- gamma_hat[j,]
    term2 <- inv_M_j(M_j, beta_sigma_out[[j]])

    d2l_d_alpha2 <- d2l_d_alpha2 + as.numeric(t(gamma_j) %*% term2 %*% gamma_j)
  }
  H_matrix[1, 1] <- d2l_d_alpha2

  # 组装 Hessian 其余块
  for (j in 1:J) {
    gamma_j <- gamma_hat[j,]
    hat_beta_j <- beta_hat_out[j,]
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

  # --- 稳健增强：solve & 方差检查 ---
  information <- tryCatch(solve(H_matrix),
                          error = function(e) NULL
  )

  if (is.null(information)) {
    # Hessian 不可逆 → 返回 NA
    return(list(variance = NA_real_, std_error = NA_real_))
  }

  variance <- information[1, 1]

  # 若方差为负或非有限 → 返回 NA
  if (!is.finite(variance) || variance < 0) {
    return(list(variance = NA_real_, std_error = NA_real_))
  }

  std_error <- sqrt(variance)
  if (!is.finite(std_error)) {
    return(list(variance = NA_real_, std_error = NA_real_))
  }

  list(variance = variance, std_error = std_error)
}

calculate_alpha_variance_internal_two <- function(
  alpha_hat,
  beta_hat_out, beta_hat_exp,
  beta_sigma_out, beta_sigma_exp,
  S_f, S_t
) {
  J <- nrow(beta_hat_out)
  S_combin <- intersect(S_f, S_t)
  S_fdt <- setdiff(S_f, S_t)
  S_tdf <- setdiff(S_t, S_f)
  l_1 <- 0
  l_2 <- 0
  l_3 <- 0
  for (j in 1:J) {
    # 第一种情况
    if (j %in% S_combin) {
      beta_sigma_out_j <- beta_sigma_out[[j]]
      beta_sigma_exp_j <- beta_sigma_exp[[j]]
      beta_hat_exp_j <- beta_hat_exp[j,]
      u_j <- beta_hat_out[j,] - alpha_hat * beta_hat_exp_j
      Sigma_j_alpha <- beta_sigma_out_j + alpha_hat^2 * beta_sigma_exp_j
      Sigma_j_alpha_inv <- solve(Sigma_j_alpha)

      A_j <- Sigma_j_alpha_inv %*%
        beta_sigma_exp_j %*%
        Sigma_j_alpha_inv

      term1 <- beta_hat_exp_j %*%
        Sigma_j_alpha_inv %*%
        beta_hat_exp_j - u_j %*% A_j %*% u_j
      term2 <- 4 * alpha_hat * beta_hat_exp_j %*% A_j %*% u_j
      term3 <- 2 *
        alpha_hat^2 *
        u_j %*%
          (A_j %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv + Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% A_j) %*%
          u_j
      l_1 <- term1 + term2 + term3 + l_1
    }
    # 第二种情况
    if (j %in% S_fdt) {
      beta_hat_out_o <- beta_hat_out[j, 1]
      beta_hat_exp_o <- beta_hat_exp[j, 1]
      u_j <- beta_hat_out_o - alpha_hat * beta_hat_exp_o
      beta_sigma_out_o <- beta_sigma_out[[j]][1, 1]
      beta_sigma_exp_o <- beta_sigma_exp[[j]][1, 1]
      v_j <- beta_sigma_out_o + alpha_hat^2 * beta_sigma_exp_o

      term1 <- beta_hat_exp_o^2 / v_j
      term2 <- (4 *
        alpha_hat *
        beta_hat_exp_o *
        beta_sigma_exp_o *
        u_j) / v_j^2
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

      A_j <- Sigma_j_alpha_inv %*%
        beta_sigma_exp_j %*%
        Sigma_j_alpha_inv

      term1 <- beta_hat_exp_j %*%
        Sigma_j_alpha_inv %*%
        beta_hat_exp_j - u_j %*% A_j %*% u_j
      term2 <- 4 * alpha_hat * beta_hat_exp_j %*% A_j %*% u_j
      term3 <- 2 *
        alpha_hat^2 *
        u_j %*%
          (A_j %*% beta_sigma_exp_j %*% Sigma_j_alpha_inv + Sigma_j_alpha_inv %*% beta_sigma_exp_j %*% A_j) %*%
          u_j
      l_3 <- term1 + term2 + term3 + l_3
    }
  }
  l <- l_1 + l_2 + l_3
  l
  se <- 1 / sqrt(l)
  return(list(variance = se^2, std_error = se))
}

calculate_alpha_variance <- function(
  alpha_hat, gamma_hat, beta_hat_out,
  beta_sigma_out, beta_sigma_exp, S_f, S_t
) {
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
  calculate_alpha_variance_internal(
    alpha_hat = alpha_hat, gamma_hat = gamma_hat,
    beta_hat_out = beta_hat_out, beta_sigma_out = beta_sigma_out,
    inv_beta_sigma_exp = inv_beta_sigma_exp,
    S_f = S_f, S_t = S_t
  )
}


calculate_bic <- function(
  alpha_hat, gamma_hat,
  beta_hat_out, beta_hat_exp,
  inv_beta_sigma_exp, # <--- 接收预计算的逆矩阵
  beta_sigma_out,
  S_f, S_t, n = 1000, p = 1, q = 1 / 3 * p
) {
  J <- nrow(beta_hat_out)
  term1 <- 0
  term2 <- 0
  a <- length(S_f)
  b <- length(S_t)
  for (j in 1:J) {
    M_j <- diag(c(
      ifelse(j %in% S_f, 1, 0),
      ifelse(j %in% S_t, 1, 0), ifelse(j %in% S_t, 1, 0)
    ))
    inv_beta_sigma_out_j <- inv_M_j(M_j, beta_sigma_out[[j]])
    inv_beta_sigma_exp_j <- inv_beta_sigma_exp[[j]]
    term1 <- 1 / 2 * ((beta_hat_exp[j,] - gamma_hat[j,]) %*%
      inv_beta_sigma_exp_j %*%
      (beta_hat_exp[j,] - gamma_hat[j,])) + term1
    term21 <- (beta_hat_out[j,] - alpha_hat * gamma_hat[j,])
    term2 <- 1 / 2 * term21 %*% inv_beta_sigma_out_j %*% term21 + term2
  }
  penalty <- log(n) * 2 * J +
    log(n) * q * (2 * (J - b)) +
    log(n) * p * ((J - a))
  bic_weight <- penalty + 2 * (term1 + term2)
  return(bic_weight)
}

# 多起点重启 + 多数一致性，仅比较 alpha 与 gamma
cml_majority_controller_ag <- function(
  beta_hat_exp, beta_hat_out, inv_beta_sigma_exp, beta_sigma_out,
  a, b,
  n_runs = 5,
  alpha_inits = NULL, # 若为 NULL，将自动生成等距起点
  max_iter = 100,
  tol = 1e-6, # 传给单次迭代函数
  alpha_tol = NULL, # 判定同一解用的容差；默认 5*tol
  gamma_tol = NULL, # 判定同一解用的容差；默认 5*tol
  majority_frac = 0.6, # “大多数”的阈值（>=0.5 合理）
  verbose = FALSE
) {
  # 默认容差：略大于算法内部 tol
  if (is.null(alpha_tol)) alpha_tol <- 5 * tol
  if (is.null(gamma_tol)) gamma_tol <- 5 * tol

  # 生成多起点
  if (is.null(alpha_inits)) {
    alpha_inits <- seq(-1, 1, length.out = n_runs)
  } else {
    n_runs <- length(alpha_inits)
  }

  res_list <- vector("list", n_runs)
  AG <- matrix(NA_real_, n_runs, 2L, dimnames = list(NULL, c("alpha", "gamma")))
  ok <- rep(FALSE, n_runs)

  # 多次运行单次迭代
  for (i in seq_len(n_runs)) {
    res <- tryCatch(
      run_cml_family_internal_rcpp(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        inv_beta_sigma_exp = inv_beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        a = a, b = b,
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
      AG[i,] <- c(a_i, g_i)
      ok[i] <- TRUE
    } else if (verbose) {
      message(sprintf("run %d failed or produced non-finite alpha/gamma", i))
    }
  }

  # 没有任何有效结果
  if (!any(ok)) {
    return(NA)
  }

  AG_ok <- AG[ok, , drop = FALSE]
  idx_ok <- which(ok)
  m <- nrow(AG_ok)

  # 单个有效结果：只有在门槛允许时才算“多数”
  if (m == 1L) {
    if (1 >= ceiling(majority_frac * 1)) {
      return(res_list[[idx_ok[1L]]])
    }
    return(NA)
  }

  # 两两相似（属于同一簇）的判定：矩形阈值 |Δalpha|≤alpha_tol & |Δgamma|≤gamma_tol
  sim <- matrix(FALSE, m, m)
  for (i in seq_len(m)) {
    for (j in i:m) {
      same_cluster <- (abs(AG_ok[i, 1] - AG_ok[j, 1]) <= alpha_tol) &&
        (abs(AG_ok[i, 2] - AG_ok[j, 2]) <= gamma_tol)
      sim[i, j] <- same_cluster
      sim[j, i] <- same_cluster
    }
  }

  # 连通分量聚类（BFS）
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
        nbrs <- which(sim[k,] & is.na(labels))
        labels[nbrs] <- lab
        q <- c(q, nbrs)
      }
    }
  }

  # 取最大簇并检查是否达到“多数”
  tab <- sort(table(labels), decreasing = TRUE)
  best_lab <- as.integer(names(tab)[1L])
  size_best <- as.integer(tab[1L])
  prop_best <- size_best / m

  if (prop_best < majority_frac) {
    if (verbose) {
      message(sprintf(
        "no majority cluster (%.1f%% < %.1f%%)",
        100 * prop_best, 100 * majority_frac
      ))
    }
    return(NA)
  }

  # 从多数簇中选代表：对 (alpha,gamma) 与簇中心的 L1 距离最小的那个 run
  idx_in <- which(labels == best_lab)
  idx_runs_in <- idx_ok[idx_in]
  center <- colMeans(AG_ok[idx_in, , drop = FALSE])
  dists <- apply(
    AG_ok[idx_in, , drop = FALSE], 1,
    function(r) sum(abs(r - center))
  )
  rep_local <- which.min(dists)
  rep_run <- idx_runs_in[rep_local]

  # 返回：代表 run 的完整结果（结构与单次函数一致：含 alpha/gamma/S_f/S_t）
  res_list[[rep_run]]
}

cml_bic_weighted_controller_ag <- function(
  beta_hat_exp,
  beta_hat_out,
  inv_beta_sigma_exp,
  beta_sigma_out,
  a, b,
  n_runs = 5,
  alpha_inits = NULL, # 若为 NULL，将自动生成等距起点
  max_iter = 100,
  tol = 1e-6, # 传给单次迭代函数
  alpha_tol = NULL, # 判定同一解用的容差；默认 5*tol
  gamma_tol = NULL, # 判定同一解用的容差；默认 5*tol
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
      run_cml_family_internal_rcpp(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        inv_beta_sigma_exp = inv_beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        a = a, b = b,
        alpha_init = alpha_inits[r], max_iter = max_iter, tol = tol
      ),
      error = function(e) NULL
    )

    if (is.null(fit_r) ||
      !is.list(fit_r) ||
      is.null(fit_r$alpha) ||
      is.null(fit_r$gamma) ||
      is.null(fit_r$S_f) ||
      is.null(fit_r$S_t)) {
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
        n = n, p = p, q = q
      ),
      error = function(e) NA_real_
    )

    cand_list[[r]] <- fit_r
    bic_vec[r] <- bic_r
  }

  ok <- which(is.finite(bic_vec))
  if (length(ok) == 0L) {
    return(NULL)
  }

  # —— 按你给的公式计算权重 —— #
  bic_ok <- bic_vec[ok]
  min_bic <- min(bic_ok)
  delta_bic <- bic_ok - min_bic
  exp_term <- exp(-0.5 * delta_bic)
  denom <- sum(exp_term)

  if (!is.finite(denom) || denom <= 0) {
    warning("BIC 权重归一化出现数值问题。")
    return(NULL)
  }

  weights <- as.numeric(exp_term / denom)

  # 加权 alpha、gamma（集合仍取最优 BIC 那组，保证后续方差/矩阵构造一致）
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


# 不调用多起点函数
tune_cml_family_ab_forloop <- function(
  a_values,
  b_values,
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1 / 3 * 1
) {
  # --- 关键优化：在所有循环开始前，只计算一次逆矩阵 ---
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
  inv_beta_sigma_out <- lapply(beta_sigma_out, solve)

  # 创建参数网格
  param_grid <- expand.grid(a = a_values, b = b_values)


  # --- 1. 初始化一个列表来存储结果 ---
  # 预先分配列表大小可以略微提高效率
  results_list <- vector("list", nrow(param_grid))

  # --- 2. 使用标准的 for 循环进行迭代 ---
  for (i in 1:nrow(param_grid)) {
    current_a <- param_grid$a[i]
    current_b <- param_grid$b[i]

    # 调用内部函数，传入预计算好的逆矩阵

    est_results <- run_cml_family_internal_rcpp(
      beta_hat_exp = beta_hat_exp,
      beta_hat_out = beta_hat_out,
      inv_beta_sigma_exp = inv_beta_sigma_exp,
      beta_sigma_out = beta_sigma_out,
      a = current_a, b = current_b,
      alpha_init = 0.0, max_iter = 100, tol = 1e-6
    )
    # 调用内部函数，传入预计算好的逆矩阵
    var_results <- calculate_alpha_variance_internal_two(
      alpha_hat = est_results$alpha,
      beta_hat_out = beta_hat_out,
      beta_hat_exp = beta_hat_exp,
      beta_sigma_out = beta_sigma_out,
      beta_sigma_exp = beta_sigma_exp,
      S_f = est_results$S_f, S_t = est_results$S_t
    )
    p_value <- 2 * pnorm(abs(est_results$alpha /
                               var_results$std_error), lower.tail = FALSE)

    bic <- calculate_bic(
      alpha_hat = est_results$alpha,
      gamma_hat = est_results$gamma,
      beta_hat_out = beta_hat_out,
      beta_hat_exp = beta_hat_exp,
      inv_beta_sigma_exp = inv_beta_sigma_exp,
      beta_sigma_out = beta_sigma_out,
      S_f = est_results$S_f, S_t = est_results$S_t,
      n = 1000, p = p, q = q
    )

    # 创建一个单行数据框作为本次循环的结果
    result_row <- data.frame(
      a = current_a, b = current_b,
      alpha_hat = est_results$alpha,
      variance = var_results$variance,
      std_error = var_results$std_error,
      p_value = p_value, bic = bic
    )

    # 将本次循环的结果存入列表
    results_list[[i]] <- result_row
  }

  # --- 3. 在循环结束后，将列表中的所有数据框高效地合并 ---
  results_df <- do.call(rbind, results_list)

  return(results_df)
}

# 依赖：cml_majority_controller_ag()
# 若你使用不同名字的控制函数，把 controller 替换掉即可
tune_cml_family_ab_forloop_robust <- function(
  a_values,
  b_values,
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1 / 3,
  # 控制函数与其参数：
  controller = cml_majority_controller_ag,
  controller_args = list(
    n_runs = 5,
    alpha_inits = NULL,
    max_iter = 100,
    tol = 1e-6,
    alpha_tol = NULL,
    gamma_tol = NULL,
    majority_frac = 0.8,
    verbose = FALSE
  )
) {
  # --- 关键优化：在所有循环开始前，只计算一次逆矩阵 ---
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)
  inv_beta_sigma_out <- lapply(beta_sigma_out, solve) # 目前未用到，保留以便扩展

  # 参数网格
  param_grid <- expand.grid(a = a_values, b = b_values)

  # 结果列表
  results_list <- vector("list", nrow(param_grid))

  # 工具函数：生成一行全 NA 的占位结果
  na_row <- function(a, b) {
    data.frame(
      a = a, b = b,
      alpha_hat = NA_real_,
      variance = NA_real_,
      std_error = NA_real_,
      p_value = NA_real_,
      bic = NA_real_
    )
  }

  # 循环
  for (i in seq_len(nrow(param_grid))) {
    current_a <- param_grid$a[i]
    current_b <- param_grid$b[i]

    # --- 用控制函数进行稳健拟合（可能返回 NA） ---
    est_results <- tryCatch(
      do.call(controller, c(list(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        inv_beta_sigma_exp = inv_beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        a = current_a, b = current_b
      ), controller_args)),
      error = function(e) NA
    )

    # 控制函数失败：直接写入 NA 行并继续
    if (is.null(est_results) ||
      (is.atomic(est_results) &&
        length(est_results) == 1L &&
        is.na(est_results)) ||
      !is.list(est_results)) {
      results_list[[i]] <- na_row(current_a, current_b)
      next
    }

    # --- 计算方差与标准误（出错则 NA） ---
    var_results <- tryCatch(
      calculate_alpha_variance_internal_two_cpp(
        alpha_hat = est_results$alpha,
        beta_hat_out = beta_hat_out,
        beta_hat_exp = beta_hat_exp,
        beta_sigma_out = beta_sigma_out,
        beta_sigma_exp = beta_sigma_exp,
        S_f = est_results$S_f, S_t = est_results$S_t
      ),
      error = function(e) list(variance = NA_real_, std_error = NA_real_)
    )

    # 提取标量 alpha 与 se（若是向量，取第一个；非有限值则 NA）
    alpha_hat <- suppressWarnings(as.numeric(est_results$alpha[1]))
    se <- suppressWarnings(as.numeric(var_results$std_error[1]))
    if (!is.finite(alpha_hat)) alpha_hat <- NA_real_
    if (!is.finite(se)) se <- NA_real_

    # p 值：se>0 且数值有限才计算
    p_value <- if (is.finite(alpha_hat) && is.finite(se) && se > 0) {
      2 * pnorm(abs(alpha_hat / se), lower.tail = FALSE)
    } else {
      NA_real_
    }

    # --- 计算 BIC（出错则 NA）---
    bic_val <- tryCatch(
      calculate_bic(
        alpha_hat = est_results$alpha,
        gamma_hat = est_results$gamma,
        beta_hat_out = beta_hat_out,
        beta_hat_exp = beta_hat_exp,
        inv_beta_sigma_exp = inv_beta_sigma_exp,
        beta_sigma_out = beta_sigma_out,
        S_f = est_results$S_f, S_t = est_results$S_t,
        n = n, p = p, q = q
      ),
      error = function(e) NA_real_
    )

    # 组装结果行
    result_row <- data.frame(
      a = current_a, b = current_b,
      alpha_hat = alpha_hat,
      variance = suppressWarnings(as.numeric(var_results$variance[1])),
      std_error = se,
      p_value = p_value,
      bic = bic_val
    )

    results_list[[i]] <- result_row
  }

  # 合并
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df
}

# 调用多起点函数

calculate_bic_weighted_estimates <- function(results_df) {
  stopifnot(is.data.frame(results_df))
  # 标记可参与加权的行：三项都非 NA（也可顺带要求是有限数）
  valid_mask <- with(
    results_df,
    !is.na(variance) &
      is.finite(variance) &
      !is.na(bic) &
      is.finite(bic) &
      !is.na(alpha_hat) &
      is.finite(alpha_hat)
  )

  # 先拷贝一份，新增 weight 列（默认 NA）
  results_with_weights <- results_df
  results_with_weights$weight <- NA_real_

  # 若没有任何可加权的行
  if (!any(valid_mask)) {
    warning("No rows have non-NA variance, bic, and alpha_hat. Returning NA for BMA.")
    return(list(
      alpha_bma = NA_real_,
      se_bma = NA_real_,
      results_with_weights = results_with_weights
    ))
  }

  # 仅在有效子集上计算 BIC 权重
  vr <- results_df[valid_mask, , drop = FALSE]
  min_bic <- min(vr$bic)
  delta_bic <- vr$bic - min_bic
  exp_term <- exp(-0.5 * delta_bic)

  denom <- sum(exp_term)
  if (!is.finite(denom) || denom <= 0) {
    warning("Numerical issue when computing BIC weights. Returning NA for BMA.")
    return(list(
      alpha_bma = NA_real_,
      se_bma = NA_real_,
      results_with_weights = results_with_weights
    ))
  }

  weights <- exp_term / denom
  # 把权重写回到完整数据框的对应行；其它行保持 NA
  results_with_weights$weight[valid_mask] <- weights

  # --- 用有效子集做 BMA 统计 ---
  alpha_bma <- sum(weights * vr$alpha_hat)

  # 模型内方差
  within_model_var <- sum(weights * vr$variance)

  # 模型间方差
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


cml_family_new <- function(
  a_values,
  b_values,
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1 / 3 * 1
) {
  term1 <- tune_cml_family_ab_forloop_robust(
    a_values,
    b_values,
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    n = 1000,
    p = p,
    q = q
  )
  term2 <- calculate_bic_weighted_estimates(term1)
  return(term2)
}


# MRcML_family_dp


MRcML_family_dp <- function(a_values,
                            b_values,
                            beta_hat_exp,
                            beta_hat_out,
                            beta_sigma_exp,
                            beta_sigma_out,
                            T = 100, # 重复次数
                            n = 1000,
                            p = 1,
                            q = 1 / 3) {

  n_snps <- nrow(beta_hat_exp)

  # 1. 拟合主模型 (加入 try 保护，防止原始数据就无法拟合)
  alpha_ma <- try(cml_family_new(
    a_values = a_values, # 建议使用传入参数而不是硬编码 1:n_snps
    b_values = b_values,
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_sigma_exp = beta_sigma_exp,
    beta_sigma_out = beta_sigma_out,
    n = n, p = p, q = q
  ), silent = TRUE)

  if (inherits(alpha_ma, "try-error")) {
    warning("主模型拟合失败，返回 NA")
    return(data.frame(alpha_bma = NA, se_bma = NA, p_value = NA))
  }

  main_alpha <- alpha_ma$alpha_bma
  main_se <- alpha_ma$se_bma

  # 2. 准备 Bootstrap 抽样
  # 优化点：使用矩阵抽样，比原版索引切片运行更快、逻辑更清晰
  idx_mat <- matrix(sample(1:n_snps, T * n_snps, replace = TRUE), nrow = T, ncol = n_snps)
  theta_list <- numeric(T)
  theta_se_list <- numeric(T)

  # 3. Bootstrap 循环
  for (i in 1:T) {
    sel_idx <- idx_mat[i,]

    # 核心优化：加入 try() 防止单词循环报错导致整个函数崩溃
    point_t <- try(cml_family_new(
      a_values = a_values,
      b_values = b_values,
      beta_hat_exp = beta_hat_exp[sel_idx, , drop = FALSE],
      beta_hat_out = beta_hat_out[sel_idx, , drop = FALSE],
      beta_sigma_exp = beta_sigma_exp[sel_idx],
      beta_sigma_out = beta_sigma_out[sel_idx],
      n = n, p = p, q = q
    ), silent = TRUE)

    if (!inherits(point_t, "try-error")) {
      theta_list[i] <- point_t$alpha_bma
      theta_se_list[i] <- point_t$se_bma
    } else {
      theta_list[i] <- NA
      theta_se_list[i] <- NA
    }
  }

  # 4. 统计量计算与 NA 处理
  good_mask <- is.finite(theta_list)
  n_success <- sum(good_mask)

  # 如果成功次数太少（< 2次），无法计算方差
  if (n_success < 2) {
    warning("Bootstrap 成功次数过少，无法计算 DP 统计量，返回 NA")
    return(data.frame(alpha_bma = NA, se_bma = NA, p_value = NA))
  }

  theta_var <- var(theta_list[good_mask])
  theta_sd <- sqrt(theta_var)
  mean_se <- mean(theta_se_list[good_mask])
  var_se <- var(theta_se_list[good_mask])

  # 5. Z 检验
  z_fenzi <- mean_se - theta_sd
  z_fenmu <- sqrt(var_se + 2 / (n_success - 1) * theta_var^2)

  # 避免分母为 0 的极端情况
  z_p_value <- if (z_fenmu > 1e-12) pnorm(abs(z_fenzi / z_fenmu), lower.tail = FALSE) else 1

  # 6. 决策树逻辑与最终判定
  final_se <- NA_real_

  if (main_se > theta_sd) {
    if (z_p_value < 0.05) {
      # 遵循你的原逻辑：此极端情况返回 NA
      return(data.frame(alpha_bma = NA, se_bma = NA, p_value = NA))
    } else {
      final_se <- main_se
    }
  } else {
    if (z_p_value < 0.05) {
      final_se <- theta_sd
    } else {
      final_se <- theta_sd
    }
  }

  # 统一在此处计算 P 值，彻底解决 SE 和 P 值不匹配的隐患
  final_p <- 2 * pnorm(abs(main_alpha / final_se), lower.tail = FALSE)

  # 7. 返回结果
  return(data.frame(
    alpha_bma = main_alpha,
    se_bma = final_se,
    p_value = final_p
  ))
}


############################################################
## 1. 小工具 / 复用部件
############################################################

# 为了稳健：优先调用 C++ 版本的方差计算，如果不可用就 fallback 到 R 版本
.get_alpha_var_results <- function(alpha_hat,
                                   beta_hat_out, beta_hat_exp,
                                   beta_sigma_out, beta_sigma_exp,
                                   S_f, S_t) {
  out <- tryCatch(
  {
    calculate_alpha_variance_internal_two_cpp(
      alpha_hat = alpha_hat,
      beta_hat_out = beta_hat_out,
      beta_hat_exp = beta_hat_exp,
      beta_sigma_out = beta_sigma_out,
      beta_sigma_exp = beta_sigma_exp,
      S_f = S_f,
      S_t = S_t
    )
  },
    error = function(e) {
      calculate_alpha_variance_internal_two(
        alpha_hat = alpha_hat,
        beta_hat_out = beta_hat_out,
        beta_hat_exp = beta_hat_exp,
        beta_sigma_out = beta_sigma_out,
        beta_sigma_exp = beta_sigma_exp,
        S_f = S_f,
        S_t = S_t
      )
    }
  )
  out
}

# 计算 (a,b) 网格下的所有估计，带多数投票控制 (controller)
# 这个是 tune_cml_family_ab_forloop_robust 的高效重写版
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
  inv_beta_sigma_exp_precomp = NULL, # 允许外部传入，避免重复求逆
  controller = cml_bic_weighted_controller_ag
) {
  # 只算一次 inv_beta_sigma_exp，后续循环一直复用
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

      # --- 1) 运行控制器，得到 (alpha, gamma, S_f, S_t)
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
          q = q, n = n,
          p = p
        ),
        error = function(e) NULL
      )

      # 如果控制器失败或没收敛，直接给 NA 行（含 S_f/S_t 的 NA list 占位）
      if (is.null(est_results) || !is.list(est_results)) {
        results_list[[idx]] <- data.frame(
          a = a_cur,
          b = b_cur,
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

      # --- 2) 方差 / SE 估计
      var_results <- tryCatch(
        .get_alpha_var_results(
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

      p_val <- if (is.finite(alpha_hat_val) &&
        is.finite(se_hat_val) &&
        se_hat_val > 0) {
        2 * pnorm(abs(alpha_hat_val / se_hat_val), lower.tail = FALSE)
      } else {
        NA_real_
      }

      # --- 3) BIC 计算
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

      # --- 4) 存结果（S_f/S_t 用 list-column 保留集合）
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


# 计算 BIC 加权 (BMA)：保持你的原函数逻辑，但它本来就不依赖 tidyverse
# 返回:
#   alpha_bma
#   se_bma
#   results_with_weights (原 df + weight 列)
calculate_bic_weighted_estimates_fast_1 <- function(results_df) {
  stopifnot(is.data.frame(results_df))

  # 兼容列名：S_f/S_t 或 s_f/s_t
  col_Sf <- if ("S_f" %in% names(results_df)) "S_f" else if ("s_f" %in% names(results_df)) "s_f" else NULL
  col_St <- if ("S_t" %in% names(results_df)) "S_t" else if ("s_t" %in% names(results_df)) "s_t" else NULL

  valid_mask <- with(
    results_df,
    !is.na(variance) &
      is.finite(variance) &
      !is.na(bic) &
      is.finite(bic) &
      !is.na(alpha_hat) &
      is.finite(alpha_hat)
  )

  # 拷贝并加 weight 列以便可视化
  results_with_weights <- results_df
  results_with_weights$weight <- NA_real_

  if (!any(valid_mask)) {
    warning("No rows have non-NA variance/bic/alpha_hat. Returning NA for BMA.")
    return(list(
      alpha_bma = NA_real_,
      se_bma = NA_real_,
      results_with_weights = results_with_weights,
      s_f_best = if (!is.null(col_Sf)) NA else NULL,
      s_t_best = if (!is.null(col_St)) NA else NULL,
      best_row_index = NA_integer_
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
      results_with_weights = results_with_weights,
      s_f_best = if (!is.null(col_Sf)) NA else NULL,
      s_t_best = if (!is.null(col_St)) NA else NULL,
      best_row_index = NA_integer_
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

  # —— 取 BIC 最小行（在原始 results_df 中的行号）——
  local_best_idx <- which.min(vr$bic) # 在子集 vr 里的位置
  global_indices <- which(valid_mask) # 映射回原表
  best_row_index <- global_indices[local_best_idx]

  # 提取最优行的 S_f / S_t（若存在该列）
  s_f_best <- if (!is.null(col_Sf)) results_df[[col_Sf]][[best_row_index]] else NULL
  s_t_best <- if (!is.null(col_St)) results_df[[col_St]][[best_row_index]] else NULL

  list(
    alpha_bma = alpha_bma,
    se_bma = se_bma,
    results_with_weights = results_with_weights,
    s_f_best = s_f_best,
    s_t_best = s_t_best,
    best_row_index = best_row_index
  )
}

calculate_bic_weighted_estimates_fast <- function(results_df) {
  stopifnot(is.data.frame(results_df))

  # 兼容列名：集合列可能是 list-column（每行一个向量），也可能是已经存好的计数
  col_Sf <- if ("S_f" %in% names(results_df)) "S_f" else if ("s_f" %in% names(results_df)) "s_f" else NULL
  col_St <- if ("S_t" %in% names(results_df)) "S_t" else if ("s_t" %in% names(results_df)) "s_t" else NULL

  # 把“集合元素”转成“基数”的小工具（对单行取基数）
  sf_count_at <- function(col, idx) {
    if (is.null(col)) {
      return(NA_integer_)
    }
    el <- col[[idx]]
    if (is.null(el)) {
      return(0L)
    }
    # 若本来就是已存好的非负整数计数
    if (length(el) == 1 &&
      is.numeric(el) &&
      is.finite(el) &&
      el >= 0 &&
      abs(el - round(el)) < 1e-8) {
      return(as.integer(round(el)))
    }
    # 常见：list-column 里存放的是一个原子向量（如 integer/numeric/character）
    if (is.atomic(el)) {
      return(sum(!is.na(el)))
    }
    # 退而求其次：是个 list，就取长度
    if (is.list(el)) {
      return(length(el))
    }
    NA_integer_
  }

  valid_mask <- with(
    results_df,
    !is.na(variance) &
      is.finite(variance) &
      !is.na(bic) &
      is.finite(bic) &
      !is.na(alpha_hat) &
      is.finite(alpha_hat)
  )

  results_with_weights <- results_df
  results_with_weights$weight <- NA_real_

  if (!any(valid_mask)) {
    warning("No rows have non-NA variance/bic/alpha_hat. Returning NA for BMA.")
    return(list(
      alpha_bma = NA_real_,
      se_bma = NA_real_,
      results_with_weights = results_with_weights,
      s_f_best = if (!is.null(col_Sf)) NA else NULL,
      s_t_best = if (!is.null(col_St)) NA else NULL,
      s_f_best_count = if (!is.null(col_Sf)) NA_integer_ else NULL,
      s_t_best_count = if (!is.null(col_St)) NA_integer_ else NULL,
      best_row_index = NA_integer_
    ))
  }

  # —— 全局最优行（在 valid_mask 内 BIC 最小）——
  vr_valid <- results_df[valid_mask, , drop = FALSE]
  local_best_idx <- which.min(vr_valid$bic)
  global_indices <- which(valid_mask)
  best_row_index <- global_indices[local_best_idx]

  # —— 取最优行的 S_f / S_t（原值 + 基数）——
  s_f_best_raw <- if (!is.null(col_Sf)) results_df[[col_Sf]][[best_row_index]] else NULL
  s_t_best_raw <- if (!is.null(col_St)) results_df[[col_St]][[best_row_index]] else NULL
  s_f_best_count <- if (!is.null(col_Sf)) sf_count_at(results_df[[col_Sf]], best_row_index) else NA_integer_
  s_t_best_count <- if (!is.null(col_St)) sf_count_at(results_df[[col_St]], best_row_index) else NA_integer_

  # —— 若 |S_f|=0，则仅保留 |S_f|==0 的行参与 BMA；NA 视为未知，保守地剔除 ——
  effective_mask <- valid_mask
  if (!is.null(col_Sf) && isTRUE(!is.na(s_f_best_count) && s_f_best_count == 0L)) {
    # 为每一行计算 |S_f|
    sf_count_vec <- vapply(
      seq_len(nrow(results_df)),
      function(i) sf_count_at(results_df[[col_Sf]], i),
      integer(1)
    )
    keep_sf0 <- (!is.na(sf_count_vec) & sf_count_vec == 0L)
    effective_mask <- valid_mask & keep_sf0

    if (!any(effective_mask)) {
      warning("After filtering (|S_f|=0 -> keep only |S_f|==0), no rows remain. Returning NA.")
      return(list(
        alpha_bma = NA_real_,
        se_bma = NA_real_,
        results_with_weights = results_with_weights,
        s_f_best = s_f_best_raw,
        s_t_best = s_t_best_raw,
        s_f_best_count = s_f_best_count,
        s_t_best_count = s_t_best_count,
        best_row_index = best_row_index
      ))
    }
  }

  # —— 在有效集合上计算 BIC 权重 ——
  vr <- results_df[effective_mask, , drop = FALSE]

  min_bic <- min(vr$bic)
  delta_bic <- vr$bic - min_bic
  exp_term <- exp(-0.5 * delta_bic)
  denom <- sum(exp_term)

  if (!is.finite(denom) || denom <= 0) {
    warning("Numerical issue computing BIC weights.")
    return(list(
      alpha_bma = NA_real_,
      se_bma = NA_real_,
      results_with_weights = results_with_weights,
      s_f_best = s_f_best_raw,
      s_t_best = s_t_best_raw,
      s_f_best_count = s_f_best_count,
      s_t_best_count = s_t_best_count,
      best_row_index = best_row_index
    ))
  }

  weights <- exp_term / denom
  results_with_weights$weight[effective_mask] <- weights

  # BMA 计算
  alpha_bma <- sum(weights * vr$alpha_hat)
  within_model_var <- sum(weights * vr$variance)
  between_model_var <- sum(weights * (vr$alpha_hat - alpha_bma)^2)
  total_var <- within_model_var + between_model_var
  if (!is.finite(total_var) || total_var < 0) total_var <- NA_real_
  se_bma <- if (is.na(total_var)) NA_real_ else sqrt(total_var)

  list(
    alpha_bma = alpha_bma,
    se_bma = se_bma,
    results_with_weights = results_with_weights,
    s_f_best = s_f_best_raw, # 原始集合（或已存计数）
    s_t_best = s_t_best_raw,
    s_f_best_count = s_f_best_count, # 基数
    s_t_best_count = s_t_best_count, # 基数
    best_row_index = best_row_index
  )
}


############################################################
## 2. 普通版入口
############################################################
# 统一处理默认 (a,b) 网格，避免普通版 / DP 版各自维护一套逻辑。
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

# 普通版：单次全样本拟合 + BIC-BMA 聚合。
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

  # 1) 预先计算 inv_beta_sigma_exp，只做一次
  inv_beta_sigma_exp <- lapply(beta_sigma_exp, solve)

  # 2) 在 (a,b) 网格上运行多数投票+拟合
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

  # 3) BIC-BMA 加权，得到 alpha_bma / se_bma
  bic_weighted <- calculate_bic_weighted_estimates_fast(grid_fit)

  # 返回和你之前一样的结构
  bic_weighted
}

# 兼容旧接口名。
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


############################################################
## 3. Bootstrap / DP 核心
############################################################
# 这个核心函数做所有数值工作，但不画图。


############################################################
## 4. MRcML_family_dp() 轻量版 (无绘图依赖)
############################################################
# DP 版：bootstrap / DP 合成，仅返回主结果与可选中间量。
# MRcML_family_dp <- function(
#   a_values,
#   b_values,
#   beta_hat_exp,
#   beta_hat_out,
#   beta_sigma_exp,
#   beta_sigma_out,
#   T = 100,
#   n = 1000,
#   p = 1,
#   q = 1 / 3,
#   return_core = FALSE
# ) {
#   core_res <- MRcML_family_dp_core(
#     a_values = a_values,
#     b_values = b_values,
#     beta_hat_exp = beta_hat_exp,
#     beta_hat_out = beta_hat_out,
#     beta_sigma_exp = beta_sigma_exp,
#     beta_sigma_out = beta_sigma_out,
#     T = T,
#     n = n,
#     p = p,
#     q = q
#   )
#
#   out <- data.frame(
#     alpha_bma = core_res$alpha_final,
#     se_bma = core_res$se_final,
#     p_value = core_res$p_final
#   )
#
#   if (isTRUE(return_core)) {
#     attr(out, "core_res") <- core_res
#   }
#
#   out
# }

MRcML_family_dp_ver2 <- function(a_values,
                                 b_values,
                                 beta_hat_exp,
                                 beta_hat_out,
                                 beta_sigma_exp,
                                 beta_sigma_out,
                                 T = 100, # 重复次数
                                 n = 1000,
                                 p = 1,
                                 q = 1 / 3 * 1) {
  # ---- Palette & Theme (统一风格) ----
  # 注意： col_secondary 和 col_accent 是用于其他元素的（如直方图、渐变中点、粉点）。
  # 原代码中 BIC 颜色由 scale_fill_gradient2 独立控制。
  # 您的要求是浅蓝和深蓝顺序变化，我将在 scale_fill_distiller 中实现，不再使用这里定义的 col_secondary 作为中点。
  col_primary <- "#4061a0" # 主色（线/等高线/密度线）
  col_secondary <- "#77bbdd" # 浅蓝（直方图/中位色/渐变中点） -> 用于直方图
  col_accent <- "#eeb8c4" # 粉（另一端渐变） -> 保留给其他图层
  col_danger <- "#d0574e" # 历史变量名，BIC 图中不再使用红色
  col_shade <- "#b3cde3" # 阴影

  theme_pub <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  # ========= 1. 主估计（不重采样）+ bootstrap-like 重采样 =========
  n_snps <- nrow(beta_hat_exp)

  # 原始整体估计，用全部 SNP
  alpha_ma <- cml_family_new(
    a_values = 1:n_snps,
    b_values = 0:n_snps,
    beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    n = n,
    p = p,
    q = q
  )

  # bagging / DP 重采样
  sample_idx <- sample(1:n_snps, T * n_snps, replace = TRUE)
  theta_list <- numeric(T)
  theta_se_list <- numeric(T)

  for (i in 1:T) {
    index <- n_snps * (i - 1) + 1
    sel <- sample_idx[index:(index + n_snps - 1)]

    point_t <- cml_family_new(
      a_values = 1:n_snps,
      b_values = 0:n_snps,
      beta_hat_exp[sel,],
      beta_hat_out[sel,],
      beta_sigma_exp[sel],
      beta_sigma_out[sel],
      n = n,
      p = p,
      q = q
    )

    theta_list[i] <- point_t$alpha_bma
    theta_se_list[i] <- point_t$se_bma
  }

  # ========= 2. DP 方差整合 / 不确定性判断，得到最终 alpha, se, p =========
  theta_var <- var(theta_list)

  z_fenzi <- mean(theta_se_list) - sqrt(theta_var)
  z_fenmu <- sqrt(var(theta_se_list) + 2 / (T - 1) * theta_var^2)
  z <- abs(z_fenzi / z_fenmu)
  z_p_value <- pnorm(z, lower.tail = FALSE)

  alpha_final <- alpha_ma$alpha_bma
  se_final <- NA_real_
  p_final <- NA_real_

  if (alpha_ma$se_bma > sqrt(theta_var)) {
    if (z_p_value < 0.05) {
      se_final <- NA_real_
      p_final <- NA_real_
    } else {
      se_final <- alpha_ma$se_bma
      p_final <- 2 * pnorm(abs(alpha_ma$alpha_bma / alpha_ma$se_bma), lower.tail = FALSE)
    }
  } else {
    if (z_p_value < 0.05) {
      se_final <- sqrt(theta_var)
      p_final <- 2 * pnorm(abs(alpha_ma$alpha_bma / sqrt(theta_var)), lower.tail = FALSE)
    } else {
      se_final <- sqrt(theta_var)
      p_final <- 2 * pnorm(abs(alpha_ma$alpha_bma / alpha_ma$se_bma), lower.tail = FALSE)
    }
  }

  result <- tibble::tibble(
    alpha_bma = alpha_final,
    se_bma = se_final,
    p_value = p_final
  )

  # ========= 3. 画图 part 1：theta_list 的分布图 (plot_density) =========
  x <- theta_list
  x <- x[is.finite(x)]
  stopifnot(is.numeric(x), length(x) > 1)

  s <- stats::sd(x)
  qs <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)

  dens <- stats::density(x, adjust = 1.0, n = 2048)
  df_dens <- data.frame(x = dens$x, y = dens$y)
  df_hist <- data.frame(x = x)
  df_shade <- subset(df_dens, x >= qs[1] & x <= qs[2])

  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)

  plot_density <- ggplot(df_hist, aes(x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      fill = col_secondary,
      alpha = 0.35,
      color = NA
    ) +
    geom_area(
      data = df_shade,
      aes(x, y),
      alpha = 0.25,
      fill = col_shade
    ) +
    geom_line(
      data = df_dens,
      aes(x, y),
      linewidth = 1.1,
      color = col_primary
    ) +
    geom_rug(alpha = 0.35, sides = "b") +
    labs(
      title = "MRcML-family-MA-DP — Sampling Distribution",
      # 去掉了均值和中位数，上一次修复保留
      subtitle = sprintf(
        "T = %d, sd = %.3f; 95%% sample interval [%.3f, %.3f]",
        length(x), s, qs[1], qs[2]
      ),
      x = NULL,
      y = "Density"
    ) +
    theme_pub

  # ========= 4. 画图 part 2：BIC landscape (plot_bic) =========
  a_tbl <- alpha_ma$results_with_weights

  abic_all <- a_tbl %>%
    dplyr::select(a, b, bic) %>%
    dplyr::mutate(
      a = as.integer(a),
      b = as.integer(b)
    )

  abic_obs <- abic_all %>%
    dplyr::filter(!is.na(bic))

  if (nrow(abic_obs) == 0) {
    plot_bic <- ggplot() +
      annotate("text",
               x = 0, y = 0,
               label = "BIC 不可用 / 无非空BIC值",
               fontface = "bold",
               size = 5,
               colour = col_primary
      ) +
      theme_void() +
      labs(title = "BIC Landscape over (a, b)")
  } else {
    best <- abic_obs %>%
      dplyr::slice_min(bic, n = 1, with_ties = FALSE)

    a_vals_full <- sort(unique(abic_all$a))
    b_vals_full <- sort(unique(abic_all$b))

    full_grid <- tidyr::expand_grid(
      a = a_vals_full,
      b = b_vals_full
    ) %>%
      dplyr::left_join(abic_all, by = c("a", "b"))

    bic_range <- range(abic_obs$bic, na.rm = TRUE)
    legend_breaks <- if (diff(bic_range) > 0) bic_range else bic_range[1]

    choose_breaks_smart <- function(uvals, max_n = 10) {
      if (length(uvals) <= max_n) {
        uvals
      } else {
        step <- ceiling(length(uvals) / max_n)
        uvals[seq(1, length(uvals), by = step)]
      }
    }

    bx <- choose_breaks_smart(a_vals_full, max_n = 10)
    by <- choose_breaks_smart(b_vals_full, max_n = 10)

    axis_len <- max(length(a_vals_full), length(b_vals_full))
    txt_sz <- if (axis_len > 20) {
      6
    } else if (axis_len > 10) {
      7
    } else {
      10
    }

    plot_bic <- ggplot() +
      geom_tile(
        data = full_grid,
        aes(x = a, y = b),
        fill = "white",
        alpha = 0.4,
        color = "grey90"
      ) +
      geom_tile(
        data = abic_obs,
        aes(x = a, y = b, fill = bic),
        color = "grey90"
      ) +
      scale_fill_gradient(
        low = col_secondary,
        high = col_primary,
        breaks = legend_breaks,
        labels = function(x) sprintf("%.2f", x),
        name = "BIC'\n(lower is better)",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = 20,
          barheight = 1
        )
      ) +
      # 最优点也保持蓝色系，避免图中出现红色
      geom_point(
        data = best,
        aes(x = a, y = b),
        size = 3.2,
        shape = 21,
        fill = col_primary,
        color = "white",
        stroke = 0.5
      ) +
      scale_x_continuous(breaks = bx) +
      scale_y_continuous(breaks = by) +
      coord_fixed() +
      labs(
        title = "BIC Landscape over (a, b)",
        x = "a",
        y = "b"
      ) +
      theme_pub +
      theme(
        legend.position = "bottom",         # 【将图例移至底部】，底部布局能彻底解决字挤在一起的问题。
        # 增加一点图例和图表的间距，以防图例过长被截断
        legend.box.margin = margin(10, 10, 10, 10),
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = txt_sz
        ),
        axis.text.y = element_text(size = txt_sz)
      )
  }

  # ========= 5. 返回 =========
  return(list(
    results = result,
    theta_list = theta_list,
    plot_density = plot_density,
    plot_bic = plot_bic
  ))
}


# 依赖：ggplot2（仅在需要画图时）
# install.packages("ggplot2")


.normalize_indices <- function(idx, n_rows) {
  if (is.null(idx)) {
    return(integer(0))
  }
  # 若是 list 只取第一个元素（常见每次只需一个最优集合）
  if (is.list(idx)) idx <- idx[[1]]
  # 若是字符形如 "1,3,5" -> 转为数字
  if (is.character(idx) &&
    length(idx) == 1L &&
    grepl(",", idx)) {
    idx <- as.numeric(strsplit(idx, ",")[[1]])
  }
  # 逻辑向量 -> which
  if (is.logical(idx)) idx <- which(idx)
  # 数值向量清洗
  idx <- as.integer(unique(idx))
  idx <- idx[is.finite(idx) &
               !is.na(idx) &
               idx >= 1L &
               idx <= n_rows]
  idx
}
