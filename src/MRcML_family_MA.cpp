#include <RcppArmadillo.h>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

// 优化的辅助函数：预计算逆矩阵，避免重复计算
inline arma::mat inv_M_j_cpp_fast(const arma::mat& M_j, const arma::mat& beta_sigma,
                                   std::vector<arma::mat>& cache, int cache_key) {
    // 检查缓存
    if (!cache[cache_key].is_empty()) {
        return cache[cache_key];
    }

    // 找到非零对角元素的索引
    arma::uvec indices_to_keep = arma::find(M_j.diag());

    if (indices_to_keep.is_empty()) {
        cache[cache_key] = arma::zeros<arma::mat>(3, 3);
        return cache[cache_key];
    }

    // 提取子矩阵并求逆
    arma::mat sub_mat = beta_sigma.submat(indices_to_keep, indices_to_keep);
    arma::mat inv_sub_mat;

    bool success = arma::inv(inv_sub_mat, sub_mat);
    if (!success) {
        inv_sub_mat = arma::pinv(sub_mat);
    }

    arma::mat inv_results = arma::zeros<arma::mat>(3, 3);
    inv_results.submat(indices_to_keep, indices_to_keep) = inv_sub_mat;

    cache[cache_key] = inv_results;
    return inv_results;
}

// [[Rcpp::export]]
Rcpp::List run_cml_family_internal_rcpp(const arma::mat& beta_hat_exp,
                                        const arma::mat& beta_hat_out,
                                        const Rcpp::List& inv_beta_sigma_exp,
                                        const Rcpp::List& beta_sigma_out,
                                        int a,
                                        int b,
                                        double alpha_init,
                                        int max_iter,
                                        double tol) {
    int J = beta_hat_exp.n_rows;

    // 初始化参数
    double alpha_k = alpha_init;
    arma::mat gamma_k = beta_hat_exp;

    // 转换 R lists 到 C++ vectors
    std::vector<arma::mat> inv_beta_sigma_exp_mats(J);
    std::vector<arma::mat> beta_sigma_out_mats(J);

    // 预计算协方差矩阵的子块和逆矩阵
    std::vector<arma::mat> inv_sigma_beta_fm(J);
    std::vector<double> sigma_beta_o_sq(J);

    for(int j = 0; j < J; ++j) {
        inv_beta_sigma_exp_mats[j] = Rcpp::as<arma::mat>(inv_beta_sigma_exp[j]);
        beta_sigma_out_mats[j] = Rcpp::as<arma::mat>(beta_sigma_out[j]);

        // 预计算常用的值
        sigma_beta_o_sq[j] = beta_sigma_out_mats[j](0, 0);
        arma::mat sigma_beta_fm = beta_sigma_out_mats[j].submat(1, 1, 2, 2);
        inv_sigma_beta_fm[j] = arma::inv_sympd(sigma_beta_fm);
    }

    // 预分配内存
    arma::vec f_values(J);
    arma::vec t_values(J);
    arma::mat gamma_k_plus_1(J, 3);

    // 缓存 M_j 对应的逆矩阵 (8种可能的组合: 000, 001, 010, 011, 100, 101, 110, 111)
    std::vector<std::vector<arma::mat>> inv_cache(J);
    for(int j = 0; j < J; ++j) {
        inv_cache[j].resize(8);
    }

    arma::uvec S_f_idx, S_t_idx;
    std::vector<bool> S_f_set(J, false), S_t_set(J, false);

    for (int k = 0; k < max_iter; ++k) {
        // --- Step 1: 更新 S_f 和 S_t ---
        double alpha_k_sq = alpha_k * alpha_k;

        // 计算 f_values 和 t_values
        for (int j = 0; j < J; ++j) {
            // 计算 f_values
            double diff_o = beta_hat_out(j, 0) - alpha_k * gamma_k(j, 0);
            f_values(j) = (diff_o * diff_o) / sigma_beta_o_sq[j];

            // 计算 t_values - 使用预计算的逆矩阵
            arma::vec d_fm(2);
            d_fm(0) = beta_hat_out(j, 1) - alpha_k * gamma_k(j, 1);
            d_fm(1) = beta_hat_out(j, 2) - alpha_k * gamma_k(j, 2);

            t_values(j) = arma::as_scalar(d_fm.t() * inv_sigma_beta_fm[j] * d_fm);
        }

        // 更新索引集合
        arma::uvec sf_indices = arma::sort_index(f_values);
        arma::uvec st_indices = arma::sort_index(t_values);

        // 重置集合
        std::fill(S_f_set.begin(), S_f_set.end(), false);
        std::fill(S_t_set.begin(), S_t_set.end(), false);

        if (a > 0) {
            S_f_idx = sf_indices.head(a);
            for(arma::uword i = 0; i < S_f_idx.n_elem; ++i) {
                S_f_set[S_f_idx(i)] = true;
            }
        }
        if (b > 0) {
            S_t_idx = st_indices.head(b);
            for(arma::uword i = 0; i < S_t_idx.n_elem; ++i) {
                S_t_set[S_t_idx(i)] = true;
            }
        }

        // --- Step 2: 更新 gamma ---
        #pragma omp parallel for if(J > 100)
        for (int j = 0; j < J; ++j) {
            // 构建 M_j 并计算 cache_key
            int in_f = S_f_set[j] ? 1 : 0;
            int in_t = S_t_set[j] ? 1 : 0;
            int cache_key = (in_f << 2) | (in_t << 1) | in_t;

            arma::mat M_j = arma::zeros<arma::mat>(3, 3);
            M_j(0, 0) = in_f;
            M_j(1, 1) = in_t;
            M_j(2, 2) = in_t;

            arma::mat inv_beta_sigma_out_j = inv_M_j_cpp_fast(M_j, beta_sigma_out_mats[j],
                                                               inv_cache[j], cache_key);

            // 使用 Woodbury 矩阵恒等式优化求逆（如果适用）
            arma::mat term1_inv = inv_beta_sigma_exp_mats[j] + alpha_k_sq * inv_beta_sigma_out_j;

            arma::vec term2 = (inv_beta_sigma_exp_mats[j] * beta_hat_exp.row(j).t()) +
                              (alpha_k * (inv_beta_sigma_out_j * beta_hat_out.row(j).t()));

            // 使用 solve 而不是显式求逆
            gamma_k_plus_1.row(j) = arma::solve(term1_inv, term2).t();
        }

        // --- Step 3: 更新 alpha ---
        double numerator = 0.0;
        double denominator = 0.0;

        for (int j = 0; j < J; ++j) {
            int in_f = S_f_set[j] ? 1 : 0;
            int in_t = S_t_set[j] ? 1 : 0;
            int cache_key = (in_f << 2) | (in_t << 1) | in_t;

            arma::mat M_j = arma::zeros<arma::mat>(3, 3);
            M_j(0, 0) = in_f;
            M_j(1, 1) = in_t;
            M_j(2, 2) = in_t;

            arma::rowvec gamma_j = gamma_k_plus_1.row(j);
            arma::mat term_middle = inv_cache[j][cache_key];

            // 优化矩阵向量乘法
            arma::vec temp = term_middle * beta_hat_out.row(j).t();
            numerator += arma::dot(gamma_j, temp);

            temp = term_middle * gamma_j.t();
            denominator += arma::dot(gamma_j, temp);
        }

        if (denominator == 0 || std::abs(denominator) < 1e-15) break;

        double alpha_k_plus_1 = numerator / denominator;

        // 检查收敛
        if (std::abs(alpha_k_plus_1 - alpha_k) < tol) {
            alpha_k = alpha_k_plus_1;
            gamma_k = gamma_k_plus_1;
            break;
        }

        alpha_k = alpha_k_plus_1;
        gamma_k = gamma_k_plus_1;
    }

    // 转换为 R 索引（1-based）
    arma::vec S_f_r, S_t_r;
    if (a > 0 && !S_f_idx.empty()) {
        S_f_r = arma::conv_to<arma::vec>::from(S_f_idx) + 1;
    }
    if (b > 0 && !S_t_idx.empty()) {
        S_t_r = arma::conv_to<arma::vec>::from(S_t_idx) + 1;
    }

    return Rcpp::List::create(
        Rcpp::Named("alpha") = alpha_k,
        Rcpp::Named("gamma") = gamma_k,
        Rcpp::Named("S_f") = S_f_r,
        Rcpp::Named("S_t") = S_t_r
    );
}