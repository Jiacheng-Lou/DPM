// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <map>
#include <cmath>
#include <limits>

using namespace Rcpp;

static double safe_variance(double d0, double d1, double d2, double mu)
{
    double v = d0 + d1 * mu + d2 * mu * mu;
    if (!R_FINITE(v) || v <= 1e-14)
        v = 1e-14;
    return v;
}

static double log_sum_exp(const NumericVector &log_probs)
{
    double max_lp = R_NegInf;

    for (int i = 0; i < log_probs.size(); ++i)
    {
        if (R_FINITE(log_probs[i]) && log_probs[i] > max_lp)
        {
            max_lp = log_probs[i];
        }
    }

    if (!R_FINITE(max_lp))
        return R_NegInf;

    double total = 0.0;
    for (int i = 0; i < log_probs.size(); ++i)
    {
        if (R_FINITE(log_probs[i]))
        {
            total += std::exp(log_probs[i] - max_lp);
        }
    }

    if (total <= 0.0 || !R_FINITE(total))
        return R_NegInf;

    return max_lp + std::log(total);
}

static int sample_from_log_probs(const NumericVector &log_probs)
{
    int H = log_probs.size();
    double lse = log_sum_exp(log_probs);

    if (!R_FINITE(lse))
    {
        // Last-resort fallback: sample uniformly from all finite entries.
        std::vector<int> valid;
        for (int h = 0; h < H; ++h)
        {
            if (R_FINITE(log_probs[h]))
                valid.push_back(h);
        }

        if (valid.empty())
        {
            return static_cast<int>(std::floor(R::runif(0.0, 1.0) * H));
        }

        int idx = static_cast<int>(std::floor(R::runif(0.0, 1.0) * valid.size()));
        if (idx >= static_cast<int>(valid.size()))
            idx = valid.size() - 1;
        return valid[idx];
    }

    double u = R::runif(0.0, 1.0);
    double cum = 0.0;

    for (int h = 0; h < H; ++h)
    {
        double prob = R_FINITE(log_probs[h]) ? std::exp(log_probs[h] - lse) : 0.0;
        cum += prob;
        if (u <= cum)
            return h;
    }

    return H - 1;
}

// [[Rcpp::export]]
List update_clusters_cpp_hetero(
    NumericMatrix x_train,
    NumericVector y_train,
    NumericVector d0,
    NumericVector d1,
    NumericVector d2,
    IntegerVector z,
    NumericMatrix gamma,
    double alpha,
    int m,
    double prior_sd)
{
    const int N = x_train.nrow();
    const int p = x_train.ncol();
    const int K_start = gamma.nrow();

    if (N <= 0 || p <= 0)
    {
        stop("x_train must have positive numbers of rows and columns.");
    }
    if (y_train.size() != N || d0.size() != N || d1.size() != N || d2.size() != N)
    {
        stop("y_train, d0, d1 and d2 must have length nrow(x_train).");
    }
    if (z.size() != N)
    {
        stop("z must have length nrow(x_train).");
    }
    if (gamma.ncol() != p)
    {
        stop("ncol(gamma) must match ncol(x_train).");
    }
    if (K_start <= 0)
    {
        stop("gamma must contain at least one active cluster.");
    }
    if (alpha <= 0.0 || m <= 0 || prior_sd <= 0.0 || !R_FINITE(prior_sd))
    {
        stop("alpha, m and prior_sd must be positive.");
    }

    // Store gamma dynamically so selected auxiliary clusters can be appended.
    std::vector<std::vector<double>> gamma_vec(K_start, std::vector<double>(p));
    for (int k = 0; k < K_start; ++k)
    {
        for (int j = 0; j < p; ++j)
        {
            gamma_vec[k][j] = gamma(k, j);
        }
    }

    // Current cluster sizes, using 0-based internal labels.
    std::vector<int> cluster_sizes(K_start, 0);
    for (int i = 0; i < N; ++i)
    {
        if (z[i] < 1 || z[i] > K_start)
        {
            stop("z contains labels outside 1:nrow(gamma).");
        }
        cluster_sizes[z[i] - 1]++;
    }

    // Sequential Gibbs update over genes.
    for (int i = 0; i < N; ++i)
    {
        int old_c = z[i] - 1;
        if (old_c < 0 || old_c >= static_cast<int>(cluster_sizes.size()))
        {
            stop("Internal label error during cluster update.");
        }

        // Temporarily remove observation i.
        cluster_sizes[old_c]--;
        z[i] = 0;

        const int current_K = gamma_vec.size();
        const int total_proposals = current_K + m;

        // Draw m auxiliary gamma vectors from the proposal/base distribution.
        std::vector<std::vector<double>> aux_gamma(m, std::vector<double>(p));
        for (int aux = 0; aux < m; ++aux)
        {
            for (int j = 0; j < p; ++j)
            {
                aux_gamma[aux][j] = R::rnorm(0.0, prior_sd);
            }
        }

        NumericVector log_probs(total_proposals);

        // Existing clusters.
        for (int k = 0; k < current_K; ++k)
        {
            if (cluster_sizes[k] <= 0)
            {
                log_probs[k] = R_NegInf;
                continue;
            }

            double mu = 0.0;
            for (int j = 0; j < p; ++j)
            {
                mu += x_train(i, j) * gamma_vec[k][j];
            }

            double v_i = safe_variance(d0[i], d1[i], d2[i], mu);
            double log_lik = R::dnorm(y_train[i], mu, std::sqrt(v_i), true);

            // The denominator N - 1 + alpha is common to all choices and cancels.
            log_probs[k] = std::log(static_cast<double>(cluster_sizes[k])) + log_lik;
        }

        // Auxiliary clusters.
        for (int aux = 0; aux < m; ++aux)
        {
            double mu = 0.0;
            for (int j = 0; j < p; ++j)
            {
                mu += x_train(i, j) * aux_gamma[aux][j];
            }

            double v_i = safe_variance(d0[i], d1[i], d2[i], mu);
            double log_lik = R::dnorm(y_train[i], mu, std::sqrt(v_i), true);

            // The denominator N - 1 + alpha is common and cancels.
            log_probs[current_K + aux] = std::log(alpha / static_cast<double>(m)) + log_lik;
        }

        int chosen = sample_from_log_probs(log_probs);

        if (chosen < current_K)
        {
            z[i] = chosen + 1;
            cluster_sizes[chosen]++;
        }
        else
        {
            int aux_idx = chosen - current_K;

            // Keep exactly the auxiliary gamma used in the assignment probability.
            gamma_vec.push_back(aux_gamma[aux_idx]);
            cluster_sizes.push_back(1);

            z[i] = gamma_vec.size(); // 1-based new label
        }
    }

    // Compress active clusters to labels 1..K_new.
    std::vector<int> active_clusters;
    for (int k = 0; k < static_cast<int>(cluster_sizes.size()); ++k)
    {
        if (cluster_sizes[k] > 0)
        {
            active_clusters.push_back(k);
        }
    }

    const int K_new = active_clusters.size();
    NumericMatrix new_gamma(K_new, p);
    IntegerVector old_cluster_map(K_new);
    std::map<int, int> index_mapping;

    for (int new_k = 0; new_k < K_new; ++new_k)
    {
        int old_internal_idx = active_clusters[new_k];
        index_mapping[old_internal_idx] = new_k + 1;

        for (int j = 0; j < p; ++j)
        {
            new_gamma(new_k, j) = gamma_vec[old_internal_idx][j];
        }

        if (old_internal_idx < K_start)
        {
            old_cluster_map[new_k] = old_internal_idx + 1;
        }
        else
        {
            old_cluster_map[new_k] = 0;
        }
    }

    IntegerVector new_z(N);
    for (int i = 0; i < N; ++i)
    {
        int old_internal_idx = z[i] - 1;

        if (index_mapping.find(old_internal_idx) == index_mapping.end())
        {
            stop("Internal error: active label not found during relabelling.");
        }

        new_z[i] = index_mapping[old_internal_idx];
    }

    return List::create(
        Named("z") = new_z,
        Named("gamma") = new_gamma,
        Named("old_cluster_map") = old_cluster_map);
}