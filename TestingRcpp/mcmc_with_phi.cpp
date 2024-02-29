#include <Rcpp.h>
#include <cmath>
#include <numeric>

using namespace Rcpp;

// Matern covariance function
double matern(double dist, double phi, double nu) {
    double kappa = sqrt(2 * nu) / phi;
    double part1 = pow(2, 1 - nu) / tgamma(nu);
    double part2 = kappa * dist;
    return part1 * pow(part2, nu) * R::bessel_k(part2, nu, 1);
}

// Compute covariance matrix K
NumericMatrix compute_covariance_matrix(NumericVector locations, double phi, double nu) {
    int n = locations.size();
    NumericMatrix K(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dist = fabs(locations[i] - locations[j]);
            K(i, j) = matern(dist, phi, nu);
        }
    }
    return K;
}

// Cholesky decomposition function
NumericMatrix cholesky(NumericMatrix A) {
    int n = A.nrow();
    NumericMatrix L(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = A(i, j);
            for (int k = 0; k < j; ++k)
                sum -= L(i, k) * L(j, k);
            if (i == j) {
                if (sum <= 0) {
                    throw std::runtime_error("Cholesky decomposition failed: Matrix is not positive definite");
                }
                L(i, j) = sqrt(sum);
            }
            else {
                L(i, j) = sum / L(j, j);
            }
        }
    }

    return L;
}

// Cholesky decomposition of the inverse function
NumericMatrix cholesky_inverse(NumericMatrix L) {
    int n = L.nrow();
    NumericMatrix L_inv(n, n);

    for (int i = 0; i < n; ++i) {
        L_inv(i, i) = 1.0 / L(i, i);
        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = i; k < j; ++k)
                sum -= L(j, k) * L_inv(k, i);
            L_inv(j, i) = sum / L(j, j);
        }
    }

    return L_inv;
}

// Log of prior distribution for phi
double phi_prior(double phi) {
    if (phi > 0 && phi < 1) {
        return 0; // log(1) = 0 for phi in (0, 1)
    }
    else {
        return -INFINITY; // Return negative infinity for log(0)
    }
}

// Log-likelihood function
double log_likelihood(double phi, double nu, NumericVector Y, NumericVector locations) {
    // Check if phi is within (0, 1)
    if (phi <= 0 || phi >= 1) return -INFINITY;

    // int n = Y.size();

    // Compute covariance matrix K
    NumericMatrix K = compute_covariance_matrix(locations, phi, nu);

    // Compute Cholesky decomposition of covariance matrix
    NumericMatrix L = cholesky(K);

    // Compute Cholesky decomposition of inverse
    NumericMatrix L_inv = cholesky_inverse(L);

    // Compute log-likelihood
    NumericVector y_vec(Y);
    NumericVector y_trans = L_inv * y_vec;
    double loglik = -0.5 * (2 * sum(log(diag(L))) +
    sum(pow(y_trans, 2)));
    return loglik;
}

// Log of target posterior function
double log_target_posterior(double phi, double nu, NumericVector Y, NumericVector locations) {

    double log_likelihood_val = log_likelihood(phi, nu, Y, locations);
    return log_likelihood_val + phi_prior(phi);
}

// Metropolis-Hastings algorithm to sample phi
// [[Rcpp::export]]
NumericVector metropolis_hastings(double init_phi,
    double nu,
    int niters,
    double h,
    NumericVector Y,
    NumericVector locations) {
    NumericVector phi_chain(niters);
    phi_chain[0] = init_phi;

    NumericVector errors = sqrt(h) * rnorm(niters, 0, 1);

    int accept_count = 0;

    for (int iter = 1; iter < niters; ++iter) {
        // Generate proposal
        double proposal_phi = phi_chain[iter - 1] + errors[iter];


        double log_acc = log_target_posterior(proposal_phi, nu, Y, locations) -
            log_target_posterior(phi_chain[iter - 1], nu, Y, locations);

        if (log(runif(1)[0]) < log_acc) {
            phi_chain[iter] = proposal_phi;
            accept_count++;
        }
        else {
            phi_chain[iter] = phi_chain[iter - 1];
        }

        // Progress update
        if (iter % (niters / 10) == 0) {
            Rcout << "Progress: " << (iter / (niters / 10)) * 10 << "%" << std::endl;
        }
    }

    Rcout << "Acceptance rate: " << static_cast<double>(accept_count) / niters << std::endl;

    return phi_chain;
}
