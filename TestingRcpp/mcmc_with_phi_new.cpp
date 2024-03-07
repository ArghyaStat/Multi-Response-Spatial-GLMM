// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <numeric>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// Matern covariance function
double matern(double dist, double phi, double nu) {
    double out;
    if (dist == 0.0) {
        out = 1.0; // Covariance at the same point is 1
    }
    else {
        double kappa = sqrt(2 * nu) / phi;
        double part1 = pow(2, 1 - nu) / tgamma(nu);
        double part2 = kappa * dist;

        // Print out values for debugging
        Rcpp::Rcout << "part2: " << part2 << ", nu: " << nu << std::endl;

        out = part1 * pow(part2, nu) * R::bessel_k(dist, nu, 1);
    }

    Rcpp::Rcout << "Matern covariance computed: " << out << std::endl;
    return out;
}

// Compute covariance matrix K
arma::mat cov_mat(arma::mat locations, double phi, double nu) {
    int n = locations.n_rows;
    arma::mat K(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dist = fabs(locations(i, 0) - locations(j, 0));
            K(i, j) = matern(dist, phi, nu);
        }
    }

    Rcpp::Rcout << "Covariance matrix K computed." << K << std::endl;

    return K;
}

// Log-likelihood function
double log_likelihood(double phi, double nu, arma::vec Y, arma::mat locations) {
    // Compute covariance matrix K
    arma::mat K = cov_mat(locations, phi, nu);


    // Compute SVD decomposition of K
    arma::vec s;
    arma::mat U, V;
    arma::svd(U, s, V, K);

    // Calculate the inverse of D
    arma::mat S_inv = arma::diagmat(1 / s);

    

    //// Compute log-likelihood
    //arma::vec Y_trans = U.t() * Y;
    //arma::vec inv_s = 1 / s;
    //double loglik = -0.5 * (arma::sum(arma::log(s)) + arma::as_scalar(Y_trans.t() * (U * arma::diagmat(inv_s) * U.t() * Y_trans)));

    // Compute log-likelihood
    // arma::vec Y_trans = U.t() * Y;
  
    double loglik = -0.5 * (arma::sum(arma::log(s)) + arma::as_scalar(Y.t() * (V.t() * S_inv * U * Y)));

    //// Compute log-likelihood
    //double loglik = -0.5 * (arma::sum(arma::log(s)) + arma::as_scalar(Y.t() * (V * arma::diagmat(1 / s) * U.t() * Y)));

    Rcpp::Rcout << "Log likelihood computed." << loglik << std::endl;

    return loglik;
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

// Log of target posterior function
double log_target_posterior(double phi, double nu, arma::vec Y, arma::mat locations) {
    double log_likelihood_val = log_likelihood(phi, nu, Y, locations);

    double out = log_likelihood_val + phi_prior(phi);
    Rcpp::Rcout << "Log target posterior computed." << out << std::endl;
    return out;
}

// Metropolis-Hastings algorithm to sample phi
// [[Rcpp::export]]
arma::vec metropolis_hastings(double init_phi,
    double nu,
    int niters,
    double h,
    arma::vec Y,
    arma::mat locations) {
    arma::vec phi_chain(niters);
    phi_chain[0] = init_phi;

    arma::vec errors = sqrt(h) * arma::randn(niters);

    int accept_count = 0;

    for (int iter = 1; iter < niters; ++iter) {
        // Generate proposal
        double proposal_phi = phi_chain[iter - 1] + errors[iter];

        double log_acc = log_target_posterior(proposal_phi, nu, Y, locations) -
            log_target_posterior(phi_chain[iter - 1], nu, Y, locations);

        if (log(R::runif(0, 1)) < log_acc) {
            phi_chain[iter] = proposal_phi;
            accept_count++;
        }
        else {
            phi_chain[iter] = phi_chain[iter - 1];
        }

        // Progress update
        if (iter % (niters / 10) == 0) {
            Rcpp::Rcout << "Progress: " << (iter / (niters / 10)) * 10 << "%" << std::endl;
        }
    }

    Rcpp::Rcout << "Acceptance rate: " << static_cast<double>(accept_count) / niters << std::endl;

    return phi_chain;
}
