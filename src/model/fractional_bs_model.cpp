#include "fractional_bs_model.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>

// Initial condition for Call options
std::vector<double> initial_condition_call(const std::vector<double>& S, double K) {
    std::cout << "Calculating initial condition for Call options...\n";
    std::vector<double> U(S.size());
    for (size_t i = 0; i < S.size(); ++i) {
        U[i] = std::max(S[i] - K, 0.0);
    }
    return U;
}

// Initial condition for Put options
std::vector<double> initial_condition_put(const std::vector<double>& S, double K) {
    std::cout << "Calculating initial condition for Put options...\n";
    std::vector<double> U(S.size());
    for (size_t i = 0; i < S.size(); ++i) {
        U[i] = std::max(K - S[i], 0.0);
    }
    return U;
}

// Compute beta coefficients for fractional derivative
std::vector<double> beta_coefficients(double alpha, int n) {
    std::cout << "Calculating beta coefficients for alpha = " << alpha << ", n = " << n << "\n";
    std::vector<double> beta(n);
    for (int j = 1; j <= n; ++j) {
        if (j == 1) {
            beta[j - 1] = std::pow(j, 1 - alpha);
        } else {
            beta[j - 1] = std::pow(j, 1 - alpha) - std::pow(j - 1, 1 - alpha);
        }
    }
    return beta;
}

// Construct tridiagonal matrix coefficients
void construct_tridiagonal_coeffs(int N, double r, double sigma, double dx, double dt,
                                  std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
    std::cout << "Constructing tridiagonal matrix coefficients...\n";
    for (int i = 0; i < N; ++i) {
        double i_sq = static_cast<double>(i * i);
        if (i > 0) {
            a[i - 1] = -0.5 * sigma * sigma * i_sq / (dx * dx) - r * i / (2 * dx);
        }
        b[i] = 1.0 / dt + sigma * sigma * i_sq / (dx * dx) + r;
        if (i < N - 1) {
            c[i] = -0.5 * sigma * sigma * i_sq / (dx * dx) + r * i / (2 * dx);
        }
    }
    std::cout << "Matrix coefficients constructed.\n";
}

// Solve tridiagonal system using Thomas algorithm
std::vector<double> solve_tridiagonal(const std::vector<double>& a, const std::vector<double>& b,
                                      const std::vector<double>& c, const std::vector<double>& rhs) {
    std::cout << "Solving tridiagonal system...\n";
    int N = rhs.size();
    std::vector<double> x(N, 0.0);
    std::vector<double> c_prime(N - 1, 0.0);
    std::vector<double> d_prime(N, 0.0);

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = rhs[0] / b[0];
    for (int i = 1; i < N; ++i) {
        double denom = b[i] - a[i - 1] * c_prime[i - 1];
        if (i < N - 1) {
            c_prime[i] = c[i] / denom;
        }
        d_prime[i] = (rhs[i] - a[i - 1] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    x[N - 1] = d_prime[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    std::cout << "Tridiagonal system solved.\n";
    return x;
}

// Solve the fractional Black-Scholes PDE
std::vector<double> fractional_black_scholes_alpha(
    const std::vector<double>& S, double alpha, double sigma, double r,
    double dx, int N, int L, double K, double T, const std::string& option_type) {
        
    std::cout << "Solving fractional Black-Scholes PDE...\n";
    double dt = T / (L - 1);

    std::vector<std::vector<double>> U(L, std::vector<double>(N, 0.0));
    if (option_type == "C") {
        U[0] = initial_condition_call(S, K);
    } else if (option_type == "P") {
        U[0] = initial_condition_put(S, K);
    }

    std::vector<double> beta = beta_coefficients(alpha, L);
    double gamma_alpha = 1.0 / std::tgamma(2 - alpha);

    std::vector<double> a(N - 1), b(N), c(N - 1);
    construct_tridiagonal_coeffs(N, r, sigma, dx, dt, a, b, c);

    for (int n = 1; n < L; ++n) {
        std::cout << "Time step: " << n << " / " << L << "\n";
        std::vector<double> rhs(N);
        for (int i = 0; i < N; ++i) {
            rhs[i] = U[n - 1][i] / dt;
        }
        U[n] = solve_tridiagonal(a, b, c, rhs);

        for (int j = 1; j <= n; ++j) {
            for (int i = 0; i < N; ++i) {
                U[n][i] += gamma_alpha * std::pow(dt, 1 - alpha) * beta[j - 1] * (U[n - j][i] - U[n - j - 1][i]);
            }
        }
    }

    std::cout << "PDE solved successfully.\n";
    return U.back();
}

// Calibration function for alpha optimization
double calibration_objective(
    double alpha,
    const std::vector<std::vector<double>>& options_data,
    const std::string& option_type,
    const std::vector<double>& S,
    double sigma, double r, double dx, int N, int L, double S0) {

    std::cout << "Calibrating alpha: " << alpha << "\n";
    std::vector<double> errors;

    for (const auto& row : options_data) {
        double K = row[0];
        double T = row[1] / 365.0;
        double market_price = row[2];

        std::vector<double> model_price = fractional_black_scholes_alpha(S, alpha, sigma, r, dx, N, L, K, T, option_type);

        auto lower = std::lower_bound(S.begin(), S.end(), S0);
        if (lower == S.end() || lower == S.begin()) {
            std::cerr << "S0 out of bounds: " << S0 << "\n";
            continue;
        }
        size_t idx = lower - S.begin();
        double model_price_at_S0 = model_price[idx];

        double normalized_model_price = model_price_at_S0 / S0;
        double error = std::pow(normalized_model_price - market_price, 2);
        errors.push_back(error);
    }

    if (errors.empty()) {
        std::cerr << "No valid options for calibration.\n";
        return std::numeric_limits<double>::infinity();
    }

    return std::accumulate(errors.begin(), errors.end(), 0.0);
}
