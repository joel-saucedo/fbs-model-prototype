#ifndef FRACTIONAL_BS_MODEL_H
#define FRACTIONAL_BS_MODEL_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <string>

// Function prototypes
std::vector<double> initial_condition_call(const std::vector<double>& S, double K);
std::vector<double> initial_condition_put(const std::vector<double>& S, double K);
std::vector<double> beta_coefficients(double alpha, int n);
void construct_tridiagonal_coeffs(int N, double r, double sigma, double dx, double dt,
                                  std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
std::vector<double> solve_tridiagonal(const std::vector<double>& a, const std::vector<double>& b,
                                      const std::vector<double>& c, const std::vector<double>& rhs);
                                      
inline std::vector<double> fractional_black_scholes_alpha(
    const std::vector<double>& S, double alpha, double sigma, double r,
    double dx, int N, int L, double K, double T, const std::string& option_type);

inline double calibration_objective(
    double alpha,
    const std::vector<std::vector<double>>& options_data,
    const std::string& option_type,
    const std::vector<double>& S,
    double sigma, double r, double dx, int N, int L, double S0);

#endif
