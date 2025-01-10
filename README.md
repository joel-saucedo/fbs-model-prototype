# Fractional Black-Scholes PDE Solver Prototype

This repository implements a **Fractional Black-Scholes PDE Solver** in **C++** with Python bindings using **SWIG**. The solver is designed for modeling advanced financial stochastic processes, including **cryptocurrency options**, by leveraging fractional calculus. A calibration mechanism optimizes the fractional order of the time derivative (alpha), enabling enhanced modeling of market anomalies like memory effects and heavy tails.

---

## Key Features

- **Fractional Black-Scholes PDE Solver**:
  - Implements the fractional-order Black-Scholes equation using Caputo fractional derivatives.
  - Numerically solves the PDE using finite difference methods with a **tridiagonal matrix system**.

- **Alpha-Order Calibration**:
  - Optimizes the alpha parameter based on market data for both call and put options.
  - Employs a calibration algorithm that minimizes error between model predictions and market prices.

- **Cryptocurrency Options Modeling**:
  - Focused on **Bitcoin options** data as a use case.
  - Extensible to other financial assets and derivative instruments.

- **Python-C++ Interoperability**:
  - Provides a Python interface for ease of use.
  - High-performance numerical computations are handled in C++ for efficiency.

---

## Algorithms and Methodology

### Fractional Black-Scholes PDE Solver

The fractional Black-Scholes PDE is solved numerically using:

1. **Initial Conditions**:
   - Call options: \( U(S, 0) = \max(S - K, 0) \).
   - Put options: \( U(S, 0) = \max(K - S, 0) \).

2. **Finite Difference Discretization**:
   - Discretizes the PDE over time and space using fixed step sizes.
   - Constructs a tridiagonal matrix to approximate spatial derivatives.

3. **Fractional Time Derivative**:
   - Implements the Caputo fractional derivative with a **beta coefficient series**:
     \[
     \beta_j = j^{1-\alpha} - (j-1)^{1-\alpha}
     \]
   - Uses these coefficients to calculate time-fractional updates to the solution grid.

4. **Numerical Solution**:
   - Solves the tridiagonal matrix system using the **Thomas algorithm**.
   - Incorporates memory effects by iterating through previous time steps with weighted contributions.

### Alpha Calibration Algorithm

The alpha parameter, governing the fractional derivative order, is optimized using:

- **Objective Function**:
  - Minimizes the squared error between market option prices and model predictions:
    \[
    \text{Error} = \sum_{i=1}^N \left( \frac{\text{ModelPrice}(S_0)}{S_0} - \text{MarketPrice} \right)^2
    \]

- **Market Data**:
  - Utilizes strike prices, time to maturity, and option mid-prices.

- **Optimization**:
  - Uses `scipy.optimize.minimize_scalar` to find the optimal alpha value.

This approach is based on research literature, with additional resources provided in the `docs/` directory.

---

## Installation

### Prerequisites

1. **Python**:
   - Version 3.10+.
   - Install dependencies using:
     ```bash
     python setup.py
     ```

2. **C++**:
   - Compiler supporting C++11 or higher (e.g., GCC, Clang).
   - SWIG (to generate Python bindings).

### Instructions

- **Build the C++ Library: Ensure SWIG is installed:**:
   ```bash
   swig -python -c++ -o src/cpp/fractional_bs_model_wrap.cxx src/cpp/fractional_bs_model.i
   g++ -shared -std=c++11 -fPIC src/cpp/fractional_bs_model_wrap.cxx src/cpp/fractional_bs_model.cpp \ 
   -o src/python/_fractional_bs_model.so \ 
   -I/usr/include/python3.10 \
   -I$(python -c "import numpy; print(numpy.get_include())")
   ```

### Refrences

- This prototype is based on fractional Black-Scholes literature (see docs/ for more information).
- Key concepts include Caputo fractional derivatives and finite difference methods for solving time-fractional PDEs.

### License

- This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

This project was made possible thanks to the resources from the following:

- **SWIG**:
  - The seamless integration between C++ and Python was achieved using [SWIG](http://www.swig.org/), a powerful tool for generating language bindings.

- **Cryptocurrency Data**:
  - The project utilized real-world Bitcoin options data from Deberit API to model and calibrate the fractional Black-Scholes PDE.





