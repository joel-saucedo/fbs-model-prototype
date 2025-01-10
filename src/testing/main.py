import fractional_bs_model
import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from pathlib import Path

# Calibration function for Python-SWIG bridge
def calibration_objective(alpha, options_data, option_type, S, sigma, r, dx, N, L, S0):
    errors = []
    for _, row in options_data.iterrows():
        K = row["strike"]
        T = row["days_to_maturity"] / 365
        market_price = row["mid_price"]

        # Call the C++ model
        try:
            model_price = fractional_bs_model.fractional_black_scholes_alpha(
                S, alpha, sigma, r, dx, N, L, K, T, option_type
            )
        except RuntimeError:
            continue

        # Interpolate model price at current spot price S0
        model_price_at_S0 = np.interp(S0, S, model_price)

        # Normalize model price by the underlying price
        normalized_model_price = model_price_at_S0 / S0

        # Calculate squared error
        error = (normalized_model_price - market_price) ** 2
        errors.append(error)

    if len(errors) == 0:
        return np.inf  # Return infinity if no valid options
    return np.sum(errors)


if __name__ == "__main__":
    # Load options data
    data_path = Path("../../../data/raw/btc_options_data.json")
    options_data = pd.read_json(data_path)

    # Model Parameters
    sigma = 0.6  # Volatility estimate
    r = 0.02     # Risk-free rate
    S_max = options_data["underlying_price"].max()  # Maximum Bitcoin price in USD
    S0 = options_data["underlying_price"].iloc[0]  # Current Bitcoin price

    # Discretization steps
    N = 200  # Number of spatial points
    L = 100  # Number of time steps
    dx = S_max / (N - 1)
    S = np.linspace(0, S_max, N)

    # Split options data into calls and puts
    call_options = options_data[options_data["option_type"] == "C"].copy()
    put_options = options_data[options_data["option_type"] == "P"].copy()

    # Calibration for calls
    result_call = minimize_scalar(
        calibration_objective,
        bounds=(0.1, 0.99),
        args=(call_options, 'C', S, sigma, r, dx, N, L, S0),
        method='bounded',
        options={'xatol': 1e-3}
    )
    alpha_call = result_call.x
    print(f"Calibrated alpha for calls: {alpha_call}")

    # Calibration for puts
    result_put = minimize_scalar(
        calibration_objective,
        bounds=(0.1, 0.99),
        args=(put_options, 'P', S, sigma, r, dx, N, L, S0),
        method='bounded',
        options={'xatol': 1e-3}
    )
    alpha_put = result_put.x
    print(f"Calibrated alpha for puts: {alpha_put}")
