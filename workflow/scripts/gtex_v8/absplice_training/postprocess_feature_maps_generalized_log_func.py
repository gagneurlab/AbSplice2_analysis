import pandas as pd
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from sklearn.isotonic import IsotonicRegression
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ------------------ Generalized Logistic Function ------------------
def generalized_logistic(x, A, K, B, Q, C, M):
    """
    Generalized logistic function (Richards curve):
    A: lower asymptote
    K: upper asymptote
    B: growth rate (negative for decreasing)
    Q: affects curve symmetry
    C: affects inflection point location
    M: curve's midpoint (inflection point)
    """
    return A + (K - A) / ((C + np.exp(-B * (x - M))) ** (1 / Q))


# ------------------ Weight Transformation Function ------------------
def transform_weights(weights, method):
    """Applies specified transformation to weights."""
    if method == "none":
        return None
    elif method == "raw":
        return weights
    elif method == "normalize":
        return weights / np.max(weights)
    elif method == "log":
        return np.log1p(weights)
    elif method == "sqrt":
        return np.sqrt(weights)
    elif method == "inverse":
        return 1 / (weights + 1)
    else:
        raise ValueError(f"Unknown weight transformation: {method}")


# ------------------ Fitting Function ------------------
def fit_generalized_logistic_with_weights(
    x, y, weights_trimmed, weight_method="none", monotonic_type='decreasing'
):
    """
    Fits a generalized logistic function with:
    Customizable weight transformation
    Monotonicity control
    """
    weights = transform_weights(weights_trimmed, weight_method)

    # Initial parameter guesses
    A_init, K_init = min(y), max(y)
    B_init = -5 if monotonic_type == 'decreasing' else 5
    Q_init, C_init, M_init = 1.0, 1.0, np.median(x)
    p0 = [A_init, K_init, B_init, Q_init, C_init, M_init]

    # Bounds for parameters
    lower_bounds = [-np.inf, -np.inf, -np.inf, 0.01, 0.01, min(x)]
    upper_bounds = [np.inf, np.inf, 0 if monotonic_type == 'decreasing' else np.inf, 10, 10, max(x)]

    # Residuals function
    def residuals(params, x_data, y_data, weights_data):
        A, K, B, Q, C, M = params
        y_pred = generalized_logistic(x_data, A, K, B, Q, C, M)
        return (y_data - y_pred) * (weights_data if weights_data is not None else 1)

    # Fit using curve_fit
    popt, _ = curve_fit(
        lambda x, A, K, B, Q, C, M: generalized_logistic(x, A, K, B, Q, C, M),
        x, y, p0=p0,
        bounds=(lower_bounds, upper_bounds),
        sigma=(1 / weights if weights is not None else None),
        absolute_sigma=True, maxfev=20000
    )

    return popt


# ------------------ Plotting Function ------------------
def plot_weight_transformation(ax, weight_method, x, y, weights):
    """Plots the fitted curve using a specific weight transformation."""
    popt = fit_generalized_logistic_with_weights(x, y, weights, weight_method=weight_method)
    x_fit = np.linspace(min(x), max(x), 500)
    y_fit = generalized_logistic(x_fit, *popt)

    ax.scatter(x, y, label='Data', color='blue', s=10)
    ax.plot(x_fit, y_fit, label='Fitted Curve', color='red', linewidth=2)
    ax.set_title(f"Weights: {weight_method}")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    ax.grid(True)
    
    
# Load the model from the pkl file
with open(snakemake.input['model'], "rb") as file:
    ebm_model = pickle.load(file)
    

# TODO: maybe trimming needs to be adapted (with y[1:-1] seems to miss out-of-range values. use y[2:])

# update the feautre functions
if 'loss_score' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('loss_score')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed = weights[1:-1]

    # weight_methods = ["none", "raw", "normalize", "log", "sqrt", "inverse"]
    weight_method = "sqrt"

    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    x_fit = np.linspace(min(x), max(x), 500)
    y_fit = generalized_logistic(x_fit, *popt)
    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    

if 'gain_score' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('gain_score')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed = weights[1:-1]

    # neglect high gain scores in fitting (not enough data to learn correct contributions)
    weights_trimmed[x_trimmed > 0.6] = 0

    weight_method = "sqrt"

    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    x_fit = np.linspace(min(x), max(x), 500)
    y_fit = generalized_logistic(x_fit, *popt)
    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    
    
if 'median_n' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('median_n')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed = weights[1:-1]

    weight_method = "none"

    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    # x_fit = np.linspace(min(x), max(x), 500)
    x_min, x_max = min(x), max(x)
    break_point = 200  # Adjust based on where you need more detail
    num_low = 400  # Higher density for small values
    num_high = 100  # Lower density for large values
    x_low = np.geomspace(x_min + 1e-6, break_point, num=num_low)  # Avoid zero in log-space
    x_high = np.linspace(break_point, x_max, num=num_high)  # Linear spacing for large values
    x_fit = np.concatenate([x_low, x_high])  # Combine both ranges
    
    y_fit = generalized_logistic(x_fit, *popt)
    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    
    
if 'median_n_pangolin' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('median_n_pangolin')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed = weights[1:-1]

    weight_method = "none"

    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    # x_fit = np.linspace(min(x), max(x), 500)
    x_min, x_max = min(x), max(x)
    break_point = 200  # Adjust based on where you need more detail
    num_low = 400  # Higher density for small values
    num_high = 100  # Lower density for large values
    x_low = np.geomspace(x_min + 1e-6, break_point, num=num_low)  # Avoid zero in log-space
    x_high = np.linspace(break_point, x_max, num=num_high)  # Linear spacing for large values
    x_fit = np.concatenate([x_low, x_high])  # Combine both ranges
    
    y_fit = generalized_logistic(x_fit, *popt)
    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    
    
if 'pangolin_tissue_score' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('pangolin_tissue_score')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed_all = weights[1:-1]
    x_fit = np.linspace(min(x), max(x), 1000)

    weight_method = "sqrt"

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed > 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_below_zero = generalized_logistic(x_fit, *popt)

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed < 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_above_zero = generalized_logistic(x_fit, *popt)

    y_fit = [*y_fit_below_zero[x_fit < 0], *y_fit_above_zero[x_fit >= 0]]

    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    
    
if 'delta_logit_psi' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('delta_logit_psi')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed_all = weights[1:-1]
    x_fit = np.linspace(min(x), max(x), 1000)

    weight_method = "sqrt"

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed > 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_below_zero = generalized_logistic(x_fit, *popt)

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed < 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_above_zero = generalized_logistic(x_fit, *popt)

    y_fit = [*y_fit_below_zero[x_fit < 0], *y_fit_above_zero[x_fit >= 0]]

    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    
    
if 'delta_psi' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('delta_psi')
    feature_idx = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[feature_idx]
    x = ebm_model.preprocessor_.col_bin_edges_[feature_idx]
    y = contrib_func_orig
    weights = ebm_model.preprocessor_.col_bin_counts_[feature_idx]

    # Trim data (as in your pipeline)
    x_trimmed = x
    y_trimmed = y[1:-1]
    weights_trimmed_all = weights[1:-1]
    x_fit = np.linspace(min(x), max(x), 1000)

    weight_method = "sqrt"

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed > 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_below_zero = generalized_logistic(x_fit, *popt)

    weights_trimmed = weights_trimmed_all.copy()
    weights_trimmed[x_trimmed < 0] = 0
    popt = fit_generalized_logistic_with_weights(x_trimmed, y_trimmed, weights_trimmed, weight_method=weight_method)
    y_fit_above_zero = generalized_logistic(x_fit, *popt)

    y_fit = [*y_fit_below_zero[x_fit < 0], *y_fit_above_zero[x_fit >= 0]]

    y_fit = np.array([0, *y_fit, y_fit[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_fit
    ebm_model.additive_terms_[feature_idx] = y_fit
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_fit))
    

# Save the model
# assert snakemake.output['model_postprocessed'] == snakemake.input['model'].replace('.pkl', '_FITTED_GENERALIZED_LOGISTIC.pkl')
with open(snakemake.output['model_postprocessed'], 'wb') as file:
    pickle.dump(ebm_model, file)