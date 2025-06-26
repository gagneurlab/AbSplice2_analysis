import pandas as pd
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from sklearn.isotonic import IsotonicRegression
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

def detect_trend_change(x, y_smooth):
    """
    Detects the best split point where the trend changes from increasing to decreasing or vice versa.
    The split point is chosen as the center of the region where the gradient changes sign.
    """
    gradient = np.gradient(y_smooth)
    sign_change_indices = np.where(np.diff(np.sign(gradient)) != 0)[0] + 1  # Find indices where gradient changes sign

    if len(sign_change_indices) == 0:
        return None  # No clear split point found

    # Choose the split point in the center of the trend change region
    split_index = sign_change_indices[len(sign_change_indices) // 2]

    return split_index

def smooth_monotone(x, y, increasing="auto", smoothing="gaussian", plot=False, split_index=None):
    """
    Smooths and enforces monotonicity on the input data, automatically detecting trend changes if needed.
    
    Parameters:
    - x: array-like, the x values of the data.
    - y: array-like, the y values of the data.
    - increasing: 'auto' (default), 'both', or boolean.
                  - 'auto' detects whether the first part should be increasing or decreasing automatically.
                  - 'both' assumes a piecewise trend with one increasing and one decreasing section.
    - smoothing: str, "gaussian" or "savgol" to select the smoothing technique.
    - plot: bool, whether to display a plot of the process.

    Returns:
    - x_final: array-like, the interpolated x values.
    - y_final: array-like, the interpolated monotonic smoothed values.
    """
    
    x = np.array(x)
    y = np.array(y)

    # Ignore the first point for smoothing and fitting
    x_filtered = x[1:]
    y_filtered = y[1:]

    # Step 1: Apply Smoothing
    if smoothing == "gaussian":
        y_smooth = gaussian_filter1d(y_filtered, sigma=2)
    elif smoothing == "savgol":
        y_smooth = savgol_filter(y_filtered, window_length=7, polyorder=3)
    else:
        raise ValueError("Invalid smoothing method. Choose 'gaussian' or 'savgol'.")

    # Determine if the first segment should be increasing or decreasing
    if increasing == "auto":
        # Detect the optimal split point
        split_index = detect_trend_change(x_filtered, y_smooth)
        increasing_first_part = np.gradient(y_smooth[:split_index]).mean() > 0
    elif isinstance(increasing, bool):
        increasing_first_part = increasing
    else:
        increasing_first_part = True  # Default to increasing

    # Apply isotonic regression in one or two parts
    if split_index is not None:
        # First part: Apply isotonic regression
        iso_reg_first = IsotonicRegression(increasing=increasing_first_part)
        y_mono_first = iso_reg_first.fit_transform(x_filtered[:split_index + 1], y_smooth[:split_index + 1])

        # Second part: Apply isotonic regression in the opposite direction
        iso_reg_second = IsotonicRegression(increasing=not increasing_first_part)
        # y_mono_second = iso_reg_second.fit_transform(x_filtered[split_index:], y_smooth[split_index:])
        y_mono_second = iso_reg_second.fit_transform(x_filtered[split_index + 1:], y_smooth[split_index + 1:])  # FIXED

        # Ensure x-values are strictly increasing
        x_mono_first = x_filtered[:split_index + 1]
        x_mono_second = x_filtered[split_index + 1:]

        # Interpolation for both segments
        interp_first = PchipInterpolator(x_mono_first, y_mono_first)
        interp_second = PchipInterpolator(x_mono_second, y_mono_second)

        # Merge interpolated results
        x_final = np.concatenate([x_mono_first, x_mono_second])
        y_final = np.concatenate([y_mono_first, y_mono_second])

    else:
        # Apply isotonic regression globally (either increasing or decreasing)
        iso_reg = IsotonicRegression(increasing=increasing_first_part)
        y_mono = iso_reg.fit_transform(x_filtered, y_smooth)

        # Ensure x-values are strictly increasing
        x_mono = x_filtered

        # Step 3: Interpolation for further smoothness
        interp = PchipInterpolator(x_mono, y_mono)
        x_final = x_mono
        y_final = interp(x_final)    

    # Re-add the first data point (optional)
    x_final = np.insert(x_final, 0, x[0])
    y_final = np.insert(y_final, 0, y[0])

    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(x, y, 'o', alpha=0.3, label="Original Data")
        plt.plot(x_filtered, y_smooth, '--', label=f"Smoothed ({smoothing})")
        plt.plot(x_final, y_final, '-', label="Monotonic Interpolation")
        if split_index is not None:
            plt.axvline(x[split_index], color='gray', linestyle='dashed', label="Split Point")
        plt.legend()
        plt.show()

    return y_final

# Example Usage:
# x_smooth, y_smooth_monotone = smooth_monotone(x, y, increasing="auto", smoothing="savgol", plot=True)


# Load the model from the pkl file
with open(snakemake.input['model'], "rb") as file:
    ebm_model = pickle.load(file)
    
# update the feautre functions
if 'gain_score' in ebm_model.feature_names:
    gain_score_index = ebm_model.feature_names.index('gain_score')
    i = gain_score_index
    contrib_func_orig = ebm_model.additive_terms_[i]
    x = range(len(contrib_func_orig))
    y = contrib_func_orig
    ebm_model.additive_terms_[i] = smooth_monotone(x, y, increasing=True, plot=False)
    
if 'loss_score' in ebm_model.feature_names:
    loss_score_index = ebm_model.feature_names.index('loss_score')
    i = loss_score_index
    contrib_func_orig = ebm_model.additive_terms_[i]
    x = range(len(contrib_func_orig))
    y = contrib_func_orig
    ebm_model.additive_terms_[i] = smooth_monotone(x, y, increasing=False, plot=False)
    
if 'pangolin_tissue_score' in ebm_model.feature_names:
    pangolin_tissue_score_index = ebm_model.feature_names.index('pangolin_tissue_score')
    i = pangolin_tissue_score_index
    contrib_func_orig = ebm_model.additive_terms_[i]
    x = range(len(contrib_func_orig))
    y = contrib_func_orig
    ebm_model.additive_terms_[i] = smooth_monotone(x, y, increasing='auto', plot=False)

# save postprocessed model    
with open(snakemake.output['model_postprocessed'], 'wb') as file:
    pickle.dump(ebm_model, file)