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


shift = int(snakemake.wildcards['shift'])

# Load the model from the pkl file
with open(snakemake.input['model'], "rb") as file:
    ebm_model = pickle.load(file)

# update the feautre functions
if 'median_n_pangolin' in ebm_model.feature_names:
    # Replace with your actual EBM model data extraction:
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

    interp_func = interp1d(x_trimmed - shift, y_trimmed, bounds_error=False, fill_value="extrapolate")
    
    # Calculate shifted y-values
    y_shifted = interp_func(x_trimmed)
    x_shifted = x

    y_shifted = np.array([0, *y_shifted, y_shifted[-1]])

    ebm_model.preprocessor_.col_bin_edges_[feature_idx] = x_shifted
    ebm_model.additive_terms_[feature_idx] = y_shifted
    ebm_model.term_standard_deviations_[feature_idx] = np.tile(np.mean(ebm_model.term_standard_deviations_[feature_idx]), len(y_shifted))    
    
# save model   
with open(snakemake.output['model_postprocessed'], 'wb') as file:
    pickle.dump(ebm_model, file)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



