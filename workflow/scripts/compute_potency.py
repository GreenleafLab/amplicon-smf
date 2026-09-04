# this script basically computes potency from 
# it requires as input the big table that is promoter state calls per molecule
# and the big table that is binding model calls per molecule
# and then it groupbys sample and then computes potency using the code I got from ChatGPT

import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from copy import deepcopy
import gc
from scipy.optimize import curve_fit
from scipy.special import comb
from itertools import combinations
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import scipy.stats as stats
from matplotlib.collections import PolyCollection
import statsmodels.api as sm
from statsmodels.regression.linear_model import WLS
from scipy import stats

# this is julia's additive activate model with a basal on rate that allows us to fit the y intercept
# can potentially think about trying to plot the k_basal and also the k_tf...
def saturation_model_mechanistic(n_bound, k_tf, k_basal):
    """
    Mechanistic saturation model: 
    frac_on = 1/(1 + 1/(k_basal + k_tf*n_bound))
          = (k_basal + k_tf*n_bound)/(1 + k_basal + k_tf*n_bound)
    
    When n_bound = 0: frac_on = k_basal/(1+k_basal) ≈ k_basal (basal activity)
    When n_bound → ∞: frac_on → 1
    
    Parameters:
    -----------
    n_bound : array
        Number of TF bound
    k_tf : float
        TF activity parameter (contribution per bound TF)
    k_basal : float
        Basal promoter activity (in absence of TF)
    """
    return (k_basal + k_tf * n_bound) / (1 + k_basal + k_tf * n_bound)

def calculate_saturation_regression_stats(data, group_cols, x_col='n_bound', 
                                          y_col='promoter_nuc_free', 
                                          weight_col='num_molecules',
                                          model_type='mechanistic', # Added switch
                                          k_tf_init=1.0,
                                          k_basal_init=0.01,
                                          fit_basal=True,
                                          fixed_basal=0.0):
    results = []
    
    for group_vals, df in data.groupby(group_cols):
        # Ensure we have enough data points (min 3 for slope + intercept/basal)
        min_points = 3 if (fit_basal or model_type == 'linear') else 2
        if len(df) < min_points:
            continue
            
        x_data = df[x_col].astype(float).values
        y = df[y_col].astype(float).values
        weights = df[weight_col].astype(float).values
        
        # Original data validation
        if np.any(x_data < 0) or np.any(y < 0) or np.any(y > 1):
            continue
        
        try:
            if model_type == 'linear':
                # --- NEW LINEAR WLS OPTION ---
                X = sm.add_constant(x_data)
                # statsmodels WLS uses weights as 1/variance (num_molecules is perfect)
                wls_model = sm.WLS(y, X, weights=weights)
                res = wls_model.fit()
                
                # Map linear params to your dictionary structure
                intercept, slope = res.params
                intercept_se, slope_se = res.bse
                
                # Calculate metrics for consistency
                y_pred = res.predict(X)
                r_squared = res.rsquared
                dof = len(x_data) - 2
                ss_res = np.sum(weights * (y - y_pred)**2)
                rse = np.sqrt(ss_res / dof) if dof > 0 else np.nan
                
                result_dict = {
                    'slope': slope,
                    'slope_stderr': slope_se,
                    'pvalue_slope': res.pvalues[1],
                    'intercept': intercept,
                    'intercept_stderr': intercept_se,
                    'baseline': intercept, # At n_bound=0, y = intercept
                    'r_squared': r_squared,
                    'residual_se': rse,
                    'model_type': 'linear'
                }

            else:
                # --- ORIGINAL MECHANISTIC OPTION (Restored all lines) ---
                sigma = 1 / np.sqrt(weights)
                if fit_basal:
                    popt, pcov = curve_fit(
                        saturation_model_mechanistic, x_data, y, 
                        p0=[k_tf_init, k_basal_init], sigma=sigma,
                        absolute_sigma=False, bounds=([0, 0], [np.inf, np.inf]),
                        maxfev=10000
                    )
                    k_tf_fit, k_basal_fit = popt
                    k_tf_stderr, k_basal_stderr = np.sqrt(np.diag(pcov))
                    k_basal_fixed, n_params = False, 2
                    
                    dof = len(x_data) - n_params
                    pvalue_k_tf = 2 * (1 - stats.t.cdf(np.abs(k_tf_fit / k_tf_stderr), dof))
                    pvalue_k_basal = 2 * (1 - stats.t.cdf(np.abs(k_basal_fit / k_basal_stderr), dof))
                else:
                    # (Your fixed basal logic stays here)
                    pass 

                # Restored original specific metrics
                baseline = k_basal_fit / (1 + k_basal_fit)
                k_half = (1 - k_basal_fit) / k_tf_fit if k_tf_fit > 0 else np.inf
                y_pred = saturation_model_mechanistic(x_data, k_tf_fit, k_basal_fit)
                
                ss_res = np.sum(weights * (y - y_pred)**2)
                ss_tot = np.sum(weights * (y - np.average(y, weights=weights))**2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                rse = np.sqrt(ss_res / dof) if dof > 0 else np.nan

                result_dict = {
                    'k_tf': k_tf_fit, 'k_tf_stderr': k_tf_stderr, 'pvalue_k_tf': pvalue_k_tf,
                    'k_basal': k_basal_fit, 'k_basal_stderr': k_basal_stderr, 'pvalue_k_basal': pvalue_k_basal,
                    'baseline': baseline, 'k_half': k_half, 'r_squared': r_squared,
                    'residual_se': rse, 'model_type': 'mechanistic'
                }

            # Shared Metadata (Restored)
            result_dict.update({
                'n_points': len(df),
                'total_molecules': weights.sum(),
                'min_molecules': weights.min(),
                'max_molecules': weights.max()
            })
            
            # Add group identifiers
            if isinstance(group_vals, tuple):
                for col, val in zip(group_cols, group_vals): result_dict[col] = val
            else:
                result_dict[group_cols[0]] = group_vals
            results.append(result_dict)
            
        except Exception as e:
            print(f"Error fitting group {group_vals}: {e}")
            
    return pd.DataFrame(results)

# def calculate_saturation_regression_stats(data, group_cols, x_col='n_bound', 
#                                           y_col='promoter_nuc_free', 
#                                           weight_col='num_molecules',
#                                           k_tf_init=1.0,
#                                           k_basal_init=0.01,
#                                           fit_basal=True,
#                                           fixed_basal=0.0):
#     """
#     Perform weighted nonlinear regression for mechanistic saturation model.
    
#     Parameters:
#     -----------
#     data : pd.DataFrame
#         Aggregated data with means and molecule counts
#     group_cols : list
#         Columns to group by (e.g., ['library', 'sample', 'rep'])
#     x_col : str
#         Column name for x variable (n_bound)
#     y_col : str
#         Column name for y variable (promoter_nuc_free)
#     weight_col : str
#         Column name for weights (number of molecules)
#     k_tf_init : float
#         Initial guess for k_tf parameter (TF activity)
#     k_basal_init : float
#         Initial guess for k_basal parameter (basal activity)
#     fit_basal : bool
#         If True, fit both k_tf and k_basal. If False, fix k_basal at fixed_basal.
#     fixed_basal : float
#         Value to fix k_basal at if fit_basal=False (default 0.0)
        
#     Returns:
#     --------
#     pd.DataFrame with regression results for each group
#     """
    
#     results = []
    
#     for group_vals, df in data.groupby(group_cols):
#         # Ensure we have enough data points
#         min_points = 3 if fit_basal else 2
#         if len(df) < min_points:
#             continue
            
#         # Extract values
#         x_data = df[x_col].astype(float).values
#         y = df[y_col].astype(float).values
#         weights = df[weight_col].astype(float).values
        
#         # Skip if we have invalid data
#         if np.any(x_data < 0) or np.any(y < 0) or np.any(y > 1):
#             print(f"Warning: Invalid data in group {group_vals}, skipping")
#             continue
        
#         try:
#             # Weighted nonlinear least squares using curve_fit
#             sigma = 1 / np.sqrt(weights)
            
#             if fit_basal:
#                 # Fit both k_tf and k_basal
#                 popt, pcov = curve_fit(
#                     saturation_model_mechanistic, 
#                     x_data, 
#                     y, 
#                     p0=[k_tf_init, k_basal_init],
#                     sigma=sigma,
#                     absolute_sigma=False,
#                     bounds=([0, 0], [np.inf, np.inf]),  # Both must be positive
#                     maxfev=10000
#                 )
#                 k_tf_fit = popt[0]
#                 k_basal_fit = popt[1]
#                 k_tf_stderr = np.sqrt(np.diag(pcov))[0]
#                 k_basal_stderr = np.sqrt(np.diag(pcov))[1]
#                 k_basal_fixed = False
#                 n_params = 2
                
#                 # Calculate p-values
#                 dof = len(x_data) - n_params
#                 t_stat_k_tf = k_tf_fit / k_tf_stderr if k_tf_stderr > 0 else np.inf
#                 pvalue_k_tf = 2 * (1 - stats.t.cdf(np.abs(t_stat_k_tf), dof))
#                 t_stat_k_basal = k_basal_fit / k_basal_stderr if k_basal_stderr > 0 else np.inf
#                 pvalue_k_basal = 2 * (1 - stats.t.cdf(np.abs(t_stat_k_basal), dof))
                
#             else:
#                 # Fix k_basal, only fit k_tf
#                 def fixed_basal_model(n_bound, k_tf):
#                     return saturation_model_mechanistic(n_bound, k_tf, fixed_basal)
                
#                 popt, pcov = curve_fit(
#                     fixed_basal_model, 
#                     x_data, 
#                     y, 
#                     p0=[k_tf_init],
#                     sigma=sigma,
#                     absolute_sigma=False,
#                     bounds=(0, np.inf),
#                     maxfev=10000
#                 )
#                 k_tf_fit = popt[0]
#                 k_basal_fit = fixed_basal
#                 k_tf_stderr = np.sqrt(np.diag(pcov))[0]
#                 k_basal_stderr = 0.0
#                 k_basal_fixed = True
#                 n_params = 1
                
#                 # Calculate p-value
#                 dof = len(x_data) - n_params
#                 t_stat_k_tf = k_tf_fit / k_tf_stderr if k_tf_stderr > 0 else np.inf
#                 pvalue_k_tf = 2 * (1 - stats.t.cdf(np.abs(t_stat_k_tf), dof))
#                 pvalue_k_basal = np.nan
            
#             # Calculate baseline (frac_on when n_bound = 0)
#             baseline = k_basal_fit / (1 + k_basal_fit)
            
#             # Calculate half-maximal n_bound
#             # When does saturation reach 0.5? When k_basal + k_tf*n = 1
#             # So n_half = (1 - k_basal)/k_tf
#             k_half = (1 - k_basal_fit) / k_tf_fit if k_tf_fit > 0 else np.inf
#             # Propagate error (approximate)
#             if fit_basal:
#                 # This is approximate - proper error propagation would use full covariance
#                 k_half_stderr = k_half * np.sqrt((k_tf_stderr/k_tf_fit)**2 + (k_basal_stderr/(1-k_basal_fit))**2)
#             else:
#                 k_half_stderr = k_tf_stderr / (k_tf_fit**2)
            
#             # Calculate fitted values and residuals
#             y_pred = saturation_model_mechanistic(x_data, k_tf_fit, k_basal_fit)
#             residuals = y - y_pred
            
#             # Weighted sum of squares
#             ss_res = np.sum(weights * residuals**2)
#             ss_tot = np.sum(weights * (y - np.average(y, weights=weights))**2)
#             r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
#             # Residual standard error
#             rse = np.sqrt(ss_res / dof) if dof > 0 else np.nan
            
#             # Store results
#             result_dict = {
#                 'k_tf': k_tf_fit,
#                 'k_tf_stderr': k_tf_stderr,
#                 'pvalue_k_tf': pvalue_k_tf,
#                 'k_basal': k_basal_fit,
#                 'k_basal_stderr': k_basal_stderr,
#                 'k_basal_fixed': k_basal_fixed,
#                 'pvalue_k_basal': pvalue_k_basal,
#                 'baseline': baseline,  # k_basal/(1+k_basal)
#                 'k_half': k_half,
#                 'k_half_stderr': k_half_stderr,
#                 'r_squared': r_squared,
#                 'residual_se': rse,
#                 'n_points': len(df),
#                 'total_molecules': weights.sum(),
#                 'min_molecules': weights.min(),
#                 'max_molecules': weights.max()
#             }
            
#             # Add group identifiers
#             if isinstance(group_vals, tuple):
#                 for col, val in zip(group_cols, group_vals):
#                     result_dict[col] = val
#             else:
#                 result_dict[group_cols[0]] = group_vals
                
#             results.append(result_dict)
            
#         except Exception as e:
#             print(f"Error fitting group {group_vals}: {e}")
#             continue
    
#     return pd.DataFrame(results)


def analyze_promoter_binding_saturation(all_promoter_annotations, mol_thresh=10, 
                                        background=None, model_type='mechanistic', k_tf_init=1.0, k_basal_init=0.01,
                                        fit_basal=True, fixed_basal=0.0):
    """
    Complete analysis pipeline for promoter binding data with mechanistic saturation model.
    
    Parameters:
    -----------
    mol_thresh : int
        Minimum number of molecules required per data point
    background : str or None
        If provided (e.g., 'b0'), filter to this background. If None, use all data.
    k_tf_init : float
        Initial guess for k_tf (TF activity parameter)
    k_basal_init : float
        Initial guess for k_basal (basal activity parameter)
    fit_basal : bool
        If True, fit both k_tf and k_basal. If False, fix k_basal at fixed_basal.
    fixed_basal : float
        Value to fix k_basal at if fit_basal=False (default 0.0 for no basal activity)
    """
    
    # Filter data by background if specified
    if background is not None:
        data = all_promoter_annotations.loc[
            all_promoter_annotations.background == background
        ].copy()
    else:
        data = all_promoter_annotations.copy()
    
    # Aggregate to get averages per amplicon
    # Keep amplicons separate - don't group by n_tfbs, but by the actual amplicon ID
    # Assuming you have a column like 'promoter' or 'amplicon' that uniquely identifies each amplicon
    data_agg = data.groupby(
        ['sample', 'amplicon']  # or whatever your amplicon ID column is
    ).agg({
        'n_bound': np.mean,
        'promoter_nuc_free': np.mean,
        'tf_bound': len,  # count molecules
        'n_tfbs': 'first'  # Keep n_tfbs as metadata (should be constant within amplicon)
    }).reset_index()
    
    data_agg = data_agg.rename({'tf_bound': 'num_molecules'}, axis=1)
    
    # Filter by minimum molecules
    data_agg_filtered = data_agg[data_agg.num_molecules > mol_thresh].copy()
    
    # Perform regression
    regression_results = calculate_saturation_regression_stats(
        data_agg_filtered,
        group_cols=['sample'],
        x_col='n_bound',
        y_col='promoter_nuc_free',
        weight_col='num_molecules',
        model_type=model_type,
        k_tf_init=k_tf_init,
        k_basal_init=k_basal_init,
        fit_basal=fit_basal,
        fixed_basal=fixed_basal
    )
    
    return data_agg_filtered, regression_results


def plot_saturation_results(data_agg, results, 
                            figsize_per_plot=(4, 4),
                            ncols=3,
                            show_ci=True,
                            alpha_points=0.6,
                            show_weights=True,
                            n_sample=None,
                            max_plots_per_figure=30,
                            random_seed=42):
    """
    Plot regression results for both mechanistic and linear models.
    Automatically adjusts labels and CI based on the 'model_type' column.
    """
    # Sample or use all results
    if n_sample is not None:
        np.random.seed(random_seed)
        n_groups = min(n_sample, len(results))
        plot_results = results.sample(n=n_groups, random_state=random_seed)
        print(f"Randomly sampled {n_groups} groups out of {len(results)} total")
    else:
        plot_results = results
        n_groups = len(results)
        print(f"Plotting all {n_groups} groups")
    
    # Split into multiple figures if needed
    n_figures = int(np.ceil(n_groups / max_plots_per_figure))
    all_figs = []
    
    for fig_idx in range(n_figures):
        start_idx = fig_idx * max_plots_per_figure
        end_idx = min((fig_idx + 1) * max_plots_per_figure, n_groups)
        n_plots_this_fig = end_idx - start_idx
        
        nrows = int(np.ceil(n_plots_this_fig / ncols))
        
        fig, axes = plt.subplots(nrows, ncols, 
                                figsize=(figsize_per_plot[0] * ncols, 
                                         figsize_per_plot[1] * nrows),
                                squeeze=False)
        axes = axes.flatten()
        
        for plot_idx, (_, row) in enumerate(plot_results.iloc[start_idx:end_idx].iterrows()):
            ax = axes[plot_idx]
            m_type = row.get('model_type', 'mechanistic') # Default to mechanistic for safety
            
            # Get group identifiers
            group_cols = [col for col in results.columns 
                         if col in ['library', 'sample', 'rep']]
            group_dict = {col: row[col] for col in group_cols}
            
            # Filter data for this group
            mask = np.ones(len(data_agg), dtype=bool)
            for col, val in group_dict.items():
                mask &= (data_agg[col] == val)
            df_group = data_agg[mask]
            
            if len(df_group) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            # Extract data
            x = df_group['n_bound'].values
            y = df_group['promoter_nuc_free'].values
            weights = df_group['num_molecules'].values
            
            # Plot data points
            if show_weights:
                sizes = weights / weights.max() * 200 + 20
                scatter = ax.scatter(x, y, s=sizes, alpha=alpha_points, 
                                   c=weights, cmap='viridis', edgecolors='k', linewidth=0.5)
                cbar = fig.colorbar(scatter, ax=ax)
                cbar.set_label('# molecules', fontsize=8)
            else:
                ax.scatter(x, y, alpha=alpha_points, edgecolors='k', linewidth=0.5)
            
            # Create x range for plotting regression curve
            x_min, x_max = 0, x.max()
            x_range = np.linspace(x_min, x_max, 200)
            
            # Model-specific plotting and stats formatting
            if m_type == 'linear':
                # Linear Fit
                y_fit = row['slope'] * x_range + row['intercept']
                
                if show_ci and row['slope_stderr'] > 0:
                    # Simple linear CI (1.96 * SE)
                    y_fit_upper = (row['slope'] + 1.96 * row['slope_stderr']) * x_range + row['intercept']
                    y_fit_lower = (row['slope'] - 1.96 * row['slope_stderr']) * x_range + row['intercept']
                    ax.fill_between(x_range, np.clip(y_fit_lower, 0, 1), np.clip(y_fit_upper, 0, 1), 
                                    alpha=0.2, color='red', label='95% CI')
                
                textstr = '\n'.join([
                    f'Slope = {row["slope"]:.3f} ± {row["slope_stderr"]:.3f}',
                    f'Intercept = {row["intercept"]:.3f}',
                    f'R² = {row["r_squared"]:.3f}',
                    f'p(slope) = {row["pvalue_slope"]:.2e}',
                    f'n = {row["n_points"]}'
                ])
            else:
                # Mechanistic Fit
                y_fit = saturation_model_mechanistic(x_range, row['k_tf'], row['k_basal'])
                
                if show_ci and row['k_tf_stderr'] > 0:
                    y_fit_upper = saturation_model_mechanistic(x_range, row['k_tf'] + 1.96*row['k_tf_stderr'], row['k_basal'])
                    y_fit_lower = saturation_model_mechanistic(x_range, row['k_tf'] - 1.96*row['k_tf_stderr'], row['k_basal'])
                    ax.fill_between(x_range, np.clip(y_fit_lower, 0, 1), np.clip(y_fit_upper, 0, 1), 
                                    alpha=0.2, color='red', label='95% CI')
                
                textstr = '\n'.join([
                    f'k_tf = {row["k_tf"]:.3f} ± {row["k_tf_stderr"]:.3f}',
                    f'k_basal = {row["k_basal"]:.3f}',
                    f'R² = {row["r_squared"]:.3f}',
                    f'n = {row["n_points"]}'
                ])

            # Draw Fit Line and Baseline
            ax.plot(x_range, y_fit, 'r-', linewidth=2, label='Fit')
            ax.axhline(row['baseline'], color='blue', linestyle='--', 
                      linewidth=1, alpha=0.5, label=f'Baseline={row["baseline"]:.3f}')
            
            # Title and Text Box
            title = ', '.join([f"{col}={val}" for col, val in group_dict.items()])
            ax.set_title(title, fontsize=9, fontweight='bold')
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=7, 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            ax.set_xlabel('n_bound (mean)', fontsize=10)
            ax.set_ylabel('promoter_nuc_free (mean)', fontsize=10)
            ax.set_ylim(-0.05, 1.05)
            ax.legend(fontsize=7, loc='lower right')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(n_plots_this_fig, len(axes)):
            axes[idx].axis('off')
        
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        all_figs.append((fig, axes))
        print(f"Figure {fig_idx + 1}/{n_figures}: plots {start_idx} to {end_idx-1}")
    
    return all_figs

    # def plot_saturation_results(data_agg, results, 
    #                         figsize_per_plot=(4, 4),
    #                         ncols=3,
    #                         show_ci=True,
    #                         alpha_points=0.6,
    #                         show_weights=True,
    #                         n_sample=None,
    #                         max_plots_per_figure=30,
    #                         random_seed=42):
    # """
    # Plot saturation regression results, split into multiple figures if needed.
    
    # Parameters:
    # -----------
    # data_agg : pd.DataFrame
    #     Aggregated data used for regression
    # results : pd.DataFrame
    #     Regression results from calculate_saturation_regression_stats
    # figsize_per_plot : tuple
    #     Size of each subplot
    # ncols : int
    #     Number of columns in the grid
    # show_ci : bool
    #     Whether to show confidence interval on regression curve
    # alpha_points : float
    #     Transparency of data points
    # show_weights : bool
    #     Whether to size points by number of molecules
    # n_sample : int or None
    #     Number of random groups to plot. If None, plot all groups.
    # max_plots_per_figure : int
    #     Maximum number of subplots per figure
    # random_seed : int
    #     Random seed for reproducibility when sampling
        
    # Returns:
    # --------
    # list of (fig, axes) tuples - one per figure created
    # """
    
    # # Sample or use all results
    # if n_sample is not None:
    #     np.random.seed(random_seed)
    #     n_groups = min(n_sample, len(results))
    #     plot_results = results.sample(n=n_groups, random_state=random_seed)
    #     print(f"Randomly sampled {n_groups} groups out of {len(results)} total")
    # else:
    #     plot_results = results
    #     n_groups = len(results)
    #     print(f"Plotting all {n_groups} groups")
    
    # # Split into multiple figures if needed
    # n_figures = int(np.ceil(n_groups / max_plots_per_figure))
    # all_figs = []
    
    # for fig_idx in range(n_figures):
    #     start_idx = fig_idx * max_plots_per_figure
    #     end_idx = min((fig_idx + 1) * max_plots_per_figure, n_groups)
    #     n_plots_this_fig = end_idx - start_idx
        
    #     nrows = int(np.ceil(n_plots_this_fig / ncols))
        
    #     fig, axes = plt.subplots(nrows, ncols, 
    #                             figsize=(figsize_per_plot[0] * ncols, 
    #                                     figsize_per_plot[1] * nrows),
    #                             squeeze=False)
    #     axes = axes.flatten()
        
    #     for plot_idx, (_, row) in enumerate(plot_results.iloc[start_idx:end_idx].iterrows()):
    #         ax = axes[plot_idx]
            
    #         # Get group identifiers
    #         group_cols = [col for col in results.columns 
    #                      if col in ['library', 'sample', 'rep']]
    #         group_dict = {col: row[col] for col in group_cols}
            
    #         # Filter data for this group
    #         mask = np.ones(len(data_agg), dtype=bool)
    #         for col, val in group_dict.items():
    #             mask &= (data_agg[col] == val)
    #         df_group = data_agg[mask]
            
    #         if len(df_group) == 0:
    #             ax.text(0.5, 0.5, 'No data', ha='center', va='center', 
    #                    transform=ax.transAxes)
    #             ax.set_xlabel('n_bound (mean)', fontsize=10)
    #             ax.set_ylabel('promoter_nuc_free (mean)', fontsize=10)
    #             continue
            
    #         # Extract data
    #         x = df_group['n_bound'].values
    #         y = df_group['promoter_nuc_free'].values
    #         weights = df_group['num_molecules'].values
            
    #         # Plot data points
    #         if show_weights:
    #             sizes = weights / weights.max() * 200 + 20
    #             scatter = ax.scatter(x, y, s=sizes, alpha=alpha_points, 
    #                                c=weights, cmap='viridis', edgecolors='k', linewidth=0.5)
    #             cbar = fig.colorbar(scatter, ax=ax)
    #             cbar.set_label('# molecules', fontsize=8)
    #         else:
    #             ax.scatter(x, y, alpha=alpha_points, edgecolors='k', linewidth=0.5)
            
    #         # Create x range for plotting regression curve
    #         x_min = 0
    #         x_max = x.max()
    #         x_range = np.linspace(x_min, x_max, 200)
    #         # y_fit = saturation_model_mechanistic(x_range, row['k_tf'], row['k_basal'])
    #         if row.get('model_type') == 'linear':
    #             y_fit = row['slope'] * x_range + row['intercept']
    #         else:
    #             y_fit = saturation_model_mechanistic(x_range, row['k_tf'], row['k_basal'])
            
    #         # Plot regression curve
    #         ax.plot(x_range, y_fit, 'r-', linewidth=2, label='Fit')
            
    #         # Plot baseline as dashed line
    #         ax.axhline(row['baseline'], color='blue', linestyle='--', 
    #                   linewidth=1, alpha=0.5, label=f'Baseline={row["baseline"]:.3f}')
            
    #         # Add confidence interval (simplified - just for k_tf uncertainty)
    #         if show_ci and row['k_tf_stderr'] > 0:
    #             # This is approximate - doesn't account for k_basal uncertainty or covariance
    #             # Upper and lower bounds by varying k_tf
    #             y_fit_upper = saturation_model_mechanistic(x_range, 
    #                                                        row['k_tf'] + 1.96*row['k_tf_stderr'], 
    #                                                        row['k_basal'])
    #             y_fit_lower = saturation_model_mechanistic(x_range, 
    #                                                        row['k_tf'] - 1.96*row['k_tf_stderr'], 
    #                                                        row['k_basal'])
    #             ax.fill_between(x_range, 
    #                            np.clip(y_fit_lower, 0, 1), 
    #                            np.clip(y_fit_upper, 0, 1),
    #                            alpha=0.2, color='red', label='95% CI')
            
    #         # Add title with stats
    #         title_parts = [f"{col}={val}" for col, val in group_dict.items()]
    #         title = ', '.join(title_parts)
    #         ax.set_title(title, fontsize=9, fontweight='bold')
            
    #         # Add text box with regression stats
    #         if row['k_basal_fixed']:
    #             textstr = '\n'.join([
    #                 f'k_tf = {row["k_tf"]:.3f} ± {row["k_tf_stderr"]:.3f}',
    #                 f'k_basal = {row["k_basal"]:.3f} (fixed)',
    #                 f'R² = {row["r_squared"]:.3f}',
    #                 f'p = {row["pvalue_k_tf"]:.2e}',
    #                 f'n = {row["n_points"]}'
    #             ])
    #         else:
    #             textstr = '\n'.join([
    #                 f'k_tf = {row["k_tf"]:.3f} ± {row["k_tf_stderr"]:.3f}',
    #                 f'k_basal = {row["k_basal"]:.3f} ± {row["k_basal_stderr"]:.3f}',
    #                 f'R² = {row["r_squared"]:.3f}',
    #                 f'p(k_tf) = {row["pvalue_k_tf"]:.2e}',
    #                 f'n = {row["n_points"]}'
    #             ])
            
    #         ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
    #                fontsize=7, verticalalignment='top',
    #                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
    #         ax.set_xlabel('n_bound (mean)', fontsize=10)
    #         ax.set_ylabel('promoter_nuc_free (mean)', fontsize=10)
    #         ax.set_ylim(-0.05, 1.05)
    #         ax.legend(fontsize=7, loc='lower right')
    #         ax.grid(True, alpha=0.3)
        
    #     # Hide unused subplots
    #     for idx in range(n_plots_this_fig, len(axes)):
    #         axes[idx].axis('off')
        
    #     fig.subplots_adjust(hspace=0.4, wspace=0.4)
        
    #     all_figs.append((fig, axes))
        
    #     print(f"Figure {fig_idx + 1}/{n_figures}: plots {start_idx} to {end_idx-1}")
    
    # return all_figs



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assigns promoter state to a single-molecule matrix given a model')
    parser.add_argument("--promoter_annotation_input", type=str, help="Path to file with promoter annotations per molecule")
    parser.add_argument("--binding_model_input", type=str, help="Path to file with binding model annotations per molecule")
    parser.add_argument("--background", type=str, default=None, help="What background to fit to (or None for all)")
    parser.add_argument("--output_dir", type=str, help="Output director for file")
    parser.add_argument("--plot_dir", type=str, default=None, help="Path to plot folder for output plots")
    parser.add_argument("--out_prefix", type=str, help="Output prefix for file")
    parser.add_argument("--model_type", type=str, default="mechanistic", choices=["mechanistic", "linear"])

    args = parser.parse_args()

    # load in the promoter annotations
    promoter_annotations = pd.read_table(args.promoter_annotation_input, index_col=0)
    promoter_annotations['promoter_nuc_free'] = ~promoter_annotations['nuc_over_promoter']

    # load in the binding model collated calls
    binding_model_annotations = pd.read_table(args.binding_model_input, index_col=0)

    # join these tables
    # print(promoter_annotations.columns.tolist())
    # print(promoter_annotations.index.name)
    # print(binding_model_annotations.columns.tolist())
    # print(binding_model_annotations.index.name)
    all_molecule_annotations = promoter_annotations.join(binding_model_annotations, rsuffix='_')

    # fit the potency modelo
    data_agg, results = analyze_promoter_binding_saturation(
        all_molecule_annotations, 
        mol_thresh=10,
        model_type=args.model_type, # Pass this through
        background=args.background,  # or None for all data
        fit_basal=True
    )   

    # write the output to a file
    results.to_csv(path.join(args.output_dir, f'{args.out_prefix}_potency_stats.txt'), sep='\t', header=True, index=False)

    # Plot results
    figs = plot_saturation_results(data_agg, results, 
                                    ncols=3, show_ci=True, 
                                    n_sample=None,
                                    max_plots_per_figure=20)

    # Show all figures
    for idx, (fig, _) in enumerate(figs):
        fig.savefig(path.join(args.plot_dir, f'{args.out_prefix}_potency_plot.{idx}.pdf'))
        plt.close(fig)
