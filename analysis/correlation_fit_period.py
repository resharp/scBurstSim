# fit fluctuating period on:
#   cosine function of correlation with exponential dampening
#   for a certain efficiency of RNA retrieval
#   assuming resolution on allele level
# based on df_corr_all calculated in correlation_labels.py

import io
import os
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np

if os.name == 'nt':
    dir_sep = "\\"
    prj_dir = r"D:\26 Battich Oudenaarden transcriptional bursts"
else:
    dir_sep = "/"
    prj_dir = "sc_runs"

gap = 0
label_1 = "EU"
label_2 = "4SU"


def dampening_fluctuation(x, period, amplitude, offset, dampening):

    y = amplitude * ((np.cos(2 * np.pi * x / period)) + offset) * np.exp(-dampening * x)

    return y


def fit_period(series):

    inital_guesses = [60, 120, 180, 240, 300]
    min_delta_corr = 0.5

    min_corr = min(series["corr"])
    max_corr = max(series["corr"])

    if np.abs(min_corr - max_corr) < min_delta_corr:
        period = -1

    else:
        expected = (120, 1, -1, 0.01)

        try:
            popt, pcov = curve_fit(dampening_fluctuation, series.len_win, series["corr"], p0=expected, maxfev=1000)
            period = popt[0]
        except RuntimeError:  # this happens mostly when exceeding maxfev
            period = -2

    # if period != -1:
    #     print(period)
    df = pd.DataFrame([[period]], columns=["period"])
    return df


eff = 1

out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)
plot_dir = out_dir + dir_sep + "correlation_labels.plots"

corr_name = "{od}{dir_sep}df_corr_all_G{gap}.csv".format(
    od=plot_dir, dir_sep=dir_sep, gap=gap)

df_corr_all = pd.read_csv(corr_name, sep=';')

# only take the deterministically fluctuating genes with a period of 2 hours
df_corr_all = df_corr_all[(df_corr_all.strategy.str.endswith("_2")) & (df_corr_all.strategy.str.contains("F"))]

df_predictions = df_corr_all.groupby(['strategy']).apply(fit_period).reset_index()

debug = True


plot_dir = r"{}{}runs{}{}".format(prj_dir, dir_sep, dir_sep, "correlation_fit_period.plots")
os.makedirs(plot_dir, exist_ok=True)

# save predictions for efficiency
filename = plot_dir + dir_sep + "period_predictions_{}.csv".format(eff)
df_predictions.to_csv(filename, sep=";", index=False)
