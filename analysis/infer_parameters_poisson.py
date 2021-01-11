# here we fit a Poisson distribution to simulated data of "real" transcripts
# based on time dependent counts in one labeling window
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame
from scipy.optimize import curve_fit
from scipy.stats import poisson

from analysis.SimulatedDistribution import *
from simulator.StrategyReader import StrategyReader

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "infer_parameters_poisson.plots"
os.makedirs(plot_dir, exist_ok=True)

logger = logging.getLogger(__name__)

logfile_name = out_dir + dir_sep + 'infer_parameters_poisson.log'
logging.basicConfig(filename=logfile_name, filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
print("logging in: {file_name}".format(file_name=logfile_name))

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

nr_cells = 3000

gap = 0

label_2 = "4SU"

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


def find_fitting_optimum_to_poisson_dis(x, y_norm, mu):

    # we do not start at 0 or 1 because somehow it gives a very small mean
    expected_guesses = np.linspace(2, mu * 1.5, 5)
    cov = []
    means = []
    for expected in expected_guesses:
        popt, pcov = curve_fit(poisson.pmf, x, y_norm, p0=expected, maxfev=3000)

        estimated_mean = round_sig(popt[0], 3)

        cov.append(pcov[0, 0])
        means.append(estimated_mean)

    # find minimum variance for guesses
    var, idx = min((var, idx) for (idx, var) in enumerate(cov))

    estimated_mean = round_sig(means[idx], 3)   # and related mean
    var = round_sig(var, 3)

    return estimated_mean, var


def plot_fit_to_poisson(x, y_norm, var, estimated_mean, strategy, mu):

    plt.title("{strategy};optimal fitting;variance={variance}".
              format(strategy=strategy, mean=estimated_mean, variance=var, real_mean=mu))
    plt.plot(x, poisson.pmf(x, estimated_mean), label="Poisson, mean={}".format(estimated_mean))
    plt.plot(x, y_norm, label="real, mean={}".format(mu))
    plt.legend()
    plot_name = plot_dir + dir_sep + "fit_to_poisson_{strategy}_{nr_cells}_{len_win}.svg".format(
        strategy=strategy, nr_cells=nr_cells, len_win=len_win)
    plt.savefig(plot_name)
    plt.close(1)


# determine fitting parameters for all strategies for this window length
def fit_poisson_for_len_win(len_win, create_plot=False):

    fitted_parameters = []

    filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(len_win=len_win, gap=gap)

    df_counts = pd.read_csv(filename_counts, sep=';')

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

    sim_dis = SimulatedDistribution(df_counts, nr_cells)

    for strategy in strategies:

        df_distribution, real_mean = sim_dis.create(strategy, label_2)

        # exclude first row from fitting. The first row includes the zero counts
        df_distribution = df_distribution.iloc[1:]

        x = df_distribution["real_count"]
        y_norm = df_distribution["chance"]

        estimated_poisson_mean, var = find_fitting_optimum_to_poisson_dis(x, y_norm, real_mean)

        params = sr.get(strategy=strategy)

        calculated_mean = round_sig((params.k_on / (params.k_on + params.k_off)) * len_win * params.k_syn, 4)

        fitted_parameters.append((len_win, strategy, estimated_poisson_mean, real_mean, calculated_mean))

        if create_plot:
            plot_fit_to_poisson(x, y_norm, var, estimated_poisson_mean, strategy, real_mean)

    df_return = DataFrame(fitted_parameters,
                          columns=["len_win", "strategy", "estimated_poisson_mean", "real_mean", "calculated_mean"])

    return df_return


csv_name = plot_dir + dir_sep + "parameter_fits.csv"
logger.info("Infer parameters based on time dependent distributions from {} cells".format(nr_cells))
logger.info("Results in {}".format(csv_name))
lengths_window = [15, 30, 45, 60, 120, 180, 240]

list_df_fitted_params = []

for len_win in lengths_window:
    logger.info("start fitting for data from window length={len_win}".format(len_win=len_win))
    df_fitted_params = fit_poisson_for_len_win(len_win, create_plot=False)
    list_df_fitted_params.append(df_fitted_params)

df_fitted_params_all_lengths = pd.concat(list_df_fitted_params)

df_fitted_params_all_lengths.sort_values(by=['strategy', 'len_win'], inplace=True)

df_fitted_params_all_lengths.to_csv(csv_name, sep=";", index=False)





