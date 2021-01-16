# here we fit a Poisson distribution to simulated data of "real" transcripts
# based on time dependent counts in one labeling window, and multiple window lengths
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame
from scipy.optimize import curve_fit
from scipy.stats import poisson
from scipy.stats import linregress
import seaborn as sns

from analysis.SimulatedDistribution import *
from simulator.StrategyReader import StrategyReader
from solution.stationary import p_stationary_dimensionless

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

nr_cells = 300

gap = 0

label_2 = "4SU"

sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
# sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
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


def find_fitting_to_stationary_distribution(x, y_norm):

    # TODO: How do you vary expected when having three parameters?
    expected = (0.1, 0.1, 10)

    popt, pcov = curve_fit(p_stationary_dimensionless, x, y_norm, p0=expected, maxfev=1000)

    # dimensionless parameters (divided by k_d)
    k_on_d = popt[0]
    k_off_d = popt[1]
    k_syn_d = popt[2]

    return k_on_d, k_off_d, k_syn_d


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
# fit to
# (i) Poisson
# (ii) stationary distribution
# also determine zero fraction, real mean, calculated mean
def fit_distribution_for_len_win(len_win, strategies, create_plot=False):

    fitted_parameters = []

    filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(len_win=len_win, gap=gap)

    df_counts = pd.read_csv(filename_counts, sep=';')

    sim_dis = SimulatedDistribution(df_counts, nr_cells)

    for strategy in strategies:

        # real mean is determined here because it is used in the fitting
        df_distribution, real_mean = sim_dis.create(strategy, label_2)

        # this seems to be a very complicated way to determine the zero fraction
        # it is easier to do a set-based operation with a group by
        if len(df_distribution) > 0:
            zero_count = df_distribution.head(1).cell_id.item()
            zero_fraction = zero_count / nr_cells

        # exclude first row from fitting. The first row includes the zero counts
        # TODO: what about the fit to the stationary distribution? "Powerlaw" seems to fit better when including zeroes
        df_distribution = df_distribution.iloc[1:]

        x = df_distribution["real_count"]
        y_norm = df_distribution["chance"]

        # fit to Poisson distribution
        estimated_poisson_mean, var = find_fitting_optimum_to_poisson_dis(x, y_norm, real_mean)

        # calculations that can be moved to set-based operations
        params = sr.get(strategy=strategy)

        calculated_mean = (params.k_on / (params.k_on + params.k_off)) * len_win * params.k_syn

        error_estimated_poisson_mean = np.abs(calculated_mean - estimated_poisson_mean) / calculated_mean
        error_real_mean = np.abs(calculated_mean - real_mean) / calculated_mean
        error_estimated_poisson_mean = round_sig(error_estimated_poisson_mean, 4)
        error_real_mean = round_sig(error_real_mean, 4)

        calculated_mean = round_sig(calculated_mean, 4)  # round here, to prevent divide by zero
        # end calculations that can be moved to set-based operations

        # fit to stationary distribution
        k_on_d_fit, k_off_d_fit, k_syn_d_fit = find_fitting_to_stationary_distribution(x, y_norm)

        # calculations that can be moved to set-based operations
        k_on_d = params.k_on/params.k_d
        k_off_d = params.k_off/params.k_d
        k_syn_d = params.k_syn/params.k_d

        error_k_on_d = np.abs(k_on_d - k_on_d_fit) / k_on_d
        error_k_off_d = np.abs(k_off_d - k_off_d_fit) / k_off_d
        error_k_syn_d = np.abs(k_syn_d - k_syn_d_fit) / k_syn_d
        error_stationary_fit = np.sqrt(error_k_on_d**2 + error_k_off_d**2 + error_k_syn_d**2)

        error_k_on_d = round_sig(error_k_on_d, 4)
        error_k_off_d = round_sig(error_k_off_d, 4)
        error_k_syn_d = round_sig(error_k_syn_d, 4)
        error_stationary_fit = round_sig(error_stationary_fit, 4)
        # end calculations that can be moved to set-based operations

        fitted_parameters.append((len_win, strategy,
                                  estimated_poisson_mean, var, error_estimated_poisson_mean,
                                  real_mean, calculated_mean, error_real_mean,
                                  k_on_d_fit, k_off_d_fit, k_syn_d_fit,
                                  error_k_on_d, error_k_off_d, error_k_syn_d, error_stationary_fit,
                                  zero_fraction))

        if create_plot:
            plot_fit_to_poisson(x, y_norm, var, estimated_poisson_mean, strategy, real_mean)

    df_return = DataFrame(fitted_parameters,
                          columns=["len_win", "strategy",
                                   "estimated_poisson_mean", "variance", "error_estimated_poisson_mean",
                                   "real_mean", "calculated_mean", "error_real_mean",
                                   "k_on_d_fit", "k_off_d_fit", "k_syn_d_fit",
                                   "error_k_on_d", "error_k_off_d", "error_k_syn_d", "error_stationary_fit",
                                   "zero_fraction"])

    return df_return


def infer_parameters(lengths_window, strategies):

    logger.info("Infer parameters based on time dependent distributions from {} cells".format(nr_cells))
    logger.info("Results in {}".format(csv_name))

    list_df_fitted_params = []

    for len_win in lengths_window:
        logger.info("start fitting for data from window length={len_win}".format(len_win=len_win))
        df_fitted_params = fit_distribution_for_len_win(len_win, strategies, create_plot=False)
        list_df_fitted_params.append(df_fitted_params)

    df_fitted_params_all_lengths = pd.concat(list_df_fitted_params)

    df_fitted_params_all_lengths.sort_values(by=['strategy', 'len_win'], inplace=True)

    df_strategies["k_on_d"] = df_strategies["k_on"] / df_strategies["k_d"]
    df_strategies["k_off_d"] = df_strategies["k_off"] / df_strategies["k_d"]
    df_strategies["k_syn_d"] = df_strategies["k_syn"] / df_strategies["k_d"]

    df_fitted_params_all_lengths = pd.merge(df_fitted_params_all_lengths, df_strategies, how="left",
                                            left_on=['strategy'],
                                            right_on=['name'])

    df_fitted_params_all_lengths.drop(columns=['name', 'coord_group', 'tm_id', 'fraction_ON' ], inplace=True)

    return df_fitted_params_all_lengths


def prepare_data(df):
    # Calculate k_syn
    # First determine active fraction
    df["active_fraction"] = 1 - df["zero_fraction"]

    # mean_rna = k_syn * active_fraction * len_win
    # => k_syn = mean_rna / (active_fraction * len_win)
    df["k_syn_fit"] = df["real_mean"] / \
                      (df["active_fraction"] * df["len_win"])

    # make categories for parameter ranges
    df["k_on_cat"] = df.k_on.apply(num_categorize)
    df["k_off_cat"] = df.k_off.apply(num_categorize)
    df["k_d_cat"] = df.k_d.apply(num_categorize)
    df["k_syn_cat"] = df.k_syn.apply(num_categorize)

    return df


def num_categorize(number):

    if number < 1e-5:
        ret_val = "< 1e-5"
    elif number < 1e-4:
        ret_val = "< 1e-4"
    elif number < 1e-3:
        ret_val = "< 1e-3"
    elif number < 1e-2:
        ret_val = "< 1e-2"
    elif number < 1e-1:
        ret_val = "< 1e-1"
    else:
        ret_val = ">= 1e-1"
    return ret_val


def lmplot_for(df, x_measure, y_measure, hue_category, len_win):
    hue = "k_on_cat"
    hue_order = ["< 1e-5", "< 1e-4", "< 1e-3", "< 1e-2", "< 1e-1", ">= 1e-1"]
    sns.set(rc={'figure.figsize': (12, 5)})
    sns.lmplot(x=x_measure, y=y_measure, data=df, fit_reg=True,
               hue=hue_category, hue_order=hue_order, legend=False)

    plt.title("{} against {} for t={}".format(y_measure, x_measure, len_win))
    ident = [0, int(max(df[x_measure]))]
    plt.plot(ident, ident, linestyle=':')

    plt.legend(title=hue)

    plt.show()
    plt.close(1)


window_lengths = [15, 30, 45, 60, 75, 90, 105, 120]

run_fitting = False
# strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

csv_name = plot_dir + dir_sep + "parameter_fits.csv"
if run_fitting:
    # replace with 100 strategies
    # strategies = ["generated_" + str(i) for i in range(1, 21)]
    strategies = ["generated_" + str(i) for i in range(1, 101)]
    df = infer_parameters(window_lengths, strategies)
    df.to_csv(csv_name, sep=";", index=False)
else:
    df = pd.read_csv(csv_name, sep=';')


logging.info("Average error estimated poisson mean: {}".
             format(np.average(df.error_estimated_poisson_mean)))

logging.info("Average error real mean: {}".
             format(np.average(df.error_real_mean)))

df = prepare_data(df)
len_win = 15
df = df[df.len_win == len_win]

# lmplot_for(df, "fraction_OFF", "zero_fraction", "k_on_cat", len_win)
# lmplot_for(df, "calculated_mean", "real_mean", "k_on_cat", len_win)
lmplot_for(df, "k_syn", "k_syn_fit", "k_on_cat", len_win)

