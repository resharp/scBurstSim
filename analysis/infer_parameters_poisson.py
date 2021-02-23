# try to infer dynamic parameters with approximation for time dependent solution
# based on time dependent counts in one labeling window, and multiple window lengths
# 1. fit Poisson distribution to simulated data of counts with 100% efficiency
# 2. does not work well: fit to stationary distribution to infer k_syn/k_d, k_on/k_d and k_off/k_d

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
    df["k_syn_fit_m"] = df["real_mean"] / \
                        (df["active_fraction"] * df["len_win"])

    df["k_syn_fit_p"] = df["estimated_poisson_mean"] / \
                        (df["active_fraction"] * df["len_win"])

    # make categories for parameter ranges
    df["k_on_cat"] = df.k_on.apply(num_categorize)
    df["k_off_cat"] = df.k_off.apply(num_categorize)
    df["k_d_cat"] = df.k_d.apply(num_categorize)
    df["k_syn_cat"] = df.k_syn.apply(num_categorize)

    # error prediction of fraction_OFF =k_off/(k_off + k_on) by using zero fraction (fraction of cells with no counts)
    df["error_zero_fraction"] = np.abs(df["fraction_OFF"] - df["zero_fraction"])/df["fraction_OFF"]

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
    if x_measure == "fraction_OFF":
        ident = [0, 1]
    plt.plot(ident, ident, linestyle=':')

    plt.legend(title=hue)

    plt.show()
    plt.close(1)


def linear_regression_for(df, window_lengths, x_parameter, y_parameter):
    lrs = []
    logging.info("Evaluating method using exact mean of simulated distribution")

    k_on_cats = df.k_on_cat.unique()
    for k_on_cat in k_on_cats:
        df_k_on = df[df.k_on_cat == k_on_cat]

        for len_win in window_lengths:
            df_t = df_k_on[df_k_on.len_win == len_win]

            slope, intercept, r_value, p_value, std_err = linregress(x=df_t[x_parameter], y=df_t[y_parameter])

            # mean absolute error in prediction of fraction_OFF by using zero fraction
            mean_error_zero = df_t.error_zero_fraction.mean()

            lrs.append([len_win, k_on_cat, slope, intercept, r_value, p_value, std_err, mean_error_zero])
            logging.info("slope = {slope} for length {len_win} and k_on_cat {k_on_cat} using {measure}".format(
                slope=round_sig(slope, 3), len_win=len_win, k_on_cat=k_on_cat, measure=y_parameter))
            # logging.info("R^2 = {r2}% f or prediction of k_syn length {len_win} using mean".format(
            #     r2=round_sig(r_value**2 * 100, 3), len_win=len_win))

    df_lr = pd.DataFrame(data=lrs, columns=('len_win', 'k_on_cat', 'slope', 'intercept', 'r_value', 'p_value', 'std_err',
                                            'mean_error_zero'))

    return df_lr


def heat_map_inferred_parameters(df_lr, parameter, measure):
    df_lr = df_lr[['len_win', 'k_on_cat', measure]]
    # df_lr = df_lr.sort_values(["len_win", "k_on_cat"])
    df_lr = df_lr.set_index(["len_win", "k_on_cat"])

    # convert from multi-index to cross-product table
    df_lr = df_lr.unstack()

    # rename columns, unstack
    df_lr.columns = ["_".join(x) for x in df_lr.columns.ravel()]

    offset_column = 1 + measure.count("_")
    df_lr.columns = ["_".join(x.split("_")[offset_column:]) for x in df_lr.columns]

    df_lr = df_lr.transpose()  # if you would like to switch columns and rows

    sns.set(rc={'figure.figsize': (12, 5)})
    ax = sns.heatmap(df_lr, cmap="seismic_r", annot=True)
    plt.title("values for {measure} for fitting {parameter}".format(measure=measure, parameter=parameter))
    plt.xlabel("Length of labeling window (minutes)")
    plt.ylabel("k_on range")
    plt.show()
    plt.close(1)


def infer_parameters_for_window_lenghts(window_lengths, run_fitting=False):

    csv_name = plot_dir + dir_sep + "parameter_fits.csv"
    logger.info("Infer parameters based on time dependent distributions from {} cells".format(nr_cells))
    logger.info("Results in {}".format(csv_name))

    if run_fitting:
        # replace with 100 strategies
        # strategies = ["generated_" + str(i) for i in range(1, 21)]
        strategies = ["generated_" + str(i) for i in range(1, 101)]

        df = infer_parameters(window_lengths, strategies)

        df.to_csv(csv_name, sep=";", index=False)
    else:
        df = pd.read_csv(csv_name, sep=';')

    return df


def evaluate_infer_corrected_mean(df, window_lengths):

    logging.info("Average error estimated poisson mean: {}".
                 format(np.average(df.error_estimated_poisson_mean)))

    logging.info("Average error real mean: {}".
                 format(np.average(df.error_real_mean)))

    df = prepare_data(df)
    len_win = 15
    df_t = df[df.len_win == len_win]

    # x_parameter = "fraction_OFF"; y_parameter = "zero_fraction"
    # x_parameter = "k_syn"; y_parameter = "k_syn_fit_m"
    x_parameter = "calculated_mean"; y_parameter = "estimated_poisson_mean"
    # y_parameter = "real_mean"

    lmplot_for(df, x_parameter, y_parameter, "k_on_cat", len_win)

    df_lr = linear_regression_for(df, window_lengths, x_parameter=x_parameter, y_parameter=y_parameter)

    # possible measures from linear fitting:
    # columns=('len_win', 'k_on_cat', 'slope', 'intercept', 'r_value', 'p_value', 'std_err', 'mean_error_zero')

    measure = "r_value"

    heat_map_inferred_parameters(df_lr, y_parameter, measure)


# TODO: set run_fitting to True if you use new data
run_fitting = False

window_lengths = [15, 30, 45, 60, 75, 90, 105, 120]

df_predictions = infer_parameters_for_window_lenghts(window_lengths, run_fitting)

evaluate_infer_corrected_mean(df_predictions, window_lengths)
