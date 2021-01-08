# here we fit a Poisson distribution to simulated data of "real" transcripts
# based on time dependent counts in one labeling window
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
plot_dir = out_dir + dir_sep + "show_distributions.plots"
os.makedirs(plot_dir, exist_ok=True)

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

strategy = "second_example"
nr_cells = 3000

len_win = 60
gap = 0

label_2 = "4SU"

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(len_win=len_win, gap=gap)

df_counts = pd.read_csv(filename_counts, sep=';')


def poisson_model(x, mu):

    return poisson.pmf(x, mu)


# naive method with one inital guess
def fit_to_poisson_and_display(strategy, x, y_norm, mu):

    expected_guesses = [mu/4]

    for expected in expected_guesses:

        popt, pcov = curve_fit(poisson.pmf, x, y_norm, p0=expected, maxfev=5000)

        estimated_mean = round_sig(popt[0], 3)
        variance = round_sig(pcov[0, 0], 3)

        plt.title("{strategy};initial guess={expected}; estimated mean={mean};variance={variance}".format(
            strategy=strategy, expected=expected, mean=estimated_mean, variance=variance))

        plt.plot(x, poisson.pmf(x, estimated_mean), label="Poisson")
        plt.plot(x, y_norm, label="real")
        plt.legend()
        plt.show()
        plt.close(1)


def find_fitting_optimum_by_changing_initial_guesses(strategy, x, y_norm, mu):

    # we do not start at 0 or 1 because somehow it gives a very small mean
    expected_guesses = np.linspace(2, mu * 1.5, 5)
    cov = []
    means = []
    for expected in expected_guesses:
        popt, pcov = curve_fit(poisson.pmf, x, y_norm, p0=expected, maxfev=1400)

        estimated_mean = round_sig(popt[0], 3)

        cov.append(pcov[0, 0])
        means.append(estimated_mean)

    # find minimum variance for guesses
    var, idx = min((var, idx) for (idx, var) in enumerate(cov))

    estimated_mean = round_sig(means[idx], 3)   # and related mean
    var = round_sig(var, 3)
    plt.title("{strategy};optimal fitting;variance={variance}".
              format(strategy=strategy, mean=estimated_mean, variance=var, real_mean=mu))
    plt.plot(x, poisson.pmf(x, estimated_mean), label="Poisson, mean={}".format(estimated_mean))
    plt.plot(x, y_norm, label="real, mean={}".format(mu))
    plt.legend()
    plt.show()
    plt.close(1)


strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

sim_dis = SimulatedDistribution(df_counts, nr_cells)

for strategy in strategies:

    df_distribution, real_mean = sim_dis.create(strategy, label_2)

    # exclude first row from fitting. The first row includes the zero counts
    df_distribution = df_distribution.iloc[1:]

    x = df_distribution["real_count"]
    y_norm = df_distribution["chance"]

    # fit_to_poisson_and_display(strategy, x, y_norm, real_mean)

    find_fitting_optimum_by_changing_initial_guesses(strategy, x, y_norm, real_mean)

