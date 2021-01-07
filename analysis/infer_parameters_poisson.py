# here we fit a Poisson distribution to dummy data and try to find a global minimum
# by changing inital conditions

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

strategy = "third_example"
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


# simplified model
def poisson_model(x, mu):

    return poisson.pmf(x, mu)


def fit_to_poisson_and_display(mu):

    expected_guesses = [mu/4]

    for expected in expected_guesses:

        popt, pcov = curve_fit(poisson_model, x, y_norm, expected)

        estimated_mean = round_sig(popt[0], 3)

        plt.title("initial guess={expected}; estimated mean={mean}".format(expected=expected, mean=estimated_mean))

        plt.plot(x, poisson.pmf(x, estimated_mean), label="Poisson")
        plt.plot(x, y_norm, label="real")
        plt.legend()
        plt.show()
        plt.close(1)


def find_fitting_optimum_by_changing_initial_guesses(mu):

    expected_guesses = np.linspace(0, mu * 1.5, 10)
    cov = []
    means = []
    for expected in expected_guesses:
        popt, pcov = curve_fit(poisson_model, x, y_norm, expected)
        estimated_mean = round_sig(popt[0], 3)
        cov.append(pcov[0])
        means.append(estimated_mean)

    plt.title("find fitting optimum by changing initial guesses")
    plt.plot(expected_guesses, cov)
    plt.xlabel("initial guess")
    plt.ylabel("variance")
    plt.show()
    plt.close(1)


sim_dis = SimulatedDistribution(df_counts, nr_cells, strategy)
df_distribution, real_mean = sim_dis.create(label_2)

x = df_distribution["real_count"]
y_norm = df_distribution["chance"]

fit_to_poisson_and_display(real_mean)

# find_fitting_optimum_by_changing_initial_guesses(real_mean)
