# plotting (properties of) distributions
# create distributions from simulations
#
# we read the counts from the df_counts_W[window_length]_G[gap_length].csv generated by main
# e.g. df_counts_W60_G0.csv
# and then join with the strategies file which has been used in main
#   e.g. the generated strategies strategies_generated.csv (use df_strategies of StrategyReader)
#   or the strategies.csv with some predefined strategies
# compare simulated (time-dependent and stationary) distributions against theoretical stationary distribution
# examine how quickly the means converge towards the means of the stationary distributions (use multiple window lengths)
import os

import matplotlib.pyplot as plt
import pandas as pd

from analysis.SimulatedDistribution import *
from simulator.StrategyReader import StrategyReader
from solution.stationary import create_distribution

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

# strategy = "bimodal"
strategy = "second_example"
nr_cells = 3000
efficiency = 100

len_win = 60  # 60, 120, 180, 240
gap = 0
label_1 = "EU"
label_2 = "4SU"

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap, eff=efficiency)

df_counts = pd.read_csv(filename_counts, sep=';')


def plot_label_distribution(label, df_counts, strategy, nr_cells):

    return plot_distribution("real_count", df_counts, strategy, nr_cells, label=label)


def plot_stationary_distribution(df_counts, strategy, nr_cells):

    return plot_distribution("count_all", df_counts, strategy, nr_cells)


# goal: compare simulated distribution against theoretical stationary distribution
# for both labeled transcripts (non stationary) and total number of transcripts (should approach stationary)
# measure may be count_all or real_count (for label)
def plot_distribution(measure, df_counts, strategy, nr_cells, label=None):

    params = sr.get(strategy=strategy)

    sim_dis = SimulatedDistribution(df_counts, nr_cells, strategy)

    df_distribution, real_mean = sim_dis.create(measure, label)

    # distribution from simulation
    plt.step(df_distribution[measure], df_distribution.chance, where="post")

    # theoretical distribution
    x_list, y_list = create_distribution(params.k_on, params.k_off, params.k_syn, params.k_d)
    plt.step(x_list, y_list, where="post", color="red")

    plt.xlabel("Number of molecules")
    if measure == "count_all":
        title = "Distribution all mRNA for strategy '{strategy}'; mean={real_mean}".\
            format(strategy=strategy, real_mean=real_mean)
        plt.ylabel("Chance P(N)")
    else:
        title = "Distribution labeled mRNA for strategy '{strategy}'; mean={real_mean}".\
            format(strategy=strategy, real_mean=real_mean)
        plt.ylim(0, 0.08)
        plt.ylabel("Chance P(N,t)")
    plt.title(title)

    if measure == "count_all":
        label = "stationary"

    plot_name = plot_dir + dir_sep + "{label}_distribution_{strategy}_{nr_cells}_{len_win}.svg".format(
        label=label, strategy=strategy, nr_cells=nr_cells, len_win=len_win)

    plt.savefig(plot_name)
    plt.close(1)


# compare simulated (time-dependent and stationary) distributions against theoretical stationary distribution
def plot_distributions(nr_cells):

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

    for strategy in strategies:

        plot_stationary_distribution(df_counts, strategy, nr_cells)
        plot_label_distribution(label_2, df_counts, strategy, nr_cells)


# calculate the mean (TODO: variance) of the number of molecules
# and see how well it predicts the dynamic parameters
# x-axis: time
# y-axis: mean
# df_counts is determined by the right time
def plot_means_against_time(label_1, label_2):
    times = [60, 120, 180, 240]

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

    for strategy in strategies:
        real_means_1 = []
        real_means_2 = []

        for len_win in times:

            filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
                len_win=len_win, gap=gap, eff=efficiency)

            df_counts = pd.read_csv(filename_counts, sep=';')

            sim_dis = SimulatedDistribution(df_counts, nr_cells, strategy)

            df_distribution, real_mean = sim_dis.create("real_count", label_1)
            real_means_1.append(real_mean)
            df_distribution, real_mean = sim_dis.create("real_count", label_2)
            real_means_2.append(real_mean)
            # print("real mean for t={len_win}: {real_mean}".format(len_win=len_win, real_mean=real_mean))

        # stationary distribution
        df_distribution, stat_mean = sim_dis.create("count_all")

        plt.plot(times, real_means_1, 'o-', label="label 1")
        plt.plot(times, real_means_2, 'o-', label="label 2")
        plt.title("time dependent mean for strategy: " + strategy + ";stationary mean={stat_mean}".format(
            stat_mean=stat_mean
        ))
        plt.xlim(0, None)
        plt.xlabel("time in minutes")
        plt.ylabel("mean nr of molecules")
        plt.legend()

        plot_name = plot_dir + dir_sep + "mean_vs_time_{strategy}_{nr_cells}.svg".format(
            strategy=strategy, nr_cells=nr_cells)
        plt.savefig(plot_name)
        plt.close(1)


# compare simulated (time-dependent and stationary) distributions against theoretical stationary distribution
plot_distributions(nr_cells)

# examine how quickly the means converge towards the means of the stationary distributions (use multiple window lengths)
plot_means_against_time(label_1, label_2)
