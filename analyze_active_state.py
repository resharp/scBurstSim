import os
import pandas as pd
import matplotlib.pyplot as plt
from simulator.StrategyReader import StrategyReader
import numpy as np
import seaborn as sns

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "analyze_active_state.plots"
os.makedirs(plot_dir, exist_ok=True)

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

len_win = 180  # 60, 100, 180
gap = 0
label_1 = "EU"
label_2 = "4SU"

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap)
df_counts = pd.read_csv(filename_counts, sep=';')


def get_parameters(strategy):
    df_strategy = df_strategies[df_strategies.name == strategy]

    return df_strategy.k_on.item(), df_strategy.k_off.item()


def plot_distribution_active_state():

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]
    # strategies = ["powerlaw"]

    for strategy in strategies:

        df_counts_label = df_counts[(df_counts.label == label_2) & (df_counts.strategy == strategy)]

        # convert to density for perc (=percentage active state)

        bins = np.arange(0, 1.01, 0.01)

        pdf = np.histogram(df_counts_label.perc, bins=bins)[0]

        bins = np.delete(bins, 0)  # remove first
        sum_pdf = sum(pdf)
        norm_pdf = pdf/sum_pdf

        plt.step(bins, norm_pdf, where="post")
        plt.ylim(0, None)
        plt.xlabel("percentage active state")
        plt.ylabel("chance")

        k_on, k_off = get_parameters(strategy)

        plt.title("fraction active period ({strategy});k_on={k_on},k_off={k_off} (win={len_win})".format(
            strategy=strategy, k_on=k_on, k_off=k_off, len_win=len_win))
        plot_name = plot_dir + dir_sep + "percentage_state_distribution_{strategy}_{len_win}.svg".\
            format(strategy=strategy, len_win=len_win)
        plt.savefig(plot_name)
        plt.close(1)


# create distributions for percentage of active state ("Hidden Markov distribution")
plot_distribution_active_state()


