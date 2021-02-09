# we use correlation between two label counts as a filter

import os
import pandas as pd
from simulator.StrategyReader import StrategyReader
from analysis.data_analysis import *
from scipy.stats import pearsonr

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
    # out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs_on_server"
else:
    dir_sep = "/"
    out_dir = "sc_runs"

efficiency = 100

gap = 0
label_1 = "EU"
label_2 = "4SU"
len_win = 60

plot_dir = out_dir + dir_sep + "correlation_labels.plots"
os.makedirs(plot_dir, exist_ok=True)

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap, eff=efficiency)


# df_counts are counts of two labeling windows (single window length)
df_counts = pd.read_csv(filename_counts, sep=';')

strategies_file = out_dir + dir_sep + "strategies_mixed.csv"
sr = StrategyReader(strategies_file)
# sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


def scatter_plot_for_allele(df_counts, strategy):
    df_counts_allele = df_counts[df_counts_12.strategy == strategy]

    corr, _ = pearsonr(df_counts_allele.norm_count_1, df_counts_allele.norm_count_2)

    corr_s = '%.3f' % corr
    # print("Pearsons correlation: {}".format(corr_s))

    plt.title("strategy: {} with correlation {} (len_win={})".format(strategy, corr_s, len_win))
    plt.scatter(df_counts_allele.norm_count_1, df_counts_allele.norm_count_2)
    plt.xlabel("normalized counts label 1")
    plt.ylabel("normalized counts label 2")

    fig_name = plot_dir + dir_sep + "scatter_labels_{strategy}_{len_win}.svg".format(strategy=strategy, len_win=len_win)
    plt.savefig(fig_name)
    plt.close(1)


def label_pearson_correlation(series):

    corr, _ = pearsonr(series.norm_count_1, series.norm_count_2)

    df = pd.DataFrame([corr], columns=["corr"])

    return df


def calculate_corr_and_save(df_counts):

    df_corr = df_counts.groupby(["strategy"]).apply(label_pearson_correlation).reset_index(). \
        sort_values("corr").drop('level_1', axis=1)

    filename_corr_labels = plot_dir + dir_sep + "corr_labels_{len_win}.csv".format(len_win=len_win)
    df_corr.to_csv(filename_corr_labels, sep=";", index=False)
    return df_corr


df_counts = normalize_counts(df_counts)
df_counts_12 = merge_label_counts(df_counts, label_1, label_2)


df_corr = calculate_corr_and_save(df_counts_12)


strategy = df_corr.head(1).strategy.item()
scatter_plot_for_allele(df_counts_12, strategy)

strategy = df_corr.tail(1).strategy.item()
scatter_plot_for_allele(df_counts_12, strategy)
