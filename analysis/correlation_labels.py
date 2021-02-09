from analysis.data_analysis import *
# we use correlation between two label counts as a filter

import os
import pandas as pd
from simulator.StrategyReader import StrategyReader
from analysis.data_analysis import *

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

plot_dir = out_dir + dir_sep + "correlation_labels.plots" + dir_sep + "len_win_{}".format(len_win)
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


df_counts = normalize_counts(df_counts)
df_counts_12 = merge_label_counts(df_counts, label_1, label_2)

# calculate Pearson's correlation for every allele like in
# https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/
