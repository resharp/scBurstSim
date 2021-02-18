# analysis of signal of label counts
# output: directory signal_labels.plots containing
#   - phase diagrams of mean counts of labels
import os
import pandas as pd
import seaborn as sns
from simulator.StrategyReader import StrategyReader
from analysis.data_analysis import *
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2

gap = 0

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    # out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
    prj_dir = r"D:\26 Battich Oudenaarden transcriptional bursts"
else:
    dir_sep = "/"
    prj_dir = "sc_runs"

gap = 0
label_1 = "EU"
label_2 = "4SU"


def run_all_counts(prj_dir, efficiencies, window_lengths):
    list_agg = []

    for eff in efficiencies:
        out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)

        for len_win in window_lengths:

            filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
                len_win=len_win, gap=gap)

            # df_counts are counts of two labeling windows (single window length)
            df_counts = pd.read_csv(filename_counts, sep=';')

            df_agg = df_counts.groupby(['allele_id', 'strategy', 'label']).agg(
                mean_count=('real_count', 'mean'),
                cell_count=('real_count', 'count')
                 ).reset_index()

            df_agg["eff"] = eff
            df_agg["len_win"] = len_win

            list_agg.append(df_agg)

    # df_counts = normalize_counts(df_counts)
    # df_counts_12 = merge_label_counts(df_counts, label_1, label_2)

    df_agg_all = pd.concat(list_agg)

    return df_agg_all


plot_dir = prj_dir + dir_sep + "runs" + dir_sep + "signal_of_labels.plots"
os.makedirs(plot_dir, exist_ok=True)


efficiencies = [1, 0.5, 0.2, 0.05]
# efficiencies = [1]
window_lengths = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195]
# window_lengths = [60, 90]

agg_file_name = "{plot_dir}{dir_sep}label_means_G{gap}.csv".format(
    plot_dir=plot_dir, dir_sep=dir_sep, gap=gap)

run_agg = False
if run_agg:

    # aggregate over cells
    df_agg = run_all_counts(prj_dir, efficiencies, window_lengths)

    # aggregate over all alleles
    df_agg2 = df_agg.groupby(['eff', 'len_win', 'label']).agg(
        mean_count=('mean_count', 'mean'),
        cell_count=('cell_count', 'mean')
    ).round(2).reset_index()

    df_agg2.to_csv(path_or_buf=agg_file_name, sep=';', index=False)
else:
    df_agg2 = pd.read_csv(agg_file_name, sep=';')

debug = True





