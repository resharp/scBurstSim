# analysis of signal of label counts
#
# input: all counts for 500 cells for all combinations of:
# efficiencies = [1, 0.5, 0.2, 0.05]
# window_lengths (minutes) = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195]
#
# output: directory signal_labels.plots containing
#   - phase diagrams of mean counts of labels, and number of cells that have non-zero counts
#   - split for categories of synthesis and decay rate (half life)
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


def run_all_counts(prj_dir, nr_cells, efficiencies, window_lengths, df_strategies):
    list_agg = []

    for eff in efficiencies:
        out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)

        for len_win in window_lengths:

            filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
                len_win=len_win, gap=gap)

            # df_counts are counts of two labeling windows (single window length)
            df_counts = pd.read_csv(filename_counts, sep=';')

            df_merged = df_counts.merge(df_strategies, how='inner',
                                        left_on=['strategy'],
                                        right_on=['name'])

            # join strategies, add transcription and decay rate, and calculate mean_count per allele
            df_agg = df_merged.groupby(['allele_id', 'strategy', 'label', 'k_syn', 'k_d']).agg(
                sum_count=('real_count', 'sum'),
                cell_count=('real_count', 'count')
                 ).reset_index()

            df_agg["mean_count"] = df_agg["sum_count"]/nr_cells
            df_agg["eff"] = eff
            df_agg["len_win"] = len_win

            list_agg.append(df_agg)

    df_agg_all = pd.concat(list_agg)

    return df_agg_all


def add_coord_group_to_strategy(df_alleles):

    df_merged = df_alleles.merge(df_strategies, how='inner',
                                 left_on=['strategy'],
                                 right_on=['name'])
    df_merged["display"] = np.where( df_merged.coord_group.isnull(),
                                     df_merged['strategy'],
                                     df_merged['strategy'] + "__cg" + df_merged.coord_group.astype(str))

    return df_merged


def signal_phase_diagram(data, measure, label, k_syn, half_life_h, title):

    data = data[['eff', 'len_win', measure]]
    data = data.sort_values("len_win")
    data = data.set_index(["eff", "len_win"])

    # convert from multi-index to cross-product table
    data = data.unstack()

    # rename columns, unstack
    data.columns = [x[1] for x in data.columns.ravel()]

    plt.figure(figsize=(12, 5))
    ax = sns.heatmap(data, cmap="Spectral_r", annot=True, fmt='g')

    plt.title(title)

    plot_dir = r"{}{}runs{}{}".format(prj_dir, dir_sep, dir_sep, "signal_of_labels.plots")
    phase_plot_name = plot_dir + dir_sep + "phase_plot_{}_{}_{}_{}.svg".format(measure, label, k_syn, half_life_h)

    plt.xlabel("window size in minutes (no gap)")
    plt.ylabel("efficiency")

    plt.savefig(phase_plot_name)
    plt.close(1)


def make_phase_diagrams(df_agg2):
    # split for transcription rates and half lives
    k_syns = df_agg2.k_syn.unique()
    half_lives = df_agg2.half_life_h.unique()

    for k_syn in k_syns:
        for half_life_h in half_lives:
            # temp: only look at label_1, since its abundance is the bottleneck in the analysis
            labels = [label_1]
            label_nr = 0
            for label in labels:
                label_nr = label_nr + 1
                data = df_agg2[
                    (df_agg2.k_syn == k_syn) & (df_agg2.half_life_h == half_life_h) & (df_agg2.label == label)]

                measure = 'mean_count'
                title = "mean transcript abundance per gene for label {}, k_syn: {}, half life: {}". \
                    format(label_nr, k_syn, half_life_h)
                signal_phase_diagram(data, measure, label, k_syn, half_life_h, title)

                measure = 'cell_count'
                title = "mean nr of cells with counts per gene, label {}, k_syn: {}, half life: {}". \
                    format(label_nr, k_syn, half_life_h)
                signal_phase_diagram(data, measure, label, k_syn, half_life_h, title)


plot_dir = prj_dir + dir_sep + "runs" + dir_sep + "signal_of_labels.plots"
os.makedirs(plot_dir, exist_ok=True)

strategies_file = prj_dir + dir_sep + "runs" + dir_sep + "strategies_mixed_new.csv"
sr = StrategyReader(strategies_file)
sr.read_strategies()
df_strategies = sr.df_strategies

efficiencies = [1, 0.5, 0.2, 0.05]
window_lengths = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195]

total_file_name = "{plot_dir}{dir_sep}total_label_means_G{gap}.csv".format(
    plot_dir=plot_dir, dir_sep=dir_sep, gap=gap)
alleles_file_name = "{plot_dir}{dir_sep}allele_label_means_G{gap}.csv".format(
    plot_dir=plot_dir, dir_sep=dir_sep, gap=gap)

run_agg = False
if run_agg:

    nr_cells = 500
    # aggregate over cells
    df_agg = run_all_counts(prj_dir, nr_cells, efficiencies, window_lengths, df_strategies)

    df_agg["period"] = df_agg.strategy.str.split("_", expand=True)[4]
    # df_agg = add_coord_group_to_strategy(df_agg)
    df_agg.to_csv(path_or_buf=alleles_file_name, sep=';', index=False)

    # aggregate over all alleles of a certain synthesis rate
    df_agg2 = df_agg.groupby(['eff', 'len_win', 'label', 'k_syn', 'k_d']).agg(
        mean_count=('mean_count', 'mean'),
        cell_count=('cell_count', 'mean')
    ).reset_index()

    df_agg2.mean_count = df_agg2.mean_count.round(2)
    df_agg2.cell_count = df_agg2.cell_count.round(0)

    df_agg2 = df_agg.groupby(['eff', 'len_win', 'label', 'k_syn', 'k_d']).agg(
        mean_count=('mean_count', 'mean'),
        cell_count=('cell_count', 'mean')
    ).reset_index()

    df_agg2["half_life_h"] = (np.log(2) / (60 * df_agg2["k_d"])).round(1).map(str)

    df_agg2.mean_count = df_agg2.mean_count.round(2)
    df_agg2.cell_count = df_agg2.cell_count.round(0)

    df_agg2.to_csv(path_or_buf=total_file_name, sep=';', index=False)
else:
    df_agg = pd.read_csv(alleles_file_name, sep=';')
    df_agg2 = pd.read_csv(total_file_name, sep=';')

make_phase_diagrams(df_agg2)

# df_agg2_eff = df_agg2[(df_agg2.eff == 0.2) & (df_agg2.label == label_1)]
# for each half life and k_syn, calculate pulse length giving maximum counts
# https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-value-in-groups-using-groupby

debug = True

