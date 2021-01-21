# feature analysis/descriptive statistics of transcript count tables
# for one experiment with one window length and two windows
#
# we would like to see correlations between sets of dynamic parameters ("strategies") and some features
# derived from the counts of labeled and unlabeled transcripts generated by simulations
#
# we read the counts from the df_counts_W[window_length]_G[gap_length].csv generated by main
# e.g. df_counts_W60_G0.csv
# and then join with the strategies file which has been used in main
#   e.g. the generated strategies strategies_generated.csv (use df_strategies of StrategyReader)
#   or the strategies.csv with some predefined strategies
#
# then for each allele we calculate a number of measures/features
# e.g. the fraction of cells for which a count is detected (or not),
# the mean of the distribution, etc, the number of counts
# make two dimensional scatter and density plot for two labels for one strategy (set of dynamical parameters)
# examine correlation between k_d and inferred k_d from means from label_1 and label_2 (assumes non-changing parameters)
# extra analysis of chance of being ON in 2nd window give chance ON or OFF in 1st window
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
plot_dir = out_dir + dir_sep + "analyze_parameters.plots"
os.makedirs(plot_dir, exist_ok=True)

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

efficiency = 100

gap = 0
label_1 = "EU"
label_2 = "4SU"

sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
# sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


def plot_scatter_k_on_k_off():
    min_value = min(np.log10(df_strategies.k_off))
    max_value = max(np.log10(df_strategies.k_off))
    ident = [min_value, max_value]
    plt.plot(ident, ident, linestyle=':')
    plt.scatter(np.log10(df_strategies.k_off), np.log10(df_strategies.k_on))
    plt.title("k_on vs k_off: inactive periods are longer than active periods")
    plt.xlabel("log10(k_off)")
    plt.ylabel("log10(k_on)")

    plot_name = plot_dir + dir_sep + "k_on_vs_k_off.svg"
    plt.savefig(plot_name)
    plt.close(1)


# Correlation between k_syn and mean count of label
def plot_mean_vs_k_syn(label_2, df_means, len_win):

    plt.scatter(df_means.k_syn, df_means.real_count)

    plt.title("Correlation between k_syn and mean count of label {label}; efficiency={eff}%".
              format(label=label_2, eff=efficiency))
    plt.xlabel("transcription rate")
    plt.ylabel("count transcripts sampled with {eff}% efficiency".format(eff=efficiency))

    plot_name = plot_dir + dir_sep + "mean_counts_vs_k_syn_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)
    plt.savefig(plot_name)
    plt.close(1)


def plot_nr_cells_vs_mean_fraction_of_active_state(label_2, df_cell_counts, len_win):
    # plots for number of cells vs mean fraction of time in active state
    ident = [0.0, 0.5]
    plt.plot(ident, ident, linestyle=":")

    plt.scatter(df_cell_counts.fraction_ON, df_cell_counts.real_count/nr_cells)

    plt.title("{label} counts (2nd window); efficiency={eff}%; time={len_win}".
              format(label=label_2, eff=efficiency, len_win=len_win))
    plt.xlabel("fraction of active time")
    plt.ylabel("fraction cells with counts ({eff}% efficiency)".format(eff=efficiency))
    plt.xlim(0, None)
    plt.ylim(0, None)

    plot_name = plot_dir + dir_sep + "nr_cells_vs_active_time_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)
    plt.savefig(plot_name)
    plt.close(1)


# calculate mean expression level per allele (=strategy)
def mean_expression_level(df_counts_unstack):

    df_allele_counts = df_counts_unstack[['allele_id', 'strategy', label_1, label_2]].\
        groupby(['allele_id', 'strategy']).mean().reset_index()

    df_allele_counts = pd.merge(df_allele_counts, df_strategies, how="left",
                                left_on=['strategy'],
                                right_on=['name'])

    len_decay = len_win + gap  # length window + gap
    # pseudo code: k_d = (np.log(<y_w2>) - np.log(<x_w1>)) / t
    df_allele_counts['k_d_predicted'] = \
        (np.log(df_allele_counts[label_2]) - np.log(df_allele_counts[label_1]))/len_decay

    df_allele_counts['k_d_error'] = df_allele_counts['k_d_predicted'] / df_allele_counts['k_d']

    return df_allele_counts


def plot_predicted_k_d(df_allele_counts, len_win):
    ident = [0.0, 0.03]
    plt.figure(figsize=(12, 6))
    plt.plot(ident, ident, linestyle=":")

    plt.scatter(df_allele_counts['k_d'], df_allele_counts['k_d_predicted'])
    plt.xlabel("real k_d")
    plt.ylabel("predicted k_d")

    plt.title("Predicted vs real k_d ({nr_cells} cells/{nr_alleles} alleles); eff={eff}%; time={len_win}".
              format(eff=efficiency, nr_cells=nr_cells, nr_alleles=nr_alleles, len_win=len_win))
    plt.legend(["diagonal (not a regression line)", "one allele"])

    plot_name = plot_dir + dir_sep + "prediction_k_d_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)

    plt.savefig(plot_name)
    plt.close(1)


def plot_error_k_d(df_allele_counts, len_win):

    plt.figure(figsize=(12, 6))
    plt.plot([0, 0.03], [1, 1], linestyle=":")

    plt.scatter(df_allele_counts['k_d'], df_allele_counts['k_d_error'])
    plt.xlabel("real k_d")
    plt.ylabel("relative error k_d")
    plt.title("Relative error of k_d predictions ({nr_cells} cells/{nr_alleles} alleles); eff={eff}%;time={len_win}".
              format(eff=efficiency, len_win=len_win, nr_alleles=nr_alleles, nr_cells=nr_cells))

    plt.legend(["no error (not a regression line)", "allele"])

    plot_name = plot_dir + dir_sep + "relative_error_k_d_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)

    plt.savefig(plot_name)
    plt.close(1)


def create_df_counts_unstack():
    # we can add strategy to the "index" because it is completely defined by key allele_id
    df_counts_unstack = df_counts.set_index(["cell_id", "allele_id", "strategy", "label"])['real_count'].unstack()
    df_counts_unstack = df_counts_unstack.reset_index().fillna(0)

    # what zeroes to include?
    # we want all the counts per cell and allele for label_1
    # where label_2 does not have a zero

    pseudocount = 0.1
    df_counts_unstack["log10_4SU"] = np.log10(df_counts_unstack[label_2] + pseudocount)
    df_counts_unstack["log10_EU"] = np.log10(df_counts_unstack[label_1] + pseudocount)

    # so we need the mean value of the EU label and the mean value of the 4SU label for every allele
    # however, we may underestimate the 1st label, because there will be no counts for the 0 counts
    # so we should divide the total of the 1st label by the number of counts of the 2nd label?
    # we tried that with the following line, but prediction is actually worse
    # df_counts_unstack = df_counts_unstack[df_counts_unstack[label_2] != 0]

    return df_counts_unstack


def joint_scatter_plot_two_labels(df_counts_unstack, len_win):
    sns.jointplot(x=df_counts_unstack[label_1],
                  y=df_counts_unstack[label_2],
                  kind='scatter', s=50, color='b')
    plot_name = plot_dir + dir_sep + "joint_scatter_plot_labels_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)

    plt.savefig(plot_name)

    plt.close(1)


# make two dimensional density plot for two labels for one strategy (set of dynamical parameters)
def joint_kde_plot_two_labels(df_counts_unstack, len_win):
    sns.set(style="white", color_codes=True)
    sns.jointplot(x=df_counts_unstack[label_1], y=df_counts_unstack[label_2], kind='kde', color="skyblue"
                  , xlim=(0, max(df_counts_unstack[label_1] + 5))
                  , ylim=(0, max(df_counts_unstack[label_2] + 5)))
    plot_name = plot_dir + dir_sep + "joint_kde_plot_labels_{eff}_{len_win}.svg".\
        format(eff=efficiency, len_win=len_win)
    plt.savefig(plot_name)

    plt.close(1)


def do_descriptive_analysis(df_counts):

    # mean expression (here zeroes are excluding because there are no rows with zero counts!)
    df_means = df_counts.groupby(['allele_id', 'strategy', 'label'])['real_count'].mean().reset_index()
    df_means = pd.merge(df_means, df_strategies, how="left",
                        left_on=['strategy'],
                        right_on=['name'])
    df_means = df_means[df_means.label == label_2]

    # cell counts (for plotting against fraction of time in active state)
    df_cell_counts = df_counts.groupby(['allele_id', 'strategy', 'label'])['real_count'].count().reset_index()
    df_cell_counts["no_count"] = nr_cells - df_cell_counts.real_count
    df_cell_counts = pd.merge(df_cell_counts, df_strategies, how="left",
                              left_on=['strategy'],
                              right_on=['name'])
    df_cell_counts = df_cell_counts[df_cell_counts.label == label_2]

    # show distribution of k_on and k_off (for comparison with paper)
    plot_scatter_k_on_k_off()

    # Correlation between k_syn and mean count of label
    plot_mean_vs_k_syn(label_2, df_means, len_win)

    # explore dependency of nr cells with counts vs k_on/(k_on + k_off)
    plot_nr_cells_vs_mean_fraction_of_active_state(label_2, df_cell_counts, len_win)

    df_counts_unstack = create_df_counts_unstack()

    # examine correlation between k_d and inferred k_d from means from label_1 and label_2
    # (assumes non-changing parameters) for all strategies
    df_allele_counts = mean_expression_level(df_counts_unstack)
    plot_predicted_k_d(df_allele_counts, len_win)
    plot_error_k_d(df_allele_counts, len_win)

    # strategy = "bimodal"
    strategy = "generated_23"

    # preparing 2-dim density plot of the two labels for a single strategy
    df_counts_unstack = df_counts_unstack[df_counts_unstack.strategy == strategy]

    # make two dimensional scatter and density plot for two labels for one strategy (set of dynamical parameters)
    joint_scatter_plot_two_labels(df_counts_unstack, len_win)
    joint_kde_plot_two_labels(df_counts_unstack, len_win)


def counts_as_proxy_for_p_01_11_t(df_counts, len_win):

    df_label_1 = df_counts[df_counts.label == label_1]
    df_label_2 = df_counts[df_counts.label == label_2]

    df_1_counts = df_label_1.groupby(["allele_id", "strategy", "cell_id"]).real_count.count().reset_index()
    df_1_counts.rename(columns={'real_count': 'count_1'}, inplace=True)

    df_2_counts = df_label_2.groupby(["allele_id", "cell_id"]).real_count.count().reset_index()
    df_2_counts.rename(columns={'real_count': 'count_2'}, inplace=True)

    df_12_counts = pd.merge(df_1_counts, df_2_counts, how="outer",
                            left_on=['allele_id', 'cell_id'],
                            right_on=['allele_id', 'cell_id']).fillna(0)

    # df_count_all = df_12_counts.groupby('allele_id').\
    #     apply(lambda g: g[['cell_id']].count()).\
    #     reset_index()

    # df_count_p10s = df_12_counts.groupby('allele_id').\
    #     apply(lambda g: g[(g.count_1 > 0) & (g.count_2 == 0)][['cell_id']].count()).\
    #     reset_index()

    df_count_p01s = df_12_counts.groupby('allele_id').\
        apply(lambda g:
                                  pd.Series(data=[g[(g.count_1 == 0)][['cell_id']].count().item(),
                                  g[(g.count_1 == 0) & (g.count_2 > 0)][['cell_id']].count().item()]).
              rename({0: "count_0", 1: "count_01"})
              ).\
        reset_index()

    df_count_p11s = df_12_counts.groupby('allele_id').\
        apply(lambda g:
                                  pd.Series(data=[g[(g.count_1 > 0)][['cell_id']].count().item(),
                                  g[(g.count_1 > 0) & (g.count_2 > 0)][['cell_id']].count().item()]).
              rename({0: "count_1", 1: "count_11"})
              ).\
        reset_index()

    df_count_p01s["len_win"] = len_win
    df_count_p11s["len_win"] = len_win

    df_count_p01s["p"] = df_count_p01s["count_01"]/df_count_p01s["count_0"]
    df_count_p11s["p"] = df_count_p11s["count_11"]/df_count_p11s["count_1"]

    return df_count_p01s, df_count_p11s


def plot_chance_vs_time(df_p11_example, y_label):

    plt.figure(figsize=(12, 6))
    plt.plot(df_p11_example.len_win, df_p11_example.p, '.-')
    plt.title("WRONG? Chance of being in state 1 at t when in state 1 at t=0 for allele {allele_id}".
              format(allele_id=allele_id))
    plt.xlabel("time (minutes)")
    plt.ylabel(y_label)
    plt.xlim(0, None)
    plt.ylim(0, None)
    plt.show()
    plt.close(1)


sns.set(style="white", color_codes=True)

window_lengths = [15, 30, 45, 60, 75, 90, 105, 120]

p_01s = []
p_11s = []

allele_id = 6  # 6

for len_win in window_lengths:

    filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
        len_win=len_win, gap=gap, eff=efficiency)

    # df_counts are counts of two labeling windows (single window length)
    df_counts = pd.read_csv(filename_counts, sep=';')

    nr_cells = len(df_counts.cell_id.unique())
    nr_alleles = len(df_counts.allele_id.unique())

    # main analysis (NB: creates lots of plots!)
    do_descriptive_analysis(df_counts)

    # temp for debugging
    # df_counts = df_counts[df_counts.allele_id == allele_id]

    # extra analysis of chance of being ON in 2nd window give chance ON or OFF in 1st window
    df_count_p01s, df_count_p11s = counts_as_proxy_for_p_01_11_t(df_counts, len_win)

    p_01s.append(df_count_p01s)
    p_11s.append(df_count_p11s)


df_p01s = pd.concat(p_01s)
df_p11s = pd.concat(p_11s)

df_p11_example = df_p11s[df_p11s.allele_id == allele_id]

plot_chance_vs_time(df_p11_example, y_label="inferred chance P_11(t))")

