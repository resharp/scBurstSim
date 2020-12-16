import os
import pandas as pd
import matplotlib.pyplot as plt
from simulator.StrategyReader import StrategyReader
import numpy as np
import seaborn as sns
# we would like to analyse how well we can infer the strategy parameters from
# the counts
# we read the counts from the df_counts.csv generated by main
# and then join with the strategies file which has been used in main
#   e.g. the generated strategies strategies_generated.csv (use df_strategies of StrategyReader)
#   or the strategies.csv
# then for each allele we calculate the fraction of cells for which a count is detected
# (and not detected)
from solution.stationary import create_distribution
from utils.utils import round_sig

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

# strategy = "bimodal"
strategy = "second_example"
nr_cells = 3000
efficiency = 100
label = "4SU"  # 2nd window
len_win = 60  # 60, 120, 180, 240
gap = 0


# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


def plot_scatter_k_on_k_off():
    max_value = max(np.log10(df_strategies.k_off))
    ident = [0, max_value]
    plt.plot(ident, ident)
    plt.scatter(np.log10(df_strategies.k_off), np.log10(df_strategies.k_on))
    plt.xlabel("log10(k_off)")
    plt.ylabel("log10(k_on)")
    plt.show()


# plot_scatter_k_on_k_off()


filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap, eff=efficiency)

df_counts = pd.read_csv(filename_counts, sep=';')

# mean expression #TODO: Check if this is the correct way to calculate mean? What about zeroes?
df_means = df_counts.groupby(['allele_id', 'strategy', 'label'])['real_count'].mean().reset_index()
df_means = pd.merge(df_means, df_strategies, how="left",
                    left_on=['strategy'],
                    right_on=['name'])
df_means = df_means[df_means.label == label]

# cell counts (for plotting against fraction of time in active state)
df_cell_counts = df_counts.groupby(['allele_id', 'strategy', 'label'])['real_count'].count().reset_index()
df_cell_counts["no_count"] = nr_cells - df_cell_counts.real_count
df_cell_counts = pd.merge(df_cell_counts, df_strategies, how="left",
                          left_on=['strategy'],
                          right_on=['name'])
df_cell_counts = df_cell_counts[df_cell_counts.label == label]


# plot for mean counts vs k_syn
def plot_mean_vs_k_syn():

    plt.scatter(df_means.k_syn, df_means.real_count)

    plt.title("Correlation between k_syn and mean count of label {label}; efficiency={eff}%".
              format(label=label, eff=efficiency))
    plt.xlabel("transcription rate")
    plt.ylabel("count transcripts sampled with {eff}% efficiency".format(eff=efficiency))

    plot_name = plot_dir + dir_sep + "mean_counts_vs_k_syn_{eff}.svg".format(eff=efficiency)
    plt.savefig(plot_name)
    plt.close(1)

    # plots for number of cells vs fraction of time in active state
    plt.scatter(df_cell_counts.fraction_ON, df_cell_counts.real_count)

    plt.title("{label} counts (2nd window); efficiency={eff}%".
              format(label=label, eff=efficiency))
    plt.xlabel("fraction of active time")
    plt.ylabel("# cells with counts ({eff}% efficiency)".format(eff=efficiency))

    plot_name = plot_dir + dir_sep + "nr_cells_vs_active_time_{eff}.svg".format(eff=efficiency)
    plt.savefig(plot_name)
    plt.close(1)


plot_mean_vs_k_syn()

# example
t = len_win + gap # length window + gap

y_w2 = 100
x_w1 = 10
k_d = (np.log(y_w2) - np.log(x_w1)) / t

# so we need the mean value of the EU label and the mean value of the 4SU label for every allele
# however, we may underestimate the 1st label, because there will be no counts for the 0 counts
# so we should divide the total of the 1st label by the number of counts of the 2nd label?
# we tried that with the following line, but prediction is actually worse
# df_counts_unstack = df_counts_unstack[df_counts_unstack[label_2] != 0]


# calculate mean expression level
def mean_expression_level():

    df_allele_counts = df_counts_unstack[['allele_id', 'strategy', label_1, label_2]].\
        groupby(['allele_id', 'strategy']).mean().reset_index()

    df_allele_counts = pd.merge(df_allele_counts, df_strategies, how="left",
                                left_on=['strategy'],
                                right_on=['name'])

    k_d = (np.log(y_w2) - np.log(x_w1)) / t

    df_allele_counts['k_d_predicted'] = (np.log(df_allele_counts[label_2]) - np.log(df_allele_counts[label_1]))/t

    df_allele_counts['k_d_error'] = df_allele_counts['k_d_predicted'] / df_allele_counts['k_d']

    return df_allele_counts


# df_allele_counts = mean_expression_level()


def plot_error_k_d():

    plt.plot([0, 0.03], [1, 1])

    plt.scatter(df_allele_counts['k_d'], df_allele_counts['k_d_error'])
    plt.xlabel("real k_d")
    plt.ylabel("relative error k_d")
    plt.title("Relative error of k_d predictions (100 cells/100 alleles); efficiency={eff}%".format(eff=efficiency))

    plt.legend(["no error (not a regression line)", "allele"])

    plot_name = plot_dir + dir_sep + "relative_error_k_d_{eff}.svg".format(eff=efficiency)

    plt.savefig(plot_name)
    plt.close(1)


def plot_predicted_k_d():
    ident = [0.0, 0.03]
    plt.plot(ident, ident)

    plt.scatter(df_allele_counts['k_d'], df_allele_counts['k_d_predicted'])
    plt.xlabel("real k_d")
    plt.ylabel("predicted k_d")

    plt.title("Predicted vs real k_d (100 cells/100 alleles); efficiency={eff}%".format(eff=efficiency))
    plt.legend(["diagonal (not a regression line)", "allele"])

    plot_name = plot_dir + dir_sep + "prediction_k_d_{eff}.svg".format(eff=efficiency)

    plt.savefig(plot_name)
    plt.close(1)


# plot_predicted_k_d()
# plot_error_k_d()


def create_df_counts_unstack():
    # we can add strategy to the "index" because it is completely defined by key allele_id
    df_counts_unstack = df_counts.set_index(["cell_id", "allele_id", "strategy", "label"])['real_count'].unstack()
    df_counts_unstack = df_counts_unstack.reset_index().fillna(0)

    # what zeroes to include?
    # we want all the counts per cell and allele for label_1
    # where label_2 does not have a zero
    label_1 = "EU";
    label_2 = "4SU"

    # we would like to make a 2-dim density plot of the two labels for a chosen strategy

    df_counts_unstack = df_counts_unstack[df_counts_unstack.strategy == strategy]

    pseudocount = 0.1
    df_counts_unstack["log10_4SU"] = np.log10(df_counts_unstack["4SU"] + pseudocount)
    df_counts_unstack["log10_EU"] = np.log10(df_counts_unstack["EU"] + pseudocount)

    return df_counts_unstack


df_counts_unstack = create_df_counts_unstack()


def joint_scatter_plot_labels():
    sns.jointplot(x=df_counts_unstack["EU"],
                  y=df_counts_unstack["4SU"],
                  kind='scatter', s=50, color='b')
    plot_name = plot_dir + dir_sep + "joint_scatter_plot_labels_{eff}.svg".format(eff=efficiency)

    plt.savefig(plot_name)

    plt.close(1)


def joint_kde_plot_labels():
    sns.set(style="white", color_codes=True)
    sns.jointplot(x=df_counts_unstack["EU"], y=df_counts_unstack["4SU"], kind='kde', color="skyblue"
                  , xlim=(0, max(df_counts_unstack["EU"] + 5))
                  , ylim=(0, max(df_counts_unstack["4SU"] + 5)))
    plot_name = plot_dir + dir_sep + "joint_kde_plot_labels_{eff}.svg".format(eff=efficiency)
    plt.savefig(plot_name)

    plt.close(1)


# joint_scatter_plot_labels()
# joint_kde_plot_labels()


def plot_label_distribution(label, df_counts, strategy, nr_cells):

    return plot_distribution("real_count", df_counts, strategy, nr_cells, filter_label=True, label=label)


def plot_stationary_distribution(df_counts, strategy, nr_cells):

    return plot_distribution("count_all", df_counts, strategy, nr_cells)


# measure may be count_all or real_count (for label)
def plot_distribution(measure, df_counts, strategy, nr_cells, filter_label=False, label=""):

    params = sr.get(strategy=strategy)

    df_distribution, real_mean = create_simulated_distribution(measure, df_counts, strategy, nr_cells,
                                                               filter_label, label)

    # distribution from simulation
    plt.step(df_distribution[measure], df_distribution.chance, where="post")

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


def create_simulated_distribution(measure, df_counts, strategy, nr_cells, filter_label, label):

    if filter_label:
        df_counts = df_counts[df_counts.label == label]

    df_allele_cell_counts = df_counts.groupby(['allele_id', 'strategy', 'cell_id'])[measure].max().reset_index()
    df_one_allele_counts = df_allele_cell_counts[df_allele_cell_counts.strategy == strategy]
    df_one_allele_counts = df_one_allele_counts.set_index('cell_id'). \
        reindex(range(1, nr_cells + 1)).fillna(0).reset_index()

    df_distribution = df_one_allele_counts.groupby(measure)['cell_id'].count().to_frame().reset_index()
    df_distribution[measure] = df_distribution[measure].astype(int)

    max_count = df_distribution[measure].max()

    df_distribution = df_distribution.set_index(measure).reindex(range(0, max_count + 1)).fillna(0).reset_index()

    # cell_id contains the number of cells with the "measure" value (measure may be count_all or real_count (for label))
    nr_cells = int(sum(df_distribution.cell_id))
    df_distribution["chance"] = df_distribution.cell_id / nr_cells

    # total_chance = df_distribution["chance"].sum() # for debugging, should add up to 1
    df_distribution['weighted'] = df_distribution[measure] * df_distribution.cell_id
    sum_weighted = df_distribution['weighted'].sum()
    real_mean = round_sig(sum_weighted / nr_cells, 4)

    return df_distribution, real_mean


def plot_distributions():

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

    for strategy in strategies:

        plot_stationary_distribution(df_counts, strategy, nr_cells)
        plot_label_distribution("4SU", df_counts, strategy, nr_cells)


plot_distributions()

# now we want to calculate the mean and the variance of the number of molecules
# and see how well it predicts the dynamic parameters
# x-axis: time
# y-axis: mean


