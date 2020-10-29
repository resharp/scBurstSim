import os
import pandas as pd
import matplotlib.pyplot as plt
from simulator.StrategyReader import StrategyReader
import numpy as np
import seaborn as sns
# we would like to analyse how well we can infer the strategy parameters from
# the counts
# we read the counts from the df_counts.csv generated by main
# and then join the generated strategies strategies_generated.csv (use df_strategies of StrategyReader)
# then for each allele we calculate the fraction of cells for which a count is detected
# (and not detected)

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "analyze_parameters.plots"
os.makedirs(plot_dir, exist_ok=True)

sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
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


nr_cells = 100
efficiency = 10
label = "4SU"  # 2nd window
len_win = 60
gap = 0

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap, eff=efficiency)

df_counts = pd.read_csv(filename_counts, sep=';')

# mean expression
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

# example
t = len_win + gap # length window + gap

y_w2 = 100
x_w1 = 10
k_d = (np.log(y_w2) - np.log(x_w1)) / t

# we can add strategy to the "index" because it is completely defined by key allele_id
df_counts_unstack = df_counts.set_index(["cell_id", "allele_id", "strategy", "label"])['real_count'].unstack()
df_counts_unstack = df_counts_unstack.reset_index().fillna(0)

# what zeroes to include?
# we want all the counts per cell and allele for label_1
# where label_2 does not have a zero
label_1 = "EU"; label_2 = "4SU"

# so we need the mean value of the EU label and the mean value of the 4SU label for every allele
# however, we may underestimate the 1st label, because there will be no counts for the 0 counts
# so we should divide the total of the 1st label by the number of counts of the 2nd label?
# we tried that with the following line, but prediction is actually worse
# df_counts_unstack = df_counts_unstack[df_counts_unstack[label_2] != 0]

# calculate mean expression level
df_allele_counts = df_counts_unstack[['allele_id', 'strategy', label_1, label_2]].\
    groupby(['allele_id', 'strategy']).mean().reset_index()

df_allele_counts = pd.merge(df_allele_counts, df_strategies, how="left",
                            left_on=['strategy'],
                            right_on=['name'])

k_d = (np.log(y_w2) - np.log(x_w1)) / t

df_allele_counts['k_d_predicted'] = (np.log(df_allele_counts[label_2]) - np.log(df_allele_counts[label_1]))/t

df_allele_counts['k_d_error'] = df_allele_counts['k_d_predicted'] / df_allele_counts['k_d']


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


plot_predicted_k_d()
plot_error_k_d()

strategy = "generated_8"
# strategy = "generated_27"
# we would like to make a 2-dim density plot of the two labels

df_counts_unstack = df_counts_unstack[df_counts_unstack.strategy == strategy]

pseudocount = 0.1
df_counts_unstack["log10_4SU"] = np.log10(df_counts_unstack["4SU"] + pseudocount)
df_counts_unstack["log10_EU"] = np.log10(df_counts_unstack["EU"] + pseudocount)

sns.jointplot(x=df_counts_unstack["4SU"],
              y=df_counts_unstack["EU"],
              kind='scatter', s=50, color='b')
plt.show()
plt.close(1)

sns.set(style="white", color_codes=True)
sns.jointplot(x=df_counts_unstack["4SU"], y=df_counts_unstack["EU"], kind='kde', color="skyblue"
              , xlim=(0, max(df_counts_unstack["4SU"] + 5))
              , ylim=(0, max(df_counts_unstack["EU"] + 5)))
plt.show()
