import os

from simulator.Experiment import *
from simulator.StrategyReader import StrategyReader
from simulator.data_analysis import *

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "cluster_alleles.plots"
os.makedirs(plot_dir, exist_ok=True)

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

efficiency = 100

gap = 0
label_1 = "EU"
label_2 = "4SU"
len_win = 60

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap, eff=efficiency)


def normalize(df_counts):

    df_mean_counts = df_counts[['allele_id', 'label', 'real_count']].\
        groupby(['allele_id', 'label']).mean().reset_index()

    df_mean_counts.rename(columns={'real_count': 'mean_count'}, inplace=True)

    df_ret = df_counts.merge(df_mean_counts, how='inner',
                             left_on=['allele_id', 'label'],
                             right_on=['allele_id', 'label'])

    df_ret["norm_count"] = df_ret["real_count"] / df_ret["mean_count"]

    return df_ret


# df_counts are counts of two labeling windows (single window length)
df_counts = pd.read_csv(filename_counts, sep=';')

df_counts = normalize(df_counts)

strategies_file = out_dir + dir_sep + "strategies_mixed.csv"
sr = StrategyReader(strategies_file)
# sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies

df_counts["real_count_log10"] = np.log10(df_counts["real_count"])
df_counts["norm_count_log10"] = np.log10(df_counts["norm_count"])
df_counts["count_all_log10"] = np.log10(df_counts["count_all"])


windows, fix_time = get_windows_and_fix_time(length_window=len_win, gap=gap)

exp_params = ExperimentParams(nr_cells=100,
                              strategies_file=strategies_file,
                              nr_syn_within_strategy=1,
                              nr_non_syn_within_strategy=1,
                              efficiency=1,
                              windows=windows, fix_time=fix_time)

for window in windows:
    # cluster map creates plot cluster_map.svg in run directory if you do not provide a plot name
    label = window[WINDOW_LABEL]

    # only use label counts for clustering
    measures = ["real_count_log10", "norm_count"]

    for measure in measures:
        # TODO: Be careful: the sum of count_all differs for different labels due to filtering in cluster_map
        # because rows may be missing for the label /allele/cell combination resulting in missing that count_all
        cluster_map(df_counts, measure=measure, label=label, exp_params=exp_params,
                    plot_name=plot_dir + dir_sep + "{measure}_cluster_map_{label}.svg".
                    format(label=label, measure=measure))
