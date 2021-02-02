import os

from pandas import DataFrame
from sklearn.cluster import AgglomerativeClustering

from simulator.Experiment import *
from simulator.StrategyReader import StrategyReader
from simulator.data_analysis import *
import scipy.cluster.hierarchy as sch

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


def make_cluster_maps():

    windows, fix_time = get_windows_and_fix_time(length_window=len_win, gap=gap)

    nr_cells = 100
    efficiency = 1

    for window in windows:
        # cluster map creates plot cluster_map.svg in run directory if you do not provide a plot name
        label = window[WINDOW_LABEL]

        # only use label counts for clustering
        measures = ["real_count_log10", "norm_count"]

        for measure in measures:
            # TODO: Be careful: the sum of count_all differs for different labels due to filtering in cluster_map
            # because rows may be missing for the label /allele/cell combination resulting in missing that count_all
            cluster_map(df_counts, measure=measure, label=label, nr_cells=nr_cells, efficiency=efficiency,
                        plot_name=plot_dir + dir_sep + "{measure}_cluster_map_{label}.svg".
                        format(label=label, measure=measure))


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    sch.dendrogram(linkage_matrix, **kwargs)


def show_other_method_by_correlation():

    df = df_counts_2_matrix.transpose()    # we transpose because we want the alleles as features here

    cor = df.corr()
    filename_corr = plot_dir + dir_sep + "corr_matrix.csv"
    cor.to_csv(filename_corr, sep=";")

    sns.heatmap(cor, square=True, xticklabels=cor.columns, yticklabels=cor.columns)
    plt.show()


def cluster_for_dendrogram(df_matrix):
    # https://towardsdatascience.com/hierarchical-clustering-explained-e58d2f936323
    # first calculate all distances by setting distance_threshold=0
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, affinity='euclidean', linkage='ward')

    model = model.fit(df_matrix)

    print("nr of clusters: {}".format(model.n_clusters_))
    # If there are n distances greated than [distance cutoff] so, when combined, n+1 clusters will be formed
    distances = model.distances_

    plt.figure(figsize=(12, 6))
    plt.title('Hierarchical Clustering Dendrogram')
    # plot the top three levels of the dendrogram
    plot_dendrogram(model, truncate_mode='level', p=3)

    plot_name = plot_dir + dir_sep + "dendrogram.svg"
    plt.xlabel("Number of points in node (or index of point if no parenthesis).")

    plt.savefig(plot_name)


def cluster_alleles(df_matrix):
    # After having visually inspected the dendrogram we decide to set the distance cutoff at
    threshold = 20

    model = AgglomerativeClustering(distance_threshold=threshold, n_clusters=None, affinity='euclidean', linkage='ward')
    # and cluster again

    model = model.fit(df_matrix)

    print("nr of clusters: {}".format(model.n_clusters_))

    df_labels = DataFrame(model.labels_)

    # now we want to make a DataFrame with the combination of alleles and cluster labels

    df_matrix = df_matrix.reset_index()

    # merge on index
    df_merge = df_labels.merge(df_matrix, left_index=True, right_index=True, how="inner")

    df_merge.rename(columns={0: 'cluster'}, inplace=True)

    filename_clusters = plot_dir + dir_sep + "clusters_{th}.csv".format(th=threshold)
    df_merge.to_csv(filename_clusters, sep=";")

    # count number of allles in clusters
    df_cluster_counts = df_merge[['cluster', 'strategy']].groupby('cluster').count().reset_index()

    print(df_cluster_counts)

    filename_cluster_counts = plot_dir + dir_sep + "cl_counts_{th}.csv".format(th=threshold)
    df_cluster_counts.to_csv(filename_cluster_counts, sep=";")


# using seaborn cluster maps
make_cluster_maps()

# show_other_method_by_correlation()

df_counts_2 = df_counts[(df_counts.label == label_2)][['strategy', 'cell_id', 'norm_count']]

df_counts_2_matrix = df_counts_2.set_index(["strategy", "cell_id"])['norm_count'].unstack()
df_counts_2_matrix = df_counts_2_matrix.fillna(0)

cluster_for_dendrogram(df_counts_2_matrix)

cluster_alleles(df_counts_2_matrix)

