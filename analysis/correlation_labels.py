# analysis of correlation between normalized label counts

# input: all counts for 500 cells for all combinations of:
# efficiencies = [1, 0.5, 0.2, 0.05]
# window_lengths (minutes) = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195]
#
# output: directory correlation_labels.plots containing
#   - a csv corr_labels_<len_win>.csv with Pearson correlation values
#   - some plots with correlation values for different categories (tran_type S/F and length of period)
#   - phase diagrams of differences in Pearson correlation p-value distribution for tran_type S and F
#       for dimensions win_len and efficiency, separately for every oscillation period
import os
import pandas as pd
import seaborn as sns
from simulator.StrategyReader import StrategyReader
from analysis.data_analysis import *
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2

efficiency = 1

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    # out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs_on_server_{}".format(efficiency)
    prj_dir = r"D:\26 Battich Oudenaarden transcriptional bursts"
else:
    dir_sep = "/"
    out_dir = "sc_runs_{}".format(efficiency)
    prj_dir = "sc_runs"

gap = 0
label_1 = "EU"
label_2 = "4SU"


def scatter_plot_for_allele(df_counts, strategy):

    df_counts_allele = df_counts[df_counts.strategy == strategy]

    corr, p_value = pearsonr(df_counts_allele.norm_count_1, df_counts_allele.norm_count_2)

    corr_s = '%.3f' % corr
    # print("Pearsons correlation: {}".format(corr_s))

    plt.title("strategy: {} with correlation {} (len_win={})".format(strategy, corr_s, len_win))
    plt.scatter(df_counts_allele.norm_count_1, df_counts_allele.norm_count_2)
    plt.xlabel("normalized counts label 1")
    plt.ylabel("normalized counts label 2")

    fig_name = plot_dir + dir_sep + "scatter_labels_{strategy}_{len_win}.svg".format(strategy=strategy, len_win=len_win)
    plt.savefig(fig_name)
    plt.close(1)


def calc_pearson_correlation(series):

    # criterium for at least 15 non-zero values in series.norm_count_1 and 15 in series.norm_count_2
    nr_label_1 = len(series[series.norm_count_1 > 0])
    nr_label_2 = len(series[series.norm_count_2 > 0])

    nr_data_points = np.min([nr_label_1, nr_label_2])

    if nr_data_points >= 15:
        corr, p_value = pearsonr(series.norm_count_1, series.norm_count_2)
    else:
        corr = np.nan
        p_value = 0

    df = pd.DataFrame([[corr, nr_data_points, p_value]], columns=["corr", "nr_data_points", "p_value"])

    return df


def calculate_corr_and_save(df_counts, len_win):

    df_corr = df_counts.groupby(["strategy"]).apply(calc_pearson_correlation).reset_index(). \
        sort_values("corr").drop('level_1', axis=1)

    df_corr = df_corr[(df_corr["corr"].notna()) & (df_corr.nr_data_points >= 15)]

    filename_corr_labels = plot_dir + dir_sep + "corr_labels_{len_win}.csv".format(len_win=len_win)

    # do we still want to save this separately?
    df_corr.to_csv(filename_corr_labels, sep=";", index=False)
    return df_corr


def add_coord_group_to_strategy(df_alleles):

    df_merged = df_alleles.merge(df_strategies, how='inner',
                                 left_on=['strategy'],
                                 right_on=['name'])
    df_merged["display"] = np.where( df_merged.coord_group.isnull(),
                                     df_merged['strategy'],
                                     df_merged['strategy'] + "__cg" + df_merged.coord_group.astype(str))

    return df_merged


def make_box_plot_for_periods(df, measure, agg_field, tran_type):

    pattern = "_" + tran_type
    df_type = df[df.strategy.str.contains(pattern)]

    plt.figure(figsize=(12, 5))
    plt.title("Pearson correlation for different periods, type={}".format(tran_type))
    sns.set(style="ticks")

    b = sns.boxplot(x=measure, y=agg_field, data=df_type,
                    palette="vlag")

    sns.set(font_scale=0.8)

    sns.swarmplot(x=measure, y=agg_field, data=df_type,
                  size=2, color=".3", linewidth=0)

    plt.xlabel("Pearson correlation between normalized counts two labels")
    plt.ylabel("Period (hours)")
    fig_name = plot_dir + dir_sep + "boxplot_correlation_{type}_{len_win}.svg".format(
        type=tran_type, len_win=len_win)
    plt.savefig(fig_name)
    plt.close(1)


def make_box_plot_for_len_win(df, period, measure, agg_field, efficiency):

    sns.set_theme(style="whitegrid")

    df = df.copy(deep=True)
    df["len_win_str"] = df.len_win.map(str)
    df["k_syn_str"] = df.k_syn.map(str)

    df_0_1 = df[df.k_syn == 0.1]
    df_1 = df[df.k_syn == 1]

    plt.figure(figsize=(12, 5))
    plt.title("Pearson correlation for different window lengths, period={}h efficiency={}".
              format(period, efficiency))

    gray_white = ['lightgrey', 'white']
    b = sns.boxplot(x=agg_field, y=measure, data=df, hue="tran_type", hue_order=["S", "F"],
                    palette=gray_white, orient="v")

    sns.set(font_scale=0.8)

    blues = ['steelblue', 'darkblue']
    reds = ['darksalmon', 'red']
    sw1 = sns.swarmplot(x=agg_field, y=measure, data=df_0_1, hue="tran_type", dodge=True,
                        size=2, palette=blues, linewidth=0, orient="v", hue_order=["S", "F"])

    sw2 = sns.swarmplot(x=agg_field, y=measure, data=df_1, hue="tran_type", dodge=True,
                        size=2, palette=reds, linewidth=0, orient="v", hue_order=["S", "F"])

    handles2, _ = sw2.get_legend_handles_labels()

    plt.xlabel("Length window (minutes)")
    plt.ylabel("Pearson correlation between normalized counts two labels")

    plt.legend(handles2, ['Stochastic (S)', 'Fluctuating (F)', 'S low', 'F low', 'S high', 'F high'])
    fig_name = plot_dir + dir_sep + "boxplot_correlation_{agg_field}_{period}.svg".format(
        agg_field=agg_field, period=period)
    plt.savefig(fig_name)
    plt.close(1)


def run_and_plot_one_correlation(df_counts, len_win):

    df_counts = normalize_counts(df_counts)
    df_counts_12 = merge_label_counts(df_counts, label_1, label_2)

    df_corr = calculate_corr_and_save(df_counts_12, len_win)

    df_corr = add_coord_group_to_strategy(df_corr)

    strategy = df_corr.head(1).strategy.item()
    scatter_plot_for_allele(df_counts_12, strategy)

    strategy = df_corr.tail(1).strategy.item()
    scatter_plot_for_allele(df_counts_12, strategy)

    df_corr["period"] = df_corr.strategy.str.split("_", expand=True)[4]

    for tran_type in ["S", "F"]:
        make_box_plot_for_periods(df_corr, measure="corr", agg_field="period", tran_type=tran_type)


# combine all correlation files into one and add an extra column for window size
def run_all_correlations(window_lengths):
    corr_list = []

    for len_win in window_lengths:

        filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
            len_win=len_win, gap=gap)

        # df_counts are counts of two labeling windows (single window length)
        df_counts = pd.read_csv(filename_counts, sep=';')

        df_counts = normalize_counts(df_counts)
        df_counts_12 = merge_label_counts(df_counts, label_1, label_2)

        df_counts_12["len_win"] = len_win

        df_corr = calculate_corr_and_save(df_counts_12, len_win)
        df_corr = add_coord_group_to_strategy(df_corr)
        df_corr["len_win"] = len_win

        corr_list.append(df_corr)

    df_corr_all = pd.concat(corr_list)
    df_corr_all["period"] = df_corr_all.strategy.str.split("_", expand=True)[4]

    return df_corr_all


def mwu_test_for_all_efficiencies_win_lens_and_periods(prj_dir, efficiencies, window_lengths, periods):
    mw_values = []
    # once we calculated the correlations for all efficiencies
    # we can compare F to S for all window sizes and all efficiencies for a certain window length

    for eff in efficiencies:

        out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)
        plot_dir = out_dir + dir_sep + "correlation_labels.plots"

        corr_name = "{plot_dir}{dir_sep}df_corr_all_G{gap}.csv".format(
            plot_dir=plot_dir, dir_sep=dir_sep, gap=gap)
        df_corr_all_eff = pd.read_csv(corr_name, sep=';')

        for period in periods:
            df_one_period = df_corr_all_eff[df_corr_all_eff.period == period]

            df_filtered = \
                df_corr_all_eff[(df_corr_all_eff.period == period) & (df_corr_all_eff.nr_data_points >= 15)]

            # now for each window length, calculate the difference in mean
            # between tran_type S and tran_type F
            for len_win in window_lengths:

                # here we could also average the number of data points for each window length
                mean_nr_data_points = int(df_one_period[df_one_period.len_win == len_win].nr_data_points.mean())

                df_s = df_filtered[(df_filtered.len_win == len_win) & (df_filtered.tran_type == "S")]
                df_f = df_filtered[(df_filtered.len_win == len_win) & (df_filtered.tran_type == "F")]

                mw_result = mannwhitneyu(x=df_s["corr"], y=df_f["corr"])

                mw_values.append([period, eff, len_win, mw_result.pvalue, mean_nr_data_points])

    df_mw_return = pd.DataFrame(mw_values, columns=["period", "eff", "len_win", "mw_pvalue", "mean_nr_data_points"])

    df_mw_return["minus_log10_mw_pvalue"] = -np.log10(df_mw_return.mw_pvalue).round(1)

    return df_mw_return


def create_phase_diagrams(df_mw_values, periods, prj_dir):

    for period in periods:
        data = df_mw_values[df_mw_values.period == period]

        measure = 'minus_log10_mw_pvalue'
        title = "-log10(pvalue) ManWU test between label correlation values of F and S for period: {}h".format(period)

        create_phase_diagram(data, 'eff', measure, title, period, prj_dir)

    measure = 'mean_nr_data_points'
    title = "mean nr of data points (=minimum of counts of label 1 and 2)"

    # for the mean number of data points we aggregate over the periods
    df_mw_agg_periods = df_mw_values.groupby(['eff', 'len_win']).agg(
        mean_nr_data_points=('mean_nr_data_points', 'mean'),
    ).round(0).reset_index()

    create_phase_diagram(df_mw_agg_periods, 'eff', measure, title, period="", prj_dir=prj_dir)


def create_phase_diagram(data, h_axis, measure, title, period, prj_dir):

    data = data[[h_axis, 'len_win', measure]]
    data = data.sort_values("len_win")
    data = data.set_index([h_axis, "len_win"])

    # convert from multi-index to cross-product table
    data = data.unstack().fillna(0)
    # data = data.transpose() # if you would like to switch columns and rows
    # rename columns, unstack
    data.columns = [x[1] for x in data.columns.ravel()]

    plt.figure(figsize=(12, 5))
    ax = sns.heatmap(data, cmap="Spectral_r", annot=True, fmt='g')
    plt.title(title)
    plot_dir = r"{}{}runs{}{}".format(prj_dir, dir_sep, dir_sep, "correlation_labels.plots")

    phase_plot_name = plot_dir + dir_sep + "phase_plot_{}_{}.svg".format(measure, period)

    plt.xlabel("window size in minutes (no gap)")

    if h_axis == "eff":
        plt.ylabel("efficiency")
    elif h_axis == "half_life_h":
        plt.ylabel("half life")
    plt.savefig(phase_plot_name)
    plt.close(1)


len_win = 120

plot_dir = out_dir + dir_sep + "correlation_labels.plots"
os.makedirs(plot_dir, exist_ok=True)

filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap)

# df_counts are counts of two labeling windows (single window length)
df_counts = pd.read_csv(filename_counts, sep=';')

strategies_file = out_dir + dir_sep + "strategies_mixed_new.csv"
sr = StrategyReader(strategies_file)
# sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies

run_and_plot_one_correlation(df_counts, len_win)

corr_name = "{od}{dir_sep}df_corr_all_G{gap}.csv".format(
    od=plot_dir, dir_sep=dir_sep, gap=gap)

window_lengths = [15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195]
run_corr = False

if run_corr:
    df_corr_all = run_all_correlations(window_lengths)
    df_corr_all.to_csv(path_or_buf=corr_name, sep=';', index=False)
else:
    df_corr_all = pd.read_csv(corr_name, sep=';')

periods = [1, 2, 3, 5, 12, 24]
for period in periods:
    # df_one_period = df_corr_all[(df_corr_all.period == period) & (df_corr_all.nr_data_points >= 30)]
    df_one_period = df_corr_all[df_corr_all.period == period]

    # for 1 efficiency: make box plot and compare between F (deterministic fluctuation) and S (stochastic behavior)
    make_box_plot_for_len_win(df_one_period, period, "corr", "len_win_str", efficiency)

# after having run correlation for all window lengths
efficiencies = [1, 0.5, 0.2, 0.05]
df_mw_values = mwu_test_for_all_efficiencies_win_lens_and_periods(prj_dir, efficiencies, window_lengths, periods)

create_phase_diagrams(df_mw_values, periods, prj_dir)

# based on one df_corr_all for one efficiency
# we want to determine the percentage of genes for which correlation can be calculated

corr_list_for_eff = []
for eff in efficiencies:
    out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)
    plot_dir = out_dir + dir_sep + "correlation_labels.plots"

    corr_name = "{plot_dir}{dir_sep}df_corr_all_G{gap}.csv".format(
        plot_dir=plot_dir, dir_sep=dir_sep, gap=gap)
    df_corr_all_eff = pd.read_csv(corr_name, sep=';')
    df_corr_all_eff["eff"] = eff
    corr_list_for_eff.append(df_corr_all_eff)

data = pd.concat(corr_list_for_eff)
data["half_life_h"] = (np.log(2)/(60 * data["k_d"])).round(1).map(str)

# count number of half life categories over all genes
df_hl_counts = data.groupby(['strategy', 'half_life_h']).count().reset_index()['half_life_h'].value_counts()

nr_genes = len(data.strategy.unique())

# data = data[data.p_value < 0.05]
data = data[(data.p_value < 0.05) & (data.nr_data_points > 30)]

df_agg = data.groupby(['eff', 'len_win', 'half_life_h']).agg(
    corr_count=('corr', 'count')
).reset_index()

df_agg["perc"] = (100 * df_agg["corr_count"] / nr_genes).round(1)

df_agg2 = df_agg.groupby(['eff', 'len_win']).agg(
    corr_count=('corr_count', 'sum')
).reset_index()

df_agg2["perc"] = (100 * df_agg2["corr_count"] / nr_genes).round(1)

create_phase_diagram(df_agg2, "eff", "perc",
                     "percentage of genes for which correlation can be calculated",
                     "", prj_dir)

# now take one row, e.g. efficiency 20%
df_agg = df_agg[df_agg.eff == 0.2]

df_agg_plus_norm = df_agg.merge(df_hl_counts, how="inner",
                                left_on='half_life_h',
                                right_index=True)

# normalize based on the total number of genes per half_life
df_agg_plus_norm['corr_count_perc'] = (100 * df_agg_plus_norm.corr_count / df_agg_plus_norm.half_life_h_y).round()

# and split per half life category
create_phase_diagram(df_agg_plus_norm, "half_life_h", "corr_count_perc",
                     "percentage of genes for which correlation can be calculated (eff=20%)",
                     "", prj_dir)
