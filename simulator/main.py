from simulator.Experiment import *
from simulator.Transcription import *
from simulator.data_analysis import *

max_minutes = 1440  # 24 hours = 1440 minutes
windows = [[400, 520, 'EU']]  # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = 550  # freeze 30 minutes after end of last window

params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)

nr_cells = 1000
nr_alleles = 1
exp = Experiment(nr_cells, nr_alleles, params, windows, freeze)

df_counts = exp.run()

print("Experiment run. Number of counts: {counts}.".format(counts=len(df_counts)))

df_counts["fraction"] = df_counts["real_count"] / (df_counts["real_count"] + df_counts["real_count_unlabeled"])

df_counts_eu = df_counts[df_counts.label == "EU"].copy(deep=True)

df_counts_eu = violin_plot_fraction(0.8, "80", df_counts_eu)

# do_kolmogorov_smirnov_tests_for_percentages_on(df_counts_eu)

df_all_arrivals = exp.df_all_arrivals

# try_out_logistic_regression(perc="50", df_counts_label=df_counts_eu)

regression_plot("perc_label_on", "fraction", df_counts_eu)

# show_distribution_real_counts(df_counts, nr_cells)
