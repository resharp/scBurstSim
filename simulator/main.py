from simulator.Experiment import *
from simulator.Transcription import *
from simulator.data_analysis import *
import os
import pandas as pd
if os.name == 'nt':
    dir_sep = "\\"
    # to do: set your own working directory for locally storing data sets
    work_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"
else:
    dir_sep = "/"
    work_dir = "."

run_sim = True # run_sim False uses locally stored data set

max_minutes = 1440  # 24 hours = 1440 minutes
windows = [[400, 460, 'EU']]  # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = windows[-1][WINDOW_END] + 30  # freeze 30 minutes after end of last window

trans_params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)

nr_cells = 1000
nr_alleles = 1

exp_params = ExperimentParams(nr_cells=nr_cells, nr_alleles=nr_alleles, windows=windows, freeze=freeze,
                              trans_params=trans_params)

exp = Experiment(exp_params)

filename = "{wd}{dir_sep}df_counts".format(wd=work_dir, dir_sep=dir_sep)
if run_sim:
    df_counts = exp.run()
    df_counts.to_csv(path_or_buf=filename, sep='\t', index=False)
else:
    df_counts = pd.read_csv(filename, sep='\t')

print("Experiment run. Number of counts: {counts}.".format(counts=len(df_counts)))

df_counts["fraction"] = df_counts["real_count"] / (df_counts["real_count"] + df_counts["real_count_unlabeled"])

df_counts_eu = df_counts[df_counts.label == "EU"].copy(deep=True)

# df_counts_eu = violin_plot_fraction(0.8, "80", df_counts_eu)

# do_kolmogorov_smirnov_tests_for_percentages_on(df_counts_eu)

df_all_arrivals = exp.df_all_arrivals

# try_out_logistic_regression(perc="50", df_counts_label=df_counts_eu)

regression_plot("perc_label_on", "fraction", df_counts_eu, exp_params)

# show_distribution_real_counts(df_counts, nr_cells)
