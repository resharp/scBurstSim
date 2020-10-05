import os
from simulator.Experiment import *
from simulator.data_analysis import *
import logging
import sys

# get the fully-qualified logger (here: `root.__main__`)
logger = logging.getLogger(__name__)

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own working directory for locally storing data sets
    work_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    work_dir = "."

start_windows = 600; length_window = 60; between_window = 15
window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
window_4su = [start_windows + length_window + between_window,
              start_windows + 2*length_window + between_window, '4SU'] # e.g. 120 minutes of EU labeling
windows = [window_eu, window_4su]
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = windows[-1][WINDOW_END] + 0  # freeze 0 minutes after end of last window

# under this run_dir we should also create a plot directory
run_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"

strategies_file = run_dir + dir_sep + "strategies.csv"
run_sim = False # setting run_sim to False results in use of locally stored data set
nr_cells = 10

# see strategy names in data\strategies.csv


def main():

    logging.basicConfig(filename=work_dir + dir_sep + 'main_scBurstSim.log', filemode='w',
                        # format='%(asctime)s - %(message)s',
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    logger.info("scBurstSim started")

    # TODO: Add argparse

    nr_cells = 100
    nr_syn_within_strategy = 10
    nr_non_syn_within_strategy = 10

    logger.info("scBurstSim started for {nr_cells} cells".format(nr_cells=nr_cells))

    exp_params = ExperimentParams(nr_cells=nr_cells,
                                  nr_syn_within_strategy=nr_syn_within_strategy,
                                  nr_non_syn_within_strategy=nr_non_syn_within_strategy,
                                  windows=windows, freeze=freeze)

    strategies_file = work_dir + dir_sep + "strategies.csv"
    exp = Experiment(exp_params, strategies_file)

    filename = "{wd}{dir_sep}df_counts".format(wd=work_dir, dir_sep=dir_sep)
    if run_sim:

        df_counts = exp.run()
        df_counts.to_csv(path_or_buf=filename, sep='\t', index=False)
    else:
        df_counts = pd.read_csv(filename, sep='\t')

    print("Experiment run. Number of counts: {counts}.".format(counts=len(df_counts)))
    logging.info("Experiment run. Number of counts: {counts}.".format(counts=len(df_counts)))

    label = "EU"
    df_counts_eu = df_counts[df_counts.label == label]

    # df_counts_eu = violin_plot_fraction(0.8, "80", df_counts_eu)

    # do_kolmogorov_smirnov_tests_for_percentages_on(df_counts_eu)

    # TODO: df_all_arrivals can be used for sampling (it still contains information on single molecule level)
    df_all_transcripts = exp.df_all_transcripts

    # what is the distribution of fractions?
    density_plot("fraction", "strategy", df_counts_eu, exp_params)

    # try_out_logistic_regression(perc="50", df_counts_label=df_counts_eu)

    # regression_plot("perc_label_on", "fraction", df_counts_eu, exp_params)

    # density_plot("perc_label_on", "strategy", df_counts_eu, exp_params)
    # regression_plot("real_count_unlabeled", "real_count", df_counts_eu, exp_params)

    # show_distribution_real_counts(df_counts, nr_cells)

    # cluster map creates plot cluster_map.svg in run directory if you do not provide a plot name
    label = "EU"
    cluster_map(df_counts, label=label, plot_name=work_dir + dir_sep + "cluster_map_{label}.svg".format(label=label))


if __name__ == "__main__":
    main()

