from simulator.Transcription import *
from simulator.transcription_plots import *

import os

# script to display single cell single allele example traces

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own working directory for locally storing data sets
    work_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"
else:
    dir_sep = "/"
    work_dir = "."

nr_days = 1
max_minutes = 1440*nr_days # 24 hours = 1440 minutes
# part of Experiment (Uridine analog windows)

# windows = [[400, 460, 'EU'], [520, 580, '4SU']] # e.g. 120 minutes of EU labeling
start_windows = 600; length_window = 120; between_window = 15
window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
window_4su = [start_windows + length_window + between_window,
              start_windows + 2*length_window + between_window, '4SU'] # e.g. 120 minutes of EU labeling
windows = [window_eu, window_4su]
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = windows[-1][WINDOW_END] + 0  # freeze 30 minutes after end of last window


def run_example(params):

    trans = Transcription(params)

    # set complete_trace=True to retrieve the complete trace of transcripts counts (for plotting)
    df_dtmc, dtmc_list = trans.run_bursts(max_minutes, windows, new_dtmc_trace=True, complete_trace=True)

    # only set when complete_trace = True:
    # df_events, df_unlabeled_events and dfs_labeled_events
    df_events = trans.df_events
    df_unlabeled_events = trans.df_unlabeled_events
    dfs_labeled_events = trans.dfs_labeled_events

    # calculate average burst size
    nr_bursts = len(df_dtmc[df_dtmc.state == "1"])
    burst_frequency = round(nr_bursts / max_minutes, 3)

    title = "strategy={name}; k_01={k_01}; k_10={k_10};k_syn={k_syn}; k_d={k_d} -> " \
            "burst freq: {freq}".format(name=params.name, freq=burst_frequency,
                                        k_01=params.k_01, k_10=params.k_10, k_syn=params.k_syn, k_d=params.k_d)
    plot_events(df_dtmc, df_events)
    plot_dynamics(title=title, df_events=df_unlabeled_events, freeze=freeze, max_minutes=max_minutes,
                  windows=windows, dfs_labeled_events=dfs_labeled_events)

    # show distribution of time in burst and of time until burst (in minutes)
    # plot_waiting_time_distribution(df_dtmc)


def run_all_strategies():
    global params
    for params in sr.select_all():
        run_example(params)


sr = StrategyReader(work_dir + dir_sep + "strategies.csv" )

# see strategy names in data\strategies.csv

# we can select a strategy by name
# params = sr.get(strategy="frequent_high")

# or retrieve a random strategy
params = sr.get_random()
run_example(params)

# or run an example of all strategies
# run_all_strategies()


