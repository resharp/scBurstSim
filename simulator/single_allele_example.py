from simulator.Transcription import *
from simulator.transcription_plots import *

import os
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
start_window = 600
windows = [[start_window, start_window + 120, 'EU']] # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = windows[-1][WINDOW_END] + 0  # freeze 30 minutes after end of last window


def run_example(params):

    trans = Transcription(params)

    df_dtmc, df_events = trans.run_bursts(max_minutes, windows)

    df_unlabeled_events = trans.df_unlabeled_events
    df_labeled_events = trans.df_labeled_events

    # calculate average burst size
    nr_bursts = len(df_dtmc[df_dtmc.state == "1"])
    burst_frequency = round(nr_bursts / max_minutes, 3)

    title = "strategy={name}; k_01={k_01}; k_10={k_10};k_syn={k_syn}; k_d={k_d} -> " \
            "burst freq: {freq}".format(name=params.name, freq=burst_frequency,
                                        k_01=params.k_01, k_10=params.k_10, k_syn=params.k_syn, k_d=params.k_d)
    plot_events(df_dtmc, df_events)
    plot_dynamics(title=title, df_events=df_unlabeled_events, freeze=freeze, max_minutes=max_minutes,
                  windows=windows, df_labeled_arrivals=df_labeled_events)

    # show distribution of time in burst and of time until burst (in minutes)
    # plot_waiting_time_distribution(df_dtmc)


def run_all_strategies():
    global params
    for params in sr.select_all():
        run_example(params)


sr = StrategyReader(work_dir + dir_sep + "strategies.csv" )

# see strategy names in data\strategies.csv
# params = sr.get(strategy="frequent")
params = sr.get_random()

# params = sr.get_random()
# run_example(params)

run_all_strategies()


