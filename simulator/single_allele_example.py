from simulator.Transcription import *
from simulator.transcription_plots import *
import numpy as np
nr_days = 1
max_minutes = 1440*nr_days # 24 hours = 1440 minutes
# part of Experiment (Uridine analog windows)

# windows = [[400, 460, 'EU'], [520, 580, '4SU']] # e.g. 120 minutes of EU labeling
start_window = 600
windows = [[start_window, start_window + 120, 'EU']] # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = windows[-1][WINDOW_END] + 0  # freeze 30 minutes after end of last window

k_10 = 0.02
k_01 = k_10
k_syn = 0.16
k_d = 0.01
half_life = int(np.log(2)/k_d); mean_life = int(1/k_d)
params = TranscriptParams(k_01=k_01, k_10=k_10, k_syn=k_syn, nr_refractions=2, k_d=k_d)

trans = Transcription(params)

df_dtmc, df_events = trans.run_bursts(max_minutes, windows)

df_unlabeled_events = trans.df_unlabeled_events
df_labeled_events = trans.df_labeled_events

# calculate average burst size
mean_burst_size = df_dtmc[df_dtmc.state == "1"].burst_size.mean().round(1)
std_burst_size = df_dtmc[df_dtmc.state == "1"].burst_size.std().round(1)
nr_bursts = len(df_dtmc[df_dtmc.state == "1"])
burst_frequency = round(nr_bursts/max_minutes, 3)

title = "k_01={k_01}; k_syn={k_syn}; k_d={k_d} -> burst size: {bs} +/- {std}; burst freq: {freq}".format(
    bs=mean_burst_size, std=std_burst_size
    , freq=burst_frequency
    , k_01=params.k_01, k_syn=params.k_syn, k_d=params.k_d)

plot_events(df_dtmc, df_events)
plot_dynamics(title=title, df_events=df_unlabeled_events, freeze=freeze, max_minutes=max_minutes,
              windows=windows, df_labeled_arrivals=df_labeled_events)

# show distribution of time in burst and of time until burst (in minutes)
# plot_waiting_time_distribution(df_dtmc)


class StrategyReader:
    filename = ""

    def __init__(self, filename):
        self.filename = filename
