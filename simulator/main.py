from simulator.Transcription import *
from simulator.transcription_plots import *

max_minutes = 1440 # 24 hours = 1440 minutes
# part of Experiment (Uridine analog windows)
windows = [[400, 520, 'EU']] # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
wash_out = 490

params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)
trans = Transcription(params)

df_dtmc, df_poisson_arrivals = trans.run_bursts(max_minutes, windows)


# we will now put the arrivals and decays in one table and sort by time
df_decays = df_poisson_arrivals[['label','decay', "count_d"]].\
    rename(columns={'decay': 'arrival', 'count_d': 'count_s'})

df_poisson_arrivals = df_poisson_arrivals[["label", "arrival", "count_s"]]
df_poisson_arrivals = df_poisson_arrivals.append(df_decays).sort_values(by="arrival")

plot_events(df_dtmc, df_poisson_arrivals)

df_labeled_arrivals = []

for window in windows:
    df_labeled = df_poisson_arrivals[df_poisson_arrivals.label == window[WINDOW_LABEL]]
    df_labeled['cum_count'] = df_labeled['count_s'].cumsum()
    df_labeled_arrivals.append(df_labeled)

df_poisson_arrivals = df_poisson_arrivals[df_poisson_arrivals.label == ""]
df_poisson_arrivals['cum_count'] = df_poisson_arrivals['count_s'].cumsum()

# calculate average burst size
mean_burst_size = df_dtmc[df_dtmc.state == "1"].burst_size.mean().round(1)
std_burst_size = df_dtmc[df_dtmc.state == "1"].burst_size.std().round(1)
nr_bursts = len(df_dtmc[df_dtmc.state == "1"])
burst_frequency = round(nr_bursts/max_minutes, 3)

title = "l_01={l_01}; k_syn={k_syn}; k_d={k_d} -> burst size: {bs} +/- {std}; burst freq: {freq}".format(
    bs=mean_burst_size, std=std_burst_size
    , freq=burst_frequency
    , l_01=params.l_01, k_syn=params.k_syn, k_d=params.k_d)

plot_dynamics(title=title, df_poisson_arrivals=df_poisson_arrivals, max_minutes=max_minutes,
              windows=windows, df_labeled_arrivals=df_labeled_arrivals )

# plot_waiting_time_distribution(df_dtmc)
