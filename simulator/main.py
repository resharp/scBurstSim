import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Transcription import Transcription

# part of Experiment (analog windows)
max_minutes = 1440 # 24 hours = 1440 minutes
windows = [[400, 520, 'EU']] # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
wash_out = 490

# part of Transcription model (basis parameters)
# two transitions
l_01 = 0.02
l_10 = 0.02
k_syn = 0.16  # k_syn = synthesis rate = transcription rate
nr_refractions = 1
# decay is strictly not part of transcription but we include it in the model
k_d = 0.01   # k_d = decay rate


trans = Transcription(l_01, l_10, k_syn, nr_refractions, k_d)


def plot_events(df, df_poisson_arrivals):
    plt.step(df["begin_time"], df["state"], where="post", color="m")

    arrival_list = df_poisson_arrivals[df_poisson_arrivals.count_s > 0]['arrival'].to_list()
    y_arrivals = [1] * len(arrival_list)

    decay_list = df_poisson_arrivals[df_poisson_arrivals.count_s < 0]['arrival'].to_list()
    y_decays = [0] * len(arrival_list)

    plt.scatter(arrival_list, y_arrivals, color='m', marker="d", s=9)
    plt.scatter(decay_list, y_decays, color='r', marker="o", s=9)


def plot_dynamics(df_poisson_arrivals, windows = [], df_label_arrivals = []):

    plt.title("l_01={l_01}; k_syn={k_syn}; k_d={k_d} -> burst size: {bs} +/- {std}; burst freq: {freq}".format(
        bs=mean_burst_size, std=std_burst_size
        , freq=burst_frequency
        , l_01 = l_01, k_syn=k_syn, k_d=k_d))

    plt.step(df_poisson_arrivals["arrival"], df_poisson_arrivals["cum_count"], where="post", color="tab:blue")

    colors = ["darkgreen", "m"]; color = 0
    for df_label_arrival in df_labeled_arrivals:
        plt.step(df_label_arrival["arrival"], df_label_arrival["cum_count"], where="post", color=colors[color])
        color = color + 1

    plt.xlim(0, max_minutes)
    plt.xlabel("minutes")
    plt.ylabel("nr of transcripts")

    for window in windows:
        start_window = window[WINDOW_START]
        end_window = window[WINDOW_END]

        plt.axvline(start_window, label='labeling window', c="r")
        plt.axvline(end_window, c="r")

    plt.legend()

    plt.show()


def plot_waiting_time_distribution(df):
    df_off = df[df.state == "0"]
    df_on = df[df.state == "1"]

    fig, ax = plt.subplots()
    ax.hist(df_on["state_time"], color="r", bins=60, label="active (ON)")
    ax.hist(df_off["state_time"], color="b", bins=60, label="silent (OFF)")
    legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

    plt.title("Distribution of ON and OFF time intervals")
    plt.show()


dtmc_list = []

state = "0"
current_time = 0
end_time = 0
poisson_arrivals = []

# this main part is part Of either Experiment or Transcription
# l_01 and l_10 are from the Transcription Model
while current_time < max_minutes:

    state_time = 0
    burst_size = 0
    if state == "0":
        l = l_01
        # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
        # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
        # this is a simple way to simulate multiple refractory states?
        alpha = nr_refractions
        for i in range(alpha):
            state_time = state_time + np.random.exponential(scale=1.0, size=None) / l
    else:
        l = l_10
        state_time = state_time + np.random.exponential(scale=1.0, size=None) / l

        # current_time is the start of active (ON) state
        # state_time is length of burst
        # now create new Poisson arrivals
        new_arrivals = trans.new_poisson_arrivals(current_time, state_time, windows)
        burst_size = len(new_arrivals)
        poisson_arrivals = poisson_arrivals + new_arrivals

    end_time = current_time + state_time
    dtmc_list.append([state, current_time, end_time, state_time, burst_size])

    current_time = end_time

    # switch state (can be done on Transcription)
    if state == "0":
        state = "1"
    else:
        state = "0"

df = pd.DataFrame(data=dtmc_list, columns=["state", "begin_time", "end_time", "state_time", "burst_size"])

df_poisson_arrivals = pd.DataFrame(poisson_arrivals, columns=["label", "arrival", "count_s", "decay", "count_d"])

# we will now put the arrivals and decays in one table and sort by time
df_decays = df_poisson_arrivals[['label','decay', "count_d"]].\
    rename(columns={'decay': 'arrival', 'count_d': 'count_s'})

df_poisson_arrivals = df_poisson_arrivals[["label", "arrival", "count_s"]]
df_poisson_arrivals = df_poisson_arrivals.append(df_decays).sort_values(by="arrival")

plot_events(df, df_poisson_arrivals)

df_labeled_arrivals = []
for window in windows:
    df_labeled = df_poisson_arrivals[ df_poisson_arrivals.label == window[WINDOW_LABEL]]
    df_labeled['cum_count'] = df_labeled['count_s'].cumsum()
    df_labeled_arrivals.append(df_labeled)

df_poisson_arrivals = df_poisson_arrivals[df_poisson_arrivals.label == ""]
df_poisson_arrivals['cum_count'] = df_poisson_arrivals['count_s'].cumsum()

debug = "True"

# calculate average burst size
mean_burst_size = df[df.state == "1"].burst_size.mean().round(1)
std_burst_size = df[df.state == "1"].burst_size.std().round(1)
nr_bursts = len(df[df.state == "1"])
burst_frequency = round(nr_bursts/max_minutes, 3)

plot_dynamics(df_poisson_arrivals, windows, df_labeled_arrivals)

# plot_waiting_time_distribution(df)
