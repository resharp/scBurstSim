import matplotlib.pyplot as plt

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


def plot_events(df_dtmc, df_events):
    plt.step(df_dtmc["begin_time"], df_dtmc["state"], where="post", color="m")

    arrival_list = df_events[df_events.count_s > 0]['arrival'].to_list()
    y_arrivals = [1] * len(arrival_list)

    decay_list = df_events[df_events.count_s < 0]['arrival'].to_list()
    y_decays = [0] * len(arrival_list)

    plt.scatter(arrival_list, y_arrivals, color='m', marker="d", s=9)
    plt.scatter(decay_list, y_decays, color='r', marker="o", s=9)


def plot_dynamics(title, df_events, freeze, max_minutes, windows=[], dfs_labeled_events=[]):

    plt.title(title)

    plt.step(df_events["arrival"], df_events["cum_count"], where="post", color="tab:blue")

    colors = ["darkgreen", "peru"]; color = 0
    for label, df_labeled_events in dfs_labeled_events:
        plt.step(df_labeled_events["arrival"], df_labeled_events["cum_count"], where="post", color=colors[color])
        color = color + 1

    plt.xlim(0, max_minutes)
    plt.xlabel("minutes")
    plt.ylabel("nr of transcripts")

    for window in windows:
        start_window = window[WINDOW_START]
        end_window = window[WINDOW_END]

        plt.axvline(start_window, label='labeling window', c="r")
        plt.axvline(end_window, c="r")

    plt.axvline(freeze, linestyle='--', color='black', label="freeze")

    plt.legend()

    plt.show()


def plot_waiting_time_distribution(df):

    df_off = df[df.state == "0"]
    df_on = df[df.state == "1"]

    fig, ax = plt.subplots()
    ax.hist(df_on["state_time"], color="r", bins=30, label="active (ON) time in burst")
    ax.hist(df_off["state_time"], color="b", bins=30, label="silent (OFF) time until burst")
    legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

    plt.xlabel("minutes")
    plt.ylabel("occurences")
    plt.title("Distribution of ON and OFF time intervals (minutes)")
    plt.show()
