import matplotlib.pyplot as plt

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


def plot_events(df_dtmc, df_poisson_arrivals):
    plt.step(df_dtmc["begin_time"], df_dtmc["state"], where="post", color="m")

    arrival_list = df_poisson_arrivals[df_poisson_arrivals.count_s > 0]['arrival'].to_list()
    y_arrivals = [1] * len(arrival_list)

    decay_list = df_poisson_arrivals[df_poisson_arrivals.count_s < 0]['arrival'].to_list()
    y_decays = [0] * len(arrival_list)

    plt.scatter(arrival_list, y_arrivals, color='m', marker="d", s=9)
    plt.scatter(decay_list, y_decays, color='r', marker="o", s=9)


def plot_dynamics(title, df_poisson_arrivals, max_minutes, windows=[], df_labeled_arrivals=[]):

    plt.title(title)

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
