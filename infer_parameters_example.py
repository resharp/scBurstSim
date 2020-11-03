import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns

from simulator.Experiment import Experiment
from simulator.Transcription import *

if os.name == 'nt':
    dir_sep = "\\"
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "infer_parameters_example.plots"
os.makedirs(plot_dir, exist_ok=True)
df_filename = "counts_infer_parameters_example.csv"

k_on = 0.005

k_offs = [k * 0.005 for k in range(1, 6)]
k_d = 0.002
k_syn = 1.6
k_eff = 0.1

window_lengths = [r*15 for r in range(1, 24)]


def p_1(t, k_on, k_off):

    p_on = k_on/(k_on + k_off)
    p_off = k_off/(k_on + k_off)

    p_1 = p_on + p_off * (1 - np.exp(-k_on * t))

    return p_1


def nr_molecules(t, k_on, k_off, k_syn, k_d, k_eff):

    p_on = k_on/(k_on + k_off)

    mol = p_on * k_syn * k_eff * t

    return mol


def show_chance_of_active_state():

    t = np.linspace(0, 900, 100)

    for k_off in k_offs:
        # sns.lineplot(x=t, y=p_1(t, k_on, k_off))
        sns.lineplot(x=t, y=p_1(t, k_on, k_off))
    plt.legend(k_offs)

    plt.ylim(0 ,1)
    plt.title("k_on={k_on};k_offs={k_offs}".format(k_on=k_on, k_offs=k_offs))
    plt.ylabel("chance of some active state (any length)")
    plt.xlabel("minutes")

    plt.vlines(x=window_lengths, ymin=0, ymax=1, linestyles='dashed', colors='black')

    plt.show()
    plt.close(1)


def show_production_of_mrna():
    t = np.linspace(0, 180, 100)

    max_y = 0
    for k_off in k_offs:
        # sns.lineplot(x=t, y=p_1(t, k_on, k_off))

        y = nr_molecules(t, k_on, k_off, k_syn, k_d, k_eff)
        sns.lineplot(x=t, y=y)

    plt.legend(k_offs)

    plt.title("k_on={k_on};k_offs={k_offs}".format(k_on=k_on, k_offs=k_offs))
    plt.ylabel("nr of molecules produced")
    plt.xlabel("minutes")

    plt.vlines(x=window_lengths, ymin=0, ymax=max(y), linestyles='dashed', colors='black')

    plt.show()
    plt.close(1)


def get_windows_and_fix_time(length_window=60, gap=0):

    start_windows = 600
    window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
    window_4su = [start_windows + length_window + gap,
                  start_windows + 2 * length_window + gap, '4SU'] # e.g. 120 minutes of EU labeling
    windows = [window_eu, window_4su]
    fix_time = windows[-1][WINDOW_END] + 0  # fix_time 0 minutes after end of last window

    return windows, fix_time


# this is the theoretical possibility
# show_chance_of_active_state()
# show_production_of_mrna()


def run_active_state_simulations(nr_runs):

    k_off = 0.02
    l_counts = []

    for w in window_lengths:
        nr_runs_active = 0
        nr_real_label = 0
        nr_signal_label = 0

        windows, fix_time = get_windows_and_fix_time(length_window=w, gap=0)

        params = TranscriptParams(k_on=k_on, k_off=k_off, nr_refractions=2,
                                  tm_id=np.nan,
                                  k_syn=k_syn, k_d=k_d,
                                  coord_group=0,
                                  name="test")

        trans = Transcription(params)

        # set complete_trace=True to retrieve the complete trace of transcripts counts (for plotting)
        for run in range(0, nr_runs):

            df_dtmc, dtmc_list = trans.run_bursts(fix_time, windows, new_dtmc_trace=True, complete_trace=True)

            label = "4SU"
            df_transcripts = trans.df_transcripts

            df_labeled_transcripts = df_transcripts[df_transcripts.label == label]
            if len(df_labeled_transcripts) > 0:
                nr_real_label = nr_real_label + 1

            len_sample = int(k_eff * len(df_labeled_transcripts))
            df_sampled = df_transcripts.sample(len_sample, replace=False)
            if len(df_sampled) > 0:
                nr_signal_label = nr_signal_label + 1

            # example of calculating percentage active
            perc = Experiment.perc_active_state(windows, df_dtmc, label)
            # print("Percentage active state: {perc}".format(perc=perc))
            if perc > 0:
                nr_runs_active = nr_runs_active + 1
        print("{label} window contains {nr_runs_active} runs with active state(s) for k_off {k_off} and window {window}".
              format(label=label, k_off=k_off, window=w, nr_runs_active=nr_runs_active))
        l_counts.append([w, nr_runs_active, nr_real_label, nr_signal_label])

    df_counts = pd.DataFrame(l_counts, columns=['window', 'active', 'real', 'signal'])
    df_counts.to_csv(out_dir + dir_sep + df_filename, sep=';', index=False)

    return df_counts


def show_plot(df_counts):

    plt.plot(df_counts.window, df_counts.active)
    plt.plot(df_counts.window, df_counts.real)
    plt.plot(df_counts.window, df_counts.signal)
    plt.xlim(0, max(window_lengths) + 15)
    plt.ylim(0, nr_runs)
    plt.xlabel("window size (minutes)")
    plt.ylabel("nr of runs with active state")
    plt.savefig(plot_dir + dir_sep + "counts.svg")
    # plt.show()
    plt.close(1)


run_sim = False
nr_runs = 1000
if run_sim:
    df_counts = run_active_state_simulations(nr_runs)
else:
    df_counts = pd.read_csv(out_dir + dir_sep + df_filename, sep=';')

show_plot(df_counts)