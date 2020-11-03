import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns

from simulator.Experiment import Experiment
from simulator.Transcription import *

k_on = 0.005

k_offs = [k * 0.005 for k in range(1, 6)]
k_d = 0.002
k_syn = 0.032
k_eff = 1

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


def run_active_state_simulations():

    nr_runs = 1000
    k_off = 0.02
    nr_active = []
    for w in window_lengths:
        nr_runs_active = 0
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

            # example of calculating percentage active
            label = "4SU"
            perc = Experiment.perc_active_state(windows, df_dtmc, label)
            # print("Percentage active state: {perc}".format(perc=perc))
            if perc > 0:
                nr_runs_active = nr_runs_active + 1
        print("{label} window contains {nr_runs_active} runs with active state(s) for k_off {k_off} and window {window}".
              format(label=label, k_off=k_off, window=w, nr_runs_active=nr_runs_active))
        nr_active.append(nr_runs_active)

    plt.plot(window_lengths, nr_active)
    plt.xlim(0, max(window_lengths) + 15)
    plt.ylim(0, nr_runs)
    plt.xlabel("window size (minutes)")
    plt.ylabel("nr of runs with active state")
    plt.show()


run_active_state_simulations()
