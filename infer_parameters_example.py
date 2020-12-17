# TO DO
# three categories of k_syn:
# only change k_on with fixed (k_off, k_syn, k_d)


import os

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

from simulator.Experiment import Experiment
from simulator.Transcription import *
import numpy as np

from utils.utils import round_sig

if os.name == 'nt':
    dir_sep = "\\"
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "infer_parameters_example.plots"
os.makedirs(plot_dir, exist_ok=True)
df_filename = "counts_infer_parameters_example.csv"

k_on = 0.01
k_off = 0.04
k_d = 0.02
k_syn = 0.2
k_eff = 0.1

# window_lengths = [r*15 for r in range(1, 24)]
window_lengths = [15, 30, 45, 60, 120, 180]

k_offs = [k * 0.005 for k in range(1, 6)]   # for some examples in theoretical plots


def p_1(t, k_on, k_off):

    p_on = k_on/(k_on + k_off)
    p_off = k_off/(k_on + k_off)

    p_1 = p_on + p_off * (1 - np.exp(-k_on * t))

    return p_1


# simplified model
def p_1_model(t, k_on, p_on, p_off):

    p_1 = p_on + p_off * (1 - np.exp(-k_on * t))

    return p_1


def nr_molecules_in_window_no_decay(t, k_on, k_off, k_syn, k_eff):

    p_on = k_on/(k_on + k_off)

    nr_mrna = p_on * k_syn * k_eff * t

    return nr_mrna


def plot_theoretical_chance_of_active_state():

    t = np.linspace(0, 400, 100)

    for k_off in k_offs:
        sns.lineplot(x=t, y=p_1(t, k_on, k_off))
    plt.legend(k_offs)

    plt.ylim(0, 1)
    plt.title("k_on={k_on}".format(k_on=k_on))
    plt.ylabel("chance of some active state (any length)")
    plt.xlabel("minutes")

    plt.vlines(x=window_lengths, ymin=0, ymax=1, linestyles='dashed', colors='black')

    plt.savefig(plot_dir + dir_sep + "theoretical_chance_active_{k_on}_{k_off}_{k_syn}.svg".format(
        k_on=k_on, k_off=k_off, k_syn=k_syn))

    plt.close(1)


def plot_production_of_mrna():
    t = np.linspace(0, 400, 100)

    max_y = 0
    for k_off in k_offs:

        y = nr_molecules_in_window_no_decay(t, k_on, k_off, k_syn, k_eff)
        sns.lineplot(x=t, y=y, label="k_off={k_off}".format(k_off=k_off))

    plt.legend()

    plt.title("k_on={k_on}".format(k_on=k_on))
    plt.ylabel("average nr of molecules produced")
    plt.xlabel("minutes")

    plt.vlines(x=window_lengths, ymin=0, ymax=max(y), linestyles='dashed', colors='black')

    plt.savefig(plot_dir + dir_sep + "theoretical_production_mrna_{k_on}_{k_off}_{k_syn}.svg".format(
        k_on=k_on, k_off=k_off, k_syn=k_syn))
    plt.close(1)


def get_windows_and_fix_time(length_window=60, gap=0):

    start_windows = 600
    window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
    window_4su = [start_windows + length_window + gap,
                  start_windows + 2 * length_window + gap, '4SU'] # e.g. 120 minutes of EU labeling
    windows = [window_eu, window_4su]
    fix_time = windows[-1][WINDOW_END] + 0  # fix_time 0 minutes after end of last window

    return windows, fix_time


def run_active_state_simulations(nr_runs):

    l_counts = []

    for w in window_lengths:
        nr_runs_active = 0
        nr_real_label = 0
        nr_signal_label = 0

        windows, fix_time = get_windows_and_fix_time(length_window=w, gap=0)

        params = TranscriptParams(k_on=k_on, k_off=k_off, nr_refractions=1,
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

            # TODO: sampling should be done differently
            # here we are taking a fixed percentage
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

    df_counts = pd.DataFrame(l_counts, columns=["window", "active", "real", "signal"])
    df_counts.to_csv(out_dir + dir_sep + df_filename, sep=';', index=False)

    return df_counts


def plot_chance_of_switching_to_active_state(df_counts, nr_runs):

    # we want to convert to
    plt.plot(df_counts.window, df_counts.active/nr_runs, label='with active state')
    plt.plot(df_counts.window, df_counts.real/nr_runs, label='with real counts')
    plt.plot(df_counts.window, df_counts.signal/nr_runs, label='with detected counts')

    plt.plot(df_counts.window, df_counts.theoretical, color="red", label="theoretical")

    plt.xlim(0, max(window_lengths) + 15)
    # plt.ylim(0, 1)
    plt.xlabel("window size (minutes)")
    plt.ylabel("nr of runs")
    plt.legend()
    plt.savefig(plot_dir + dir_sep + "counts_{k_on}_{k_off}_{k_syn}.svg".format(
        k_on=k_on, k_off=k_off, k_syn=k_syn))
    plt.close(1)


run_sim = False
nr_runs = 500
if run_sim:
    df_counts = run_active_state_simulations(nr_runs)
else:
    df_counts = pd.read_csv(out_dir + dir_sep + df_filename, sep=';')


def fit_to_model_p1():

    expected = (0.1, 0.5, 0.5)

    # divide by nr_runs for getting chance
    popt, pcov = curve_fit(p_1_model, df_counts.window, df_counts.active / nr_runs, expected)
    popt_active = popt
    error_k_on_active = abs(popt_active[0] / k_on - 1) * 100

    popt, pcov = curve_fit(p_1_model, df_counts.window, df_counts.real / nr_runs, expected)
    popt_real = popt
    error_k_on_real = abs(popt_real[0] / k_on - 1) * 100

    popt, pcov = curve_fit(p_1_model, df_counts.window, df_counts.signal / nr_runs, expected)
    popt_signal = popt
    error_k_on_signal = abs(popt_signal[0] / k_on - 1) * 100

    print("fitting to hidden state:   k_on={k_on}; error={error}%".format(
        k_on=round_sig(popt_active[0], 4), error=round_sig(error_k_on_active, 3)))
    print("fitting to real counts:    k_on={k_on}; error={error}%".format(
        k_on=round_sig(popt_real[0], 4), error=round_sig(error_k_on_real, 3)))
    print("fitting to sampled counts: k_on={k_on}; error={error}%".format(
        k_on=round_sig(popt_signal[0]), error=round_sig(error_k_on_signal, 3)))


plot_theoretical_chance_of_active_state()
plot_production_of_mrna()

df_counts["theoretical"] = p_1(df_counts["window"], k_on, k_off)

plot_chance_of_switching_to_active_state(df_counts, nr_runs)

fit_to_model_p1()
