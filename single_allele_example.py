import logging
import os

from simulator.Experiment import Experiment
from simulator.StrategyReader import StrategyReader
from simulator.Transcription import *
from simulator.transcription_plots import *

from solution.stationary import create_distribution
from utils.utils import round_sig
# goal:
# 1. plot single cell single allele example traces
# 2. quickly create stationary distribution for a single allele example (by sampling snapshots from single trace)
#
# single_allele_example.py uses the Transcription class that is also used by Experiment.py
if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = ""

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

plot_dir = out_dir + dir_sep + "single_allele_example.plots"
os.makedirs(plot_dir, exist_ok=True)

nr_days = 1
max_minutes = 1440*nr_days  # 24 hours = 1440 minutes

logger = logging.getLogger(__name__)
out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
logging.basicConfig(filename=out_dir + dir_sep + 'single_allele_example.log', filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv" )


# windows = [[400, 460, 'EU'], [520, 580, '4SU']] # e.g. 120 minutes of EU labeling
def get_windows_and_fix_time(length_window=60, gap=0):

    start_windows = 600
    window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
    window_4su = [start_windows + length_window + gap,
                  start_windows + 2 * length_window + gap, '4SU'] # e.g. 120 minutes of EU labeling
    windows = [window_eu, window_4su]
    fix_time = windows[-1][WINDOW_END] + 0  # fix_time 0 minutes after end of last window

    return windows, fix_time


windows, fix_time = get_windows_and_fix_time(length_window=60, gap=0)


def run_example(params):

    trans = Transcription(params)

    # set complete_trace=True to retrieve the complete trace of transcripts counts (for plotting)
    df_dtmc, dtmc_list = trans.run_bursts(max_minutes, windows, new_dtmc_trace=True, complete_trace=True)

    # example of calculating percentage active
    label = "4SU"
    perc = Experiment.perc_active_state(windows, df_dtmc, label)
    print("Percentage active state: {perc}".format(perc=perc))
    if perc > 0:
        print("{label} window contains active state(s)".format(label=label))

    # only set when complete_trace = True:
    # df_events, df_unlabeled_events and dfs_labeled_events
    df_events = trans.df_events
    df_unlabeled_events = trans.df_unlabeled_events
    dfs_labeled_events = trans.dfs_labeled_events

    # calculate average burst size
    nr_bursts = len(df_dtmc[df_dtmc.state == "1"])
    burst_frequency = round(nr_bursts / max_minutes, 3)

    title = "strategy={name}; k_on={k_on}; k_off={k_off};k_syn={k_syn}; k_d={k_d} -> " \
            "burst freq: {freq}".format(name=params.name, freq=burst_frequency,
                                        k_on=params.k_on, k_off=params.k_off, k_syn=params.k_syn, k_d=params.k_d)
    plot_events(df_dtmc, df_events)
    plot_dynamics(title=title, df_events=df_unlabeled_events, fix_time=fix_time, max_minutes=max_minutes,
                  windows=windows, dfs_labeled_events=dfs_labeled_events)

    # show distribution of time in burst and of time until burst (in minutes)
    # plot_waiting_time_distribution(df_dtmc)


def run_all_strategies():
    global params
    for params in sr.select_all():
        run_example(params)


# alternative way to create a distribution from a single trace
def run_distribution(params, interval, nr_snapshots):

    max_minutes = interval * nr_snapshots

    trans = Transcription(params)

    # set complete_trace=True to retrieve the complete trace of transcripts counts (for plotting)
    # here we set windows = [], because we will create a distribution from snapshots from a single trace
    df_dtmc, dtmc_list = trans.run_bursts(max_minutes, windows=[], new_dtmc_trace=True, complete_trace=True)

    df_unlabeled_events = trans.df_unlabeled_events

    interval = 100  # minutes between snapshots
    fix_times = [i * interval for i in range(0, int(max_minutes/interval))]
    counts = []
    for fix_time in fix_times:
        df = df_unlabeled_events[df_unlabeled_events.arrival < fix_time]
        if len(df) > 0:
            last = df.tail(1)
            count = last.cum_count.item()
        else:
            count = 0
        counts.append(count)

    df = pd.DataFrame(data=counts, columns=["mrna"]).reset_index()
    df.rename(columns={'index': 'snapshot'}, inplace=True)
    mean_mrna = round_sig(df.mrna.mean(), 4)

    df_distribution = df.groupby('mrna')['snapshot'].count().to_frame().reset_index()

    df_distribution['mrna'] = df_distribution.mrna.astype(int)
    max_count = df_distribution.mrna.max()
    df_distribution = df_distribution.set_index('mrna').reindex(range(0, max_count + 1)).fillna(0).reset_index()

    nr_snapshots = max_minutes / interval
    df_distribution["chance"] = df_distribution.snapshot / nr_snapshots

    plt.title("mean nr of RNA for strategy {strategy}: {mean}".format(mean=mean_mrna, strategy=strategy))

    # distribution from simulation
    plt.step(df_distribution.mrna, df_distribution.chance, where="post")

    # theoretical distribution
    x_list, y_list = create_distribution(params.k_on, params.k_off, params.k_syn, params.k_d)
    plt.step(x_list, y_list, where="post", color="red")

    plot_name = plot_dir + dir_sep + "distribution_{strategy}.svg".format(strategy=strategy)
    plt.savefig(plot_name)
    plt.close(1)


# we can select a strategy by name
# see strategy names in data\strategies.csv
strategy = "second_example"
params = sr.get(strategy=strategy)

# or retrieve a random strategy
# params = sr.get_random()

run_example(params)

# or run an example of all strategies (NB: be sure you have a small strategy file!)
# run_all_strategies()

interval = 100
nr_snapshots = 300

# run_distribution(params, interval, nr_snapshots)

