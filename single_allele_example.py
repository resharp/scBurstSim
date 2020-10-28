import logging
import os

from simulator.StrategyReader import StrategyReader
from simulator.Transcription import *
from simulator.transcription_plots import *

# script to display single cell single allele example traces

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = ""

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

nr_days = 1
max_minutes = 1440*nr_days  # 24 hours = 1440 minutes

# windows = [[400, 460, 'EU'], [520, 580, '4SU']] # e.g. 120 minutes of EU labeling
start_windows = 600; length_window = 60; between_window = 0
window_eu = [start_windows, start_windows + length_window, 'EU'] # e.g. 120 minutes of EU labeling
window_4su = [start_windows + length_window + between_window,
              start_windows + 2*length_window + between_window, '4SU'] # e.g. 120 minutes of EU labeling
windows = [window_eu, window_4su]
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
fix_time = windows[-1][WINDOW_END] + 0  # fix_time 30 minutes after end of last window


def run_example(params):

    trans = Transcription(params)

    # set complete_trace=True to retrieve the complete trace of transcripts counts (for plotting)
    df_dtmc, dtmc_list = trans.run_bursts(max_minutes, windows, new_dtmc_trace=True, complete_trace=True)

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


logger = logging.getLogger(__name__)
out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
logging.basicConfig(filename=out_dir + dir_sep + 'single_allele_example.log', filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv" )

# see strategy names in data\strategies.csv

# we can select a strategy by name
# params = sr.get(strategy="generated_95")

# or retrieve a random strategy
params = sr.get("one_example")

run_example(params)

# or run an example of all strategies (NB: be sure you have a small strategy file!)
# run_all_strategies()
