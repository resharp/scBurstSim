import matplotlib.pyplot as plt
import logging
import os
import sys

from simulator.StrategyGenerator import *
from simulator.StrategyReader import StrategyReader

# script to generate strategies

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"

plot_dir = out_dir + dir_sep + "generate_strategies.plots"
print("creating plot_dir: {plot_dir}".format(plot_dir=plot_dir))
os.makedirs(plot_dir, exist_ok=True)

logger = logging.getLogger(__name__)

logfile_name = out_dir + dir_sep + 'generate_strategies.log'
logging.basicConfig(filename=logfile_name, filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

print("logging in: {file_name}".format(file_name=logfile_name))

# ranges of 10^k are taken from plots generated in read_larsson2019.py
range_k_on_exp = [-3.5, 0.3]
range_k_off_exp = [-2.5, 2.7]
range_k_syn_exp = [-1, 3]
range_k_d_exp = [-1.4, 0.2]

# # generating strategies (per minute)
# range_k_on = [0.005, 0.1]
# range_k_off = [0.005, 0.1]
# range_k_syn = [0.016, 1.6]
# range_k_d = [0.0019, 0.023]

filename = out_dir + dir_sep + "strategies_mixed_new.csv"
sg = StrategyGenerator(range_k_on_exp=range_k_on_exp, range_k_off_exp=range_k_off_exp,
                       range_k_syn_exp=range_k_syn_exp, range_k_d_exp=range_k_d_exp,
                       filename=filename)

k_on_first = False
# sg.generate_and_write_strategies(100, k_on_first)

sg.generate_mixed_set()


def show_distribution(df_strategies, measure):

    units = "minutes"
    if measure == "nr_steady_state_log10":
        units = ""
    plt.title("distribution {measure}({units})".format(measure=measure, units=units))
    plt.hist(df_strategies[measure], bins=30)
    plt.xlabel(measure)
    plot_name = plot_dir + dir_sep + "distrib_{measure}_k_on_first_{k_on_first}.svg".format(
        measure=measure, k_on_first=k_on_first)
    plt.savefig(plot_name)
    plt.close(1)


def read_strategies_and_prepare_data():

    sr = StrategyReader(filename)

    strategies = sr.select_all()

    df_strategies = sr.df_strategies

    df_strategies["burst_time_log10"] = np.log10(1 / df_strategies["k_on"])
    df_strategies["silent_time_log10"] = np.log10(1 / df_strategies["k_off"])
    df_strategies["half_life_log10"] = np.log10(np.log(2)/df_strategies.k_d)
    df_strategies["chance_on"] = df_strategies.k_on / (df_strategies.k_off + df_strategies.k_on)

    df_strategies["chance_on_log10"] = np.log10(df_strategies["chance_on"])
    df_strategies["k_syn_effective"] = df_strategies.k_syn * df_strategies.chance_on
    df_strategies["k_syn_effective_log10"] = np.log10(df_strategies.k_syn_effective)
    df_strategies["nr_steady_state_log10"] = np.log10(df_strategies.k_syn_effective / df_strategies.k_d)

    measures = ["burst_time_log10", "silent_time_log10", "half_life_log10"
                , "chance_on_log10", "chance_on_log10", "k_syn_effective_log10", "nr_steady_state_log10"]

    for measure in measures:
        show_distribution(df_strategies, measure)


# read_strategies_and_prepare_data()

