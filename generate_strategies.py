import matplotlib.pyplot as plt
import logging
import os

from simulator.StrategyGenerator import *
from simulator.StrategyReader import StrategyReader

# script to generate strategies

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = ""

logger = logging.getLogger(__name__)
logging.basicConfig(filename=out_dir + dir_sep + 'generate_strategies.log', filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)

# generating strategies

range_k_on = [0.005, 0.1]
range_k_off = [0.005, 0.1]
range_k_syn = [0.016, 1.6]
range_k_d = [0.0019, 0.023]

filename = out_dir + dir_sep + "strategies_generated.csv"
sg = StrategyGenerator(range_k_on=range_k_on, range_k_off=range_k_off, range_k_syn=range_k_syn, range_k_d=range_k_d,
                       filename=filename)

sg.generate_and_write_strategies(100)

sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )

strategies = sr.select_all()

values = []
for params in strategies:

    chance_on = params.k_01 / (params.k_10 + params.k_01)
    k_syn_cor = params.k_syn * chance_on
    ss_mrna = k_syn_cor / params.k_d

    half_life = np.log(2)/params.k_d

    values.append(half_life)

plt.title("distribution half-lives (minutes)")
plt.hist(values, bins=30)
plt.ylim(1)
plt.show()

