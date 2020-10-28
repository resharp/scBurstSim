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

# ranges of 10^k are taken from plots generated in read_larsson2019.py
range_k_on_exp = [-3.5, 0.3]
range_k_off_exp = [-2.5, 2.7]
range_k_syn_exp = [-1, 3]
range_k_d_exp = [-1.4, 0.2]

# convert to rates per minute by taking the 10^k and dividing by 60
# range_k_on = [10**k/60 for k in range_k_on_exp]
# range_k_off = [10**k/60 for k in range_k_off_exp]
# range_k_syn = [10**k/60 for k in range_k_syn_exp]
# range_k_d = [10**k/60 for k in range_k_d_exp]

# generating strategies (per minute)
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

    chance_on = params.k_on / (params.k_off + params.k_on)
    k_syn_cor = params.k_syn * chance_on
    ss_mrna = k_syn_cor / params.k_d

    half_life = np.log(2)/params.k_d

    values.append(half_life)

plt.title("distribution half-lives (minutes)")
plt.hist(values, bins=30)
plt.ylim(1)
plt.show()

