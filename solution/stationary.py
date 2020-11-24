# stationary solution for P(n) taken from Raj, PLOS Biology 2006
# lambda = k_on
# gamma = k_off
# mu = k_syn
# delta = k_d

from scipy.special import gamma, hyp1f1
from utils.utils import round_sig
import matplotlib.pyplot as plt
import os

# what would we like to put in here?
from simulator.StrategyReader import StrategyReader

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"
out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
plot_dir = out_dir + dir_sep + "stationary.plots"
os.makedirs(plot_dir, exist_ok=True)


def p_stationary(n, k_on, k_off, k_syn, k_d):

    l = k_on/k_d
    g = k_off/k_d
    mu = k_syn/k_d

    part1 = gamma(l + n)/gamma(n + 1)
    part1 = part1 / gamma(l + g + n)

    part2 = gamma(l + g) / gamma(l)
    part3 = mu**n
    part4 = hyp1f1(l+n, l+g+n, -mu)

    ret_val = part1 * part2 * part3 * part4

    return ret_val


# approximation for mean based on high k_syn as compared to other parameters
def mean_with_high_k_syn(k_on, k_off, k_syn, k_d):

    ret_val = (k_syn / k_d) * k_on / (k_on + k_off)
    return ret_val


def plot_distribution(strategy, k_on, k_off, k_syn, k_d):

    y_list = []
    x_list = []
    x = 0
    y = 1
    while y > 1e-5:
        y = p_stationary(x, k_on, k_off, k_syn, k_d)

        x_list.append(x)
        y_list.append(y)
        x = x + 1

    weighted = [(x*y) for (x, y) in zip(x_list, y_list)]

    mean_approx = mean_with_high_k_syn(k_on, k_off, k_syn, k_d)
    mean_approx = round_sig(mean_approx, 4)
    mean_real = round_sig(sum(weighted), 4)

    plt.step(x_list, y_list, where="post")
    plt.title("PMF for strategy '{strategy}'; mean approx={mean_approx}; real mean={mean_real}".format(
        mean_approx=mean_approx, mean_real=mean_real, strategy=strategy
    ))
    plot_name = plot_dir + dir_sep + "stationary_{strategy}.svg".format(strategy=strategy)
    plt.savefig(plot_name)
    plt.close(1)


k_on = 0.01
k_off = 0.01
k_syn = 2
k_d = 0.02

plot_distribution("example", k_on, k_off, k_syn, k_d)


def plot_strategies_from_file():

    sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
    sr.read_strategies()

    df_strategies = sr.df_strategies

    strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]

    for strategy in strategies:
        df = df_strategies[df_strategies.name == strategy]
        plot_distribution(strategy, df.k_on.item(), df.k_off.item(), df.k_syn.item(), df.k_d.item())


plot_strategies_from_file()
