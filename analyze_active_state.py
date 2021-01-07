import os
import pandas as pd
import matplotlib.pyplot as plt
from simulator.StrategyReader import StrategyReader
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import beta

if os.name == 'nt':
    dir_sep = "\\"
    # TODO: set your own out directory
    out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"
else:
    dir_sep = "/"
    out_dir = "sc_runs"
plot_dir = out_dir + dir_sep + "analyze_active_state.plots"
os.makedirs(plot_dir, exist_ok=True)

in_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"

len_win = 180  # 60, 100, 180
gap = 0
label_1 = "EU"
label_2 = "4SU"

# sr = StrategyReader(out_dir + dir_sep + "strategies_generated.csv" )
sr = StrategyReader(in_dir + dir_sep + "strategies.csv")
sr.read_strategies()
df_strategies = sr.df_strategies


filename_counts = out_dir + dir_sep + "df_counts_W{len_win}_G{gap}.csv".format(
    len_win=len_win, gap=gap)
df_counts = pd.read_csv(filename_counts, sep=';')


def get_parameters(strategy):
    df_strategy = df_strategies[df_strategies.name == strategy]

    return df_strategy.k_on.item(), df_strategy.k_off.item()


def create_distribution(df_counts, strategy):
    df_counts_label = df_counts[(df_counts.label == label_2) & (df_counts.strategy == strategy)]

    # convert to density for perc (=percentage active state)

    bins = np.linspace(0, 1, 100)

    pdf = np.histogram(df_counts_label.perc, bins=bins)[0]

    bins = np.delete(bins, 0)  # remove first
    sum_pdf = sum(pdf)
    norm_pdf = (pdf / sum_pdf) * len(bins)

    return bins, norm_pdf


def plot_distribution_active_state(strategies):

    for strategy in strategies:

        bins, norm_pdf = create_distribution(df_counts, strategy)

        # plt.step(bins, norm_pdf, where="post")
        plt.plot(bins, norm_pdf)

        plt.ylim(0, None)
        plt.xlabel("percentage active state")
        plt.ylabel("chance")

        k_on, k_off = get_parameters(strategy)

        plt.title("fraction active period ({strategy});k_on={k_on},k_off={k_off} (win={len_win})".format(
            strategy=strategy, k_on=k_on, k_off=k_off, len_win=len_win))
        plot_name = plot_dir + dir_sep + "percentage_state_distribution_{strategy}_{len_win}.svg".\
            format(strategy=strategy, len_win=len_win)
        plt.savefig(plot_name)
        plt.close(1)


def beta_model(x, a, b):
    return beta.pdf(x, a, b)


# fit to beta distribution
# TODO: find a way to avoid local optima, the result is very dependent on "expected" (initial conditions of
# fitting procedure, e.g. for t=180 and strategy="third_example" you only find a good fit if you expect
# both a > 1 and b > 1; the peak at perc=0 and perc=1 are not fitted)
def fit_distribution_active_state_to_beta_distribution(strategy):

    bins, norm_pdf = create_distribution(df_counts, strategy)

    expected = (1.01, 1.01)

    popt, pcov = curve_fit(beta.pdf, bins, norm_pdf, expected)

    a = popt[0]
    b = popt[1]

    # a = 0.5
    # b = 0.5

    # x = np.linspace(beta.ppf(0.01, a, b), beta.ppf(0.99, a, b), 100)
    x = np.linspace(0, 1, 100)
    plt.plot(x, beta.pdf(x, a, b), 'r-', lw=2, alpha=0.6, label='beta pdf')
    plt.plot(bins, norm_pdf, 'b-', lw=2, alpha=0.6, label='norm pdf')
    plt.title("estimated a: {a}, b: {b}".format(a=a, b=b))
    plt.show()
    plt.close(1)


# build intuition: create distributions for percentage of active state ("active state fraction distribution")
strategies = ["first_example", "second_example", "third_example", "bimodal", "powerlaw"]
plot_distribution_active_state(strategies)

# explore: try to fit to beta distribution
strategy = "third_example"
fit_distribution_active_state_to_beta_distribution(strategy)


