# non stationary solution for P(n, t) taken from Dattani et al, J R Soc Interface 2016
# https://doi.org/10.1098/rsif.2016.0833

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.poch.html
# from scipy.special import gamma, hyp1f1, poch
import scipy.special as sc
import numpy as np
import math
import matplotlib.pyplot as plt
import logging
import sys
import os

from utils.utils import round_sig

fac = np.math.factorial  # in Python even a function is an object

out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs\non_stationary.plots"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"
plot_dir = out_dir + dir_sep + "non_stationary.plots"
os.makedirs(plot_dir, exist_ok=True)
logging.basicConfig(filename=out_dir + dir_sep + 'nonstationary.log', filemode='w',
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    level=logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def p_stationary(n, k_on, k_off, k_syn, k_d):

    k_on = k_on / k_d
    k_off = k_off / k_d
    k_syn = k_syn / k_d

    part1 = k_syn**n/fac(n)
    part2 = sc.poch(k_on, n)/sc.poch(k_on + k_off, n)
    part3 = sc.hyp1f1(k_on + n, k_on + k_off + n, - k_syn)

    ret_val = part1 * part2 * part3

    return ret_val


def create_distribution(k_on, k_off, k_syn, k_d):
    y_list = []
    x_list = []
    x = 0
    y = 1

    max = int(1.3 * k_syn / k_d)
    for n in range(0, max):   # this alternative is for debugging N=0 and/or N=1
    # while y > 1e-6:
        y = p_stationary(x, k_on, k_off, k_syn, k_d)

        x_list.append(x)
        y_list.append(y)
        x = x + 1
    return x_list, y_list


def plot_distribution(k_on, k_off, k_syn, k_d):

    x_list, y_list = create_distribution(k_on, k_off, k_syn, k_d)

    plt.step(x_list, y_list, where="post")
    plt.show()
    plt.close(1)


# t is times in minutes
def p_non_stationary(n, t, k_on, k_off, k_syn, k_d):

    k_on = k_on / k_d
    k_off = k_off / k_d
    mu = k_syn / k_d    # change to symbol mu for better comparison with paper

    factor_a = mu**n / fac(n)
    # logging.info("factor_a: {}".format(factor_a))
    sum_a = 0

    for r in range(0, n+1):

        logging.debug("--- start part_a for r={r}".format(r=r))
        part_a1 = math.comb(n, r)
        logging.debug("math.comb(n, r)             : {}".format(part_a1))
        part_a2 = (-1)**r * sc.poch(-k_on, r) * sc.poch(k_on, n-r) * np.exp(-r*t)
        logging.debug("(-1)**r * sc.poch(-k_on, r) * sc.poch(k_on, n-r) * np.exp(-r*t): {}".format(part_a2))
        # there is a problem with this part if k_on + k_off = 1 and r is non-zero
        part_a3 = 1/sc.poch(1 - k_on - k_off, r)
        logging.debug("1/sc.poch(1 - k_on - k_off, r): {}".format(part_a3))
        part_a4 = 1/sc.poch(k_on + k_off, n-r)
        logging.debug("1/sc.poch(k_on + k_off, n-r): {}".format(part_a4))

        part_a5 = sc.hyp1f1(-k_on + r, 1 - k_on - k_off + r, mu * np.exp(-t))
        logging.debug("sc.hyp1f1(-k_on + r, 1 - k_on - k_off + r, mu * np.exp(-t): {}".format(part_a5))
        part_a6 = sc.hyp1f1(k_on + n - r, k_on + k_off + n - r, -mu)
        logging.debug(part_a6)

        part_a = part_a1 * part_a2 * part_a3 * part_a4 * part_a5 * part_a6
        logging.debug("--- part_a: {part_a} for r={r}".format(part_a=part_a, r=r))

        sum_a = sum_a + part_a

    sum_a = sum_a * factor_a

    fac_b1 = k_on * mu**(n + 1) * np.exp(-(k_on + k_off) * t)
    fac_b2 = 1 / (
            (k_on + k_off) *
            (1 - k_on - k_off) *
            fac(n)
            )
    factor_b = fac_b1 * fac_b2
    sum_b = 0
    for r in range(0, n+1):

        logging.debug("--- start part_b for r={r}".format(r=r))
        part_b = part_of_sum(k_on, k_off, mu, n, r, t)
        logging.debug("--- part_b: {part_b} for r={r}".format(part_b=part_b, r=r))
        sum_b = sum_b + part_b

    sum_b = factor_b * sum_b

    fac_c1 = k_on * mu**n * np.exp(-(k_on + k_off) * t)
    fac_c2 = fac_b2

    factor_c = fac_c1 * fac_c2
    sum_c = 0
    for r in range(0, n):

        # the sum part for C seems to be identical to the sum part of B
        # with n replaced by n-1
        # except for the n! in the denominator
        # is it a typo? or is it to be expected?

        logging.debug("--- start part_c for r={r}".format(r=r))
        part_c = part_of_sum(k_on, k_off, mu, n-1, r, t)

        logging.debug("--- part_c: {part_c} for r={r}".format(part_c=part_c, r=r))
        sum_c = sum_c + part_c

    # NB: minus sign in front of the third term!
    sum_c = - factor_c * sum_c

    total_sum = sum_a + sum_b + sum_c

    logging.debug("sum_a for n={n}: {sum_a}".format(n=n, sum_a=sum_a))
    logging.debug("sum_b for n={n}: {sum_b}".format(n=n, sum_b=sum_b))
    logging.debug("sum_c for n={n}: {sum_c}".format(n=n, sum_c=sum_c))
    logging.info("total_sum for n={n}: {total_sum}".format(n=n, total_sum=total_sum))

    return factor_a, total_sum


def part_of_sum(k_on, k_off, mu, n, r, t):

    part_b1 = math.comb(n, r)
    logging.debug(part_b1)
    part_b2 = (-1) ** r * sc.poch(-k_off, r) * sc.poch(1 - k_off, n - r) * np.exp(-r * t)
    logging.debug(part_b2)
    part_b3 = 1 / sc.poch(1 + k_on + k_off, r)
    logging.debug(part_b3)
    part_b4 = 1 / sc.poch(2 - k_on - k_off, n - r)
    logging.debug(part_b4)
    part_b5 = sc.hyp1f1(k_off + r,
                        1 + k_on + k_off + r,
                        mu * np.exp(-t))
    logging.debug(part_b5)

    # NB: some terms in the Pochhammer's reappear in 1F1
    # and may cancel each other and this may avoid dividing by zero?
    part_b6 = sc.hyp1f1(1 - k_off + n - r,
                        2 - k_on - k_off + n - r,
                        -mu)
    logging.debug(part_b6)
    part_b = part_b1 * part_b2 * part_b3 * part_b4 * part_b5 * part_b6

    logging.debug(part_b)

    return part_b


def create_non_stationary_distribution(time, k_on, k_off, k_syn, k_d):

    x_list = []
    y_list = []

    # e.g. for n in range(0, 120):
    total_p = 0
    n = 0
    p_ns = 1  # initialize at large value
    # while p_ns > 1e-6:

    max = int(1.3 * k_syn / k_d)
    for n in range(0, max):   # this alternative is for debugging N=0 and/or N=1
        logging.info("*********************")
        logging.info("**** start for n={}".format(n))
        logging.info("*********************")
        factor_a, p_ns = p_non_stationary(n, time, k_on, k_off, k_syn, k_d)

        x_list.append(n)
        y_list.append(p_ns)
        # y_list.append(factor_a)   # for debugging factor_a
        total_p = total_p + p_ns
        n = n + 1

    logging.info("*********************")
    logging.info("total_p for all n: {}".format(total_p))
    logging.info("---")

    parameters = "t={time}; k_on={k_on}; k_off={k_off}; k_syn={k_syn}; k_d={k_d}".format(
        time=time, k_on=k_on, k_off=k_off, k_syn=k_syn, k_d=k_d
    )

    logging.info(parameters)

    return x_list, y_list


k_on = 0.151
k_off = 0.151
k_syn = 2
k_d = 0.02

# plot_distribution("test", k_on, k_off, k_syn, k_d)

times = [2.3, 2.5] + list(range(3, 6, 1))
# times = [0.8] + list(range(1, 4, 1))
# times = [0.25, 0.3, 0.4]
# times = list(range(4, 6))
# times = [0]
# times = [2, 3]
times = [1, 2, 3]


show_plot = True
show_stationary = True
for time in times:

    x_list, y_list = create_non_stationary_distribution(time, k_on, k_off, k_syn, k_d)

    sum_p = round_sig(sum(y_list), 4)

    if show_plot:
        plt.step(x_list, y_list, where="post", label="t={time}; sum_p={sum_p}".format(time=time, sum_p=sum_p))

x_list_stat, y_list_stat = create_distribution(k_on, k_off, k_syn, k_d)

if show_plot and show_stationary:
    sum_p = round_sig(sum(y_list_stat), 4)
    plt.step(x_list_stat, y_list_stat, where="post",
             label="stationary; sum_p={sum_p}".format(sum_p=sum_p)
             , color="black")

if show_plot:

    plt.legend(title='Time:')
    parameters = "k_on={k_on}; k_off={k_off}; k_syn={k_syn}; k_d={k_d}".format(
        k_on=k_on, k_off=k_off, k_syn=k_syn, k_d=k_d
    )
    plt.title(parameters)
    plt.show()
    plt.close(1)

