# time-dependent solutions for P_00(t), P_01(t), and P_10(t) and P_01(t)
# taken from Anderson 2017 Lecture Notes on Stochastic Processes with Applications in Biology
# P_00(t) is the chance of being in the inactive state at t=t, given being in the inactive state at t=0
# or P_00(t) = P(x_t = 0 | x_0 = 0) with x=0 is being inactive and x=1 is being inactive
# restrictions  P_00(t) + P_01(t) = 1
# and           P_10(t) + P_11(t) = 1
import numpy as np
import matplotlib.pyplot as plt


def P_00(k_on, k_off, t):
    la = k_on
    mu = k_off

    ret_val = mu/(mu + la) + la/(mu + la) * np.exp(-(mu + la)*t)

    return ret_val


def P_10(k_on, k_off, t):
    la = k_on
    mu = k_off

    ret_val = mu/(mu + la) - mu/(mu + la) * np.exp(-(mu + la)*t)

    return ret_val

k_on = 0.01
k_off = 0.04

wins = [j*15 for j in range(0, 9)]

p_00s = []
p_10s = []

for t in wins:

    p_00s.append(P_00(k_on, k_off, t))
    p_10s.append(P_10(k_on, k_off, t))

p_01s = [1 - p for p in p_00s]
p_11s = [1 - p for p in p_10s]


plt.plot(wins, p_00s, 'o-', label="P_00")
plt.plot(wins, p_10s, 'o-', label="P_10")
# plt.plot(wins, p_01s, 'o-', label="P_01")
# plt.plot(wins, p_11s, 'o-', label="P_11")
title = "k_on={k_on};k_off={k_off}. Chances of being in state x_t given state x_0, with e.g. P_00 = P(x_t=0 | x_0=0)".\
    format(k_on=k_on, k_off=k_off)
plt.title(title)
plt.ylim(0, 1)
plt.legend()
plt.xlabel("time in minutes")
plt.ylabel("chance")
plt.show()
plt.close(1)

