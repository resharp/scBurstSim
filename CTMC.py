import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# two transitions
l_01 = 0.1
l_10 = 0.1
k_syn = 1.6  # k_syn = synthesis rate = transcription rate
max_minutes = 600
nr_refractions = 1

def new_poisson_arrivals(start_time, state_time) -> list:
    poisson_list = []

    k_syn = 1.6

    last_arrival = 0

    while last_arrival < state_time:
        arrival = - np.log(np.random.rand()) / k_syn
        last_arrival = last_arrival + arrival
        if last_arrival < state_time:
            poisson_list.append([start_time + last_arrival, 1])

    return poisson_list


dtmc_list = []
# second list for doubling all data points in order to show block plot
dtmc_list_plot = []

state = "0"
current_time = 0
poisson_arrivals = []

while current_time < max_minutes:

    state_time = 0
    burst_size = 0
    if state == "0":
        l = l_01
        # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
        # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
        # this is a simple way to simulate multiple refractory states?
        alpha = nr_refractions
        for i in range(alpha):
            state_time = state_time + np.random.exponential(scale=1.0, size=None) / l
    else:
        l = l_10
        state_time = state_time + np.random.exponential(scale=1.0, size=None) / l

        # current_time is the start of active (ON) state
        # state_time is length of burst
        # now create new Poisson arrivals
        new_arrivals = new_poisson_arrivals(current_time, state_time)
        burst_size = len(new_arrivals)
        poisson_arrivals = poisson_arrivals + new_arrivals

    end_time = current_time + state_time
    dtmc_list.append([state, current_time, end_time, state_time, burst_size])
    dtmc_list_plot.append([state, current_time, end_time, state_time, burst_size])

    current_time = end_time
    # switch state
    if state == "0":
        state = "1"
    else:
        state = "0"

    #temp trick: also add new state for showing blocks in plot
    dtmc_list_plot.append([state, current_time, end_time, state_time, burst_size])

df = pd.DataFrame(data=dtmc_list, columns=["state", "begin_time", "end_time", "state_time", "burst_size"])

df_plot = pd.DataFrame(data=dtmc_list_plot, columns=["state", "begin_time", "end_time", "state_time", "burst_size"])

df_poisson_arrivals = pd.DataFrame(poisson_arrivals, columns=["arrival", "count"])
df_poisson_arrivals['cum_count'] = df_poisson_arrivals['count'].cumsum()

debug = "True"

# calculate average burst size
mean_burst_size = df[df.state == "1"].burst_size.mean().round(1)
std_burst_size = df[df.state == "1"].burst_size.std().round(1)
nr_bursts = len(df[df.state == "1"])

plt.title("Burst size: {bs} +/- {std}; # bursts: {nr}".format(
    bs=mean_burst_size, std=std_burst_size, nr=nr_bursts ))
plt.plot(df_poisson_arrivals["arrival"], df_poisson_arrivals["cum_count"])
plt.show()

debug = "True"

# plt.plot(df_plot["begin_time"], df_plot["state"])
# plt.show()

df_off = df[df.state == "0"]
df_on = df[df.state == "1"]

fig, ax = plt.subplots()
ax.hist(df_on["state_time"], color="r", bins=60, label="active (ON)")
ax.hist(df_off["state_time"], color="b", bins=60, label="silent (OFF)")
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.title("Distribution of ON and OFF time intervals")
plt.show()

