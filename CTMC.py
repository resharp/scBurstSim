import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# two transitions
l_01 = 0.02
l_10 = 0.02
k_syn = 0.16  # k_syn = synthesis rate = transcription rate
k_d = 0.01   # k_d = decay rate
max_minutes = 1440 # 24 hours = 1440 minutes
nr_refractions = 1


def new_poisson_arrivals(start_time, interval) -> list:
    poisson_list = []

    last_arrival = 0

    while last_arrival < interval:
        arrival = np.random.exponential(scale=1.0, size=None) / k_syn
        last_arrival = last_arrival + arrival

        # we immediately determine the decay time once a transcript comes into existence
        decay = np.random.exponential(scale=1.0, size=None) / k_d
        decay_time = start_time + last_arrival + decay

        if last_arrival < interval:
            poisson_list.append([start_time + last_arrival, 1, decay_time, -1])

    return poisson_list


dtmc_list = []

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

    current_time = end_time
    # switch state
    if state == "0":
        state = "1"
    else:
        state = "0"

df = pd.DataFrame(data=dtmc_list, columns=["state", "begin_time", "end_time", "state_time", "burst_size"])

df_poisson_arrivals = pd.DataFrame(poisson_arrivals, columns=["arrival", "count_s", "decay", "count_d"])

# we will now put the arrivals and decays in one table and sort by time
df_decays = df_poisson_arrivals[['decay', "count_d"]].\
    rename(columns={'decay': 'arrival', 'count_d': 'count_s'})

df_poisson_arrivals = df_poisson_arrivals[["arrival", "count_s"]]
df_poisson_arrivals = df_poisson_arrivals.append(df_decays).sort_values(by="arrival")

df_poisson_arrivals['cum_count'] = df_poisson_arrivals['count_s'].cumsum()

debug = "True"

# calculate average burst size
mean_burst_size = df[df.state == "1"].burst_size.mean().round(1)
std_burst_size = df[df.state == "1"].burst_size.std().round(1)
nr_bursts = len(df[df.state == "1"])
burst_frequency = round(nr_bursts/max_minutes, 3)

plt.title("Burst size: {bs} +/- {std}; # burst frequency: {freq}".format(
    bs=mean_burst_size, std=std_burst_size, freq=burst_frequency))
plt.step(df_poisson_arrivals["arrival"], df_poisson_arrivals["cum_count"])
plt.xlim(0, max_minutes)
plt.xlabel("minutes")
plt.ylabel("nr of transcripts")
plt.show()

debug = "True"

# plt.plot(df_plot["begin_time"], df_plot["state"])
# plt.show()

# df_off = df[df.state == "0"]
# df_on = df[df.state == "1"]
#
# fig, ax = plt.subplots()
# ax.hist(df_on["state_time"], color="r", bins=60, label="active (ON)")
# ax.hist(df_off["state_time"], color="b", bins=60, label="silent (OFF)")
# legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
#
# plt.title("Distribution of ON and OFF time intervals")
# plt.show()

