import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# two transitions
l_01 = 1    # make this smaller to decrease burst frequency
l_10 = 1

dtmc_list = []
# second list for doubling all data points in order to show block plot
dtmc_list_plot = []

state = "0"
current_time = 0

for i in range(5000):

    waiting_time = 0
    if state == "0":
        l = l_01
        # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
        # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
        # this is a simple way to simulate multiple refractory states?
        alpha = 6
        for i in range(alpha):
            waiting_time = waiting_time + np.random.exponential(scale=1.0, size=None) / l
    else:
        l = l_10
        waiting_time = waiting_time + np.random.exponential(scale=1.0, size=None) / l

    current_time = current_time + waiting_time

    dtmc_list.append([state, current_time, waiting_time])
    dtmc_list_plot.append([state, current_time, waiting_time])
    # switch state
    if state == "0":
        state = "1"
    else:
        state = "0"

    #temp trick: also add new state for showing blocks in plot
    dtmc_list_plot.append([state, current_time])

df = pd.DataFrame(data=dtmc_list, columns=["state", "time", "waiting_time"])

df_plot = pd.DataFrame(data=dtmc_list_plot, columns=["state", "time", "waiting_time"])

# plt.plot(df_plot["time"], df_plot["state"])
# plt.show()

df_off = df[df.state == "0"]
df_on = df[df.state == "1"]

fig, ax = plt.subplots()
ax.hist(df_on["waiting_time"], color="r", bins=60, label="active (ON)")
ax.hist(df_off["waiting_time"], color="b", bins=60, label="silent (OFF)")
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.title("Distribution of ON and OFF time intervals")
plt.show()


debug = "True"