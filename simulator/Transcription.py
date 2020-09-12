import numpy as np
import pandas as pd

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class Transcription:
    l_01 = 0
    l_10 = 0
    k_syn = 0  # k_syn = synthesis rate = transcription rate
    nr_refractions = 0
    # decay is strictly not part of transcription but we include it in the model
    k_d = 0.01  # k_d = decay rate
    state = "0"

    df_dtmc = None
    df_poisson_arrivals = None

    def __init__(self, l_01, l_10, k_syn, nr_refractions, k_d):

        self.l_01 = l_01
        self.l_10 = l_10
        self.k_syn = k_syn
        self.k_d = k_d
        self.nr_refractions = nr_refractions

        self.state = "0"

    def new_poisson_arrivals(self, start_time, interval, windows = []) -> list:
        poisson_list = []

        last_arrival = 0
        label = ""
        while last_arrival < interval:
            arrival = np.random.exponential(scale=1.0, size=None) / self.k_syn
            last_arrival = last_arrival + arrival

            # we immediately determine the decay time once a transcript comes into existence
            # to do: check if this is the right distribution for a death process
            decay = np.random.exponential(scale=1.0, size=None) / self.k_d
            decay_time = start_time + last_arrival + decay

            arrival_time = start_time + last_arrival

            if last_arrival < interval:
                for window in windows:
                    if window[WINDOW_START] < arrival_time < window[WINDOW_END]:
                        label = window[WINDOW_LABEL]
                    else:
                        label = ""
                poisson_list.append([label, arrival_time, 1, decay_time, -1])

        return poisson_list

    def switch_state(self):
        if self.state == "0":
            self.state = "1"
        else:
            self.state = "0"

    def run_bursts(self, max_minutes, windows):

        dtmc_list = []

        current_time = 0
        poisson_arrivals = []

        # this main part is part Of either Experiment or Transcription
        # l_01 and l_10 are from the Transcription Model
        while current_time < max_minutes:

            state_time = 0
            burst_size = 0

            # this state switching should be moved to Transcription
            if self.state == "0":
                l = self.l_01
                # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
                # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
                # this is a simple way to simulate multiple refractory states?
                alpha = self.nr_refractions
                for i in range(alpha):
                    state_time = state_time + np.random.exponential(scale=1.0, size=None) / l
            else:
                l = self.l_10
                state_time = state_time + np.random.exponential(scale=1.0, size=None) / l

                # current_time is the start of active (ON) state
                # state_time is length of burst
                # now create new Poisson arrivals
                new_arrivals = self.new_poisson_arrivals(current_time, state_time, windows)
                burst_size = len(new_arrivals)
                poisson_arrivals = poisson_arrivals + new_arrivals

            end_time = current_time + state_time
            dtmc_list.append([self.state, current_time, end_time, state_time, burst_size])

            current_time = end_time

            self.switch_state()

        self.df_dtmc = pd.DataFrame(data=dtmc_list,
                                    columns=["state", "begin_time", "end_time", "state_time", "burst_size"])

        self.df_poisson_arrivals = pd.DataFrame(poisson_arrivals,
                                                columns=["label", "arrival", "count_s", "decay", "count_d"])

        return self.df_dtmc, self.df_poisson_arrivals
