import numpy as np
import pandas as pd
from typing import NamedTuple

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class TranscriptParams(NamedTuple):
    l_01: int
    l_10: int
    k_syn: int              # k_syn = synthesis rate = transcription rate
    nr_refractions: int
    k_d: int                # decay is strictly not part of transcription but we include it in the model


class Transcription:

    params = None

    state = "0"

    df_dtmc = None
    df_poisson_arrivals = None
    df_events = None

    df_labeled_events = []
    df_unlabeled_events = None

    def __init__(self, params):

        self.params = params

        self.state = "0"

    def new_poisson_arrivals(self, start_time, interval, windows=[]) -> list:
        poisson_list = []

        last_arrival = 0
        label = ""
        while last_arrival < interval:
            arrival = np.random.exponential(scale=1.0, size=None) / self.params.k_syn
            last_arrival = last_arrival + arrival

            # we immediately determine the decay time once a transcript comes into existence
            # TODO: check if this is the right distribution for a death process
            decay = np.random.exponential(scale=1.0, size=None) / self.params.k_d
            decay_time = start_time + last_arrival + decay

            arrival_time = start_time + last_arrival

            if last_arrival < interval:
                for window in windows:
                    if window[WINDOW_START] < arrival_time < window[WINDOW_END]:
                        label = window[WINDOW_LABEL]
                        break
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

            if self.state == "0":
                l = self.params.l_01
                # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
                # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
                # this is a simple way to simulate multiple refractory states?
                alpha = self.params.nr_refractions
                for i in range(alpha):
                    state_time = state_time + np.random.exponential(scale=1.0, size=None) / l
            else:
                l = self.params.l_10
                alpha = self.params.nr_refractions
                for i in range(alpha):
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

        # we will now put the arrivals and decays in one table self.df_events and sort by time ..
        self.df_events = self.sort_events()

        # .. enable cumulative sums
        self.sum_labeled_events(windows)

        self.sum_unlabeled_events()

        return self.df_dtmc, self.df_events

    # here we put (time of) arrivals and decays in the same column to sort and to enable cumulative sums
    # however we lose the single molecule information this way
    def sort_events(self):

        df_decays = self.df_poisson_arrivals[['label','decay', "count_d"]].\
            rename(columns={'decay': 'arrival', 'count_d': 'count_s'})

        df_poisson_arrivals = self.df_poisson_arrivals[["label", "arrival", "count_s"]]
        df_events = df_poisson_arrivals.append(df_decays).sort_values(by="arrival")

        return df_events

    def sum_labeled_events(self, windows):

        self.df_labeled_events = []
        for window in windows:
            df_labeled = (self.df_events[self.df_events.label == window[WINDOW_LABEL]]).\
                copy(deep=True)
            df_labeled['cum_count'] = df_labeled['count_s'].cumsum()
            self.df_labeled_events.append([window[WINDOW_LABEL], df_labeled])

    def sum_unlabeled_events(self):
        self.df_unlabeled_events = (self.df_events[self.df_events.label == ""]).copy(deep=True)
        self.df_unlabeled_events['cum_count'] = self.df_unlabeled_events['count_s'].cumsum()
