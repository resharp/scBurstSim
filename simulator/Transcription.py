import numpy as np
import pandas as pd
from typing import NamedTuple

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class TranscriptParams(NamedTuple):
    # transcription matrix
    k_01: float             # part of transition matrix
    k_10: float             # part of transition matrix
    nr_refractions: int

    # high/low expression
    k_syn: float            # k_syn = synthesis rate = transcription rate
    k_d: float              # decay is strictly not part of transcription but we include it in the model

    # coordination (prerequisite: same transition matrix)
    coord_group: int

    name: str


class Transcription:

    params = None

    state = "0"

    df_dtmc = None
    df_transcripts = None
    df_events = None

    dfs_labeled_events = []
    df_unlabeled_events = None

    dtmc_list = []  # dtmc list should be stored in Transcription

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

    def run_bursts(self, max_minutes, windows, new_dtmc_trace=True) -> (pd.DataFrame, pd.DataFrame):

        if new_dtmc_trace:
            self.dtmc_list = self.create_dtmc_list(max_minutes)

        self.df_dtmc = pd.DataFrame(data=self.dtmc_list,
                                    columns=["state", "begin_time", "end_time", "state_time"])

        self.df_transcripts = self.create_transcripts(windows)

        # we will now put the arrivals and decays in one table self.df_events and sort by time ..
        self.df_events = self.sort_events()

        # .. enable cumulative sums
        self.sum_labeled_events(windows)

        self.sum_unlabeled_events()

        return self.df_dtmc, self.df_events

    def create_dtmc_list(self, max_minutes):
        dtmc_list = []
        current_time = 0
        while current_time < max_minutes:

            state_time = 0

            if self.state == "0":
                k = self.params.k_01
                # we could get a peaked distribution of waiting times by repeating (setting alpha > 1)
                # see alpha in https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
                # this is a simple way to simulate multiple refractory states?
                alpha = self.params.nr_refractions
                k = k * alpha
                for i in range(alpha):
                    state_time = state_time + np.random.exponential(scale=1.0, size=None) / k
            else:
                k = self.params.k_10
                state_time = np.random.exponential(scale=1.0, size=None) / k

            end_time = current_time + state_time

            dtmc = [self.state, current_time, end_time, state_time]
            dtmc_list.append(dtmc)

            current_time = end_time

            self.switch_state()
        return dtmc_list

    def create_transcripts(self, windows):
        poisson_arrivals = []
        for dtmc in self.dtmc_list:
            state = dtmc[0]
            if state == "1":
                current_time = dtmc[1]
                state_time = dtmc[3]
                new_arrivals = self.new_poisson_arrivals(current_time, state_time, windows)
                poisson_arrivals = poisson_arrivals + new_arrivals
        df_transcripts = pd.DataFrame(poisson_arrivals,
                                      columns=["label", "arrival", "count_s", "decay", "count_d"])
        return df_transcripts

    # here we put (time of) arrivals and decays in the same column to sort and to enable cumulative sums
    # however we lose the single molecule information this way
    def sort_events(self):

        df_decays = self.df_transcripts[['label', 'decay', "count_d"]].\
            rename(columns={'decay': 'arrival', 'count_d': 'count_s'})

        df_poisson_arrivals = self.df_transcripts[["label", "arrival", "count_s"]]
        df_events = df_poisson_arrivals.append(df_decays).sort_values(by="arrival")

        return df_events

    def sum_labeled_events(self, windows):

        self.dfs_labeled_events = []
        for window in windows:
            df_labeled = (self.df_events[self.df_events.label == window[WINDOW_LABEL]]).\
                copy(deep=True)
            df_labeled['cum_count'] = df_labeled['count_s'].cumsum()
            self.dfs_labeled_events.append([window[WINDOW_LABEL], df_labeled])

    def sum_unlabeled_events(self):
        self.df_unlabeled_events = (self.df_events[self.df_events.label == ""]).copy(deep=True)
        self.df_unlabeled_events['cum_count'] = self.df_unlabeled_events['count_s'].cumsum()


class StrategyReader:
    filename = ""
    df_strategies = None

    def __init__(self, filename):
        self.filename = filename

    def select_all(self):

        self.read_strategies()

        params_list = [TranscriptParams(k_01=item.k_01, k_10=item.k_10, nr_refractions=2,
                                        k_syn=item.k_syn, k_d=item.k_d,
                                        coord_group=item.coord_group,
                                        name=item['name'])
                       for id_dummy, item in self.df_strategies.iterrows()]
        return params_list

    def get(self, strategy):

        self.read_strategies()

        df_strategy = self.df_strategies[self.df_strategies.name == strategy]
        if len(df_strategy) == 0:
            raise RuntimeError("{strategy} is no valid strategy name! "
                               "See strategy names in file {filename}".
                               format(strategy=strategy, filename=self.filename))
        params = self.convert_to_params(df_strategy)
        return params

    def get_random(self):
        self.read_strategies()

        df_strategy = self.df_strategies.sample(1)

        return self.convert_to_params(df_strategy)

    def read_strategies(self):
        if self.df_strategies is None:
            self.df_strategies = pd.read_csv(self.filename, sep=";", comment="#")

    @staticmethod
    def convert_to_params(df_strategy):
        params = TranscriptParams(k_01=df_strategy.k_01.item(), k_10=df_strategy.k_10.item(), nr_refractions=2,
                                  k_syn=df_strategy.k_syn.item(), k_d=df_strategy.k_d.item(),
                                  coord_group=df_strategy.coord_group.item(),
                                  name=df_strategy.name.item())
        return params
