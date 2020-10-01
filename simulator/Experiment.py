from simulator.Transcription import *
import pandas as pd
from typing import NamedTuple

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class ExperimentParams(NamedTuple):
    nr_cells: int
    nr_coordinated_groups: int
    nr_trace_copies: int
    windows: list
    freeze: int


class Experiment:

    params = None

    df_all_transcripts = None
    strategies = []
    strategies_file = ""
    trans = None

    trace_id = 0

    transcripts = []
    counts = []

    def __init__(self, params, strategies_file):

        self.params = params
        self.strategies_file = strategies_file

    def read_strategies(self):
        sr = StrategyReader(self.strategies_file)
        self.strategies = sr.select_all()

    # run returns count matrix + true positives for gene being ON or OFF during labeling window(s)
    # returns (per label so we can easily generalize to two types of labels):
    # cell_id ; allele_id ; label; perc_label_on; real_count; real_count_unlabeled
    def run(self) -> pd.DataFrame:

        self.read_strategies()

        # first determine all unique alleles
        alleles = []
        allele_id = 0
        allele_group_id = 0
        for i_cg in range(self.params.nr_coordinated_groups):
            allele_group_id = allele_group_id + 1
            for i_tc in range(self.params.nr_trace_copies):
                allele_id = allele_id + 1
                alleles.append((allele_group_id, allele_id))

        for params in self.strategies:
            self.trans = Transcription(params)

            for i_c in range(self.params.nr_cells):
                cell_id = i_c + 1
                old_group_id = -1
                for allele_group_id, allele_id in alleles:

                    # within every cell the coordinated allele groups should have the same dtmc trace
                    # only change trace when there is a new group
                    if allele_group_id == old_group_id:
                        new_dtmc_trace = False
                    else:
                        new_dtmc_trace = True

                    self.run_and_count(allele_group_id, allele_id, cell_id, new_dtmc_trace=new_dtmc_trace)
                    old_group_id = allele_group_id

        df_counts = pd.DataFrame(
            self.counts, columns=["cell_id", "allele_group_id", "allele_id", "trace_id",
                                  "strategy", "label", "perc_label_on", "real_count", "real_count_unlabeled"])

        self.df_all_transcripts = pd.concat(self.transcripts)

        df_counts["fraction"] = df_counts["real_count"] / (df_counts["real_count"] + df_counts["real_count_unlabeled"])
        df_counts["strategy_group"] = df_counts.strategy + "_" + df_counts.allele_group_id.map(str)

        return df_counts

    # run and count is on allele_id level
    def run_and_count(self, allele_group_id, allele_id, cell_id, new_dtmc_trace=False):

        if new_dtmc_trace:
            self.trace_id = self.trace_id + 1

        df_dtmc, df_events = self.trans.run_bursts(max_minutes=self.params.freeze,
                                                   windows=self.params.windows,
                                                   new_dtmc_trace=new_dtmc_trace)
        df_events["cell_id"] = cell_id
        df_events["allele_id"] = allele_id

        # save real transcripts
        df_transcripts = self.trans.df_transcripts
        df_transcripts["cell_id"] = cell_id
        df_transcripts["allele_id"] = allele_id
        self.transcripts.append(df_transcripts)

        # now count on the real transcripts (preprocessed in Transcription)
        # TODO: refactor and count on the sampled transcripts
        dfs_labeled_events = self.trans.dfs_labeled_events
        df_unlabeled_events = self.trans.df_unlabeled_events

        # calculate unlabeled counts
        count_un = self.calculate_unlabeled_count(df_unlabeled_events)
        self.counts.append([cell_id, allele_group_id, allele_id, self.trace_id,
                            self.trans.params.name, "", 0, count_un, count_un])

        # calculate labeled counts
        for window in self.params.windows:
            label = window[WINDOW_LABEL]
            count = self.calculate_count(label, dfs_labeled_events)
            perc_burst_time = self.perc_burst_time(df_dtmc, label)
            self.counts.append([cell_id, allele_group_id, allele_id, self.trace_id,
                                self.trans.params.name, label, perc_burst_time, count, count_un])

    def calculate_count(self, par_label, dfs_labeled_events):

        cum_count = 0
        detected = False

        for label, df_labeled_events in dfs_labeled_events:
            if label == par_label:
                detected = True
                if len(df_labeled_events) > 0:
                    df_before_freeze = df_labeled_events[df_labeled_events.arrival < self.params.freeze]
                    if len(df_before_freeze) > 0:
                        cum_count = df_before_freeze.iloc[-1]["cum_count"]
                    else:
                        cum_count = 0
                else:
                    cum_count = 0
        if not detected:
            cum_count = 0

        return cum_count

    def calculate_unlabeled_count(self, df_unlabeled_events):

        cum_count_unlabeled = 0
        if len(df_unlabeled_events) > 0:
            df_before_freeze = df_unlabeled_events[df_unlabeled_events.arrival < self.params.freeze]
            if len(df_before_freeze) > 0:
                if len(df_before_freeze) > 0:
                    cum_count_unlabeled = df_before_freeze.iloc[-1]["cum_count"]
                else:
                    cum_count_unlabeled = 0
        else:
            cum_count_unlabeled = 0
        return cum_count_unlabeled

    def perc_burst_time(self, df_dtmc, label) -> int:

        for window in self.params.windows:
            if window[WINDOW_LABEL] == label:
                start = window[WINDOW_START]
                end = window[WINDOW_END]

        interval = end - start

        active_states = df_dtmc[(df_dtmc["state"] == "1") &
                                (df_dtmc["end_time"] >= start) &
                                (df_dtmc["begin_time"] <= end)].copy(deep=True)

        active_states["begin_time"] = np.maximum(active_states["begin_time"], start)
        active_states["end_time"] = np.minimum(active_states["end_time"], end)
        active_states["state_time"] = active_states["end_time"] - active_states["begin_time"]

        if len(active_states) > 0:
            sum_state_time = active_states["state_time"].sum()
            perc = sum_state_time / interval
        else:
            perc = 0

        return perc

