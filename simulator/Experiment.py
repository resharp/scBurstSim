from simulator.Transcription import *
import pandas as pd
import math
from typing import NamedTuple

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class ExperimentParams(NamedTuple):
    nr_cells: int
    nr_coordinated_groups: int
    nr_syn_within_strategy: int
    nr_non_syn_within_strategy: int
    windows: list
    freeze: int


class Experiment:

    params = None

    df_all_transcripts = None
    strategies = []
    strategies_file = ""
    tran = None

    trace_id = 0

    transcripts = []
    transcripts_sampled = []
    counts = []
    cell_counts = []

    def __init__(self, params, strategies_file):

        self.params = params
        self.strategies_file = strategies_file

    def read_strategies(self):
        sr = StrategyReader(self.strategies_file)
        self.strategies = sr.select_all()

    # run_alt should return count matrix + percentage of time of gene being ON or OFF during labeling window(s)
    # returns (per label so we can easily generalize to any number of labels):
    # cell_id ; allele_id ; label; perc_label_on; real_count; real_count_unlabeled
    def run_alt(self) -> pd.DataFrame:

        # run_alt is a complete different implementation of run() to be able
        # - to split run and count
        # - insert a sampling step between run and count on cell level
        # - share DTMC traces for different (k_syn, k_d) within a coordination group

        self.read_strategies()

        trans = []  # array of Transcription instances
        alleles = []
        allele_id = 0
        for params in self.strategies:
            tran = Transcription(params)
            trans.append(tran)

            tm_id = params.tm_id
            strategy = params.name
            group_id = params.coord_group
            if not math.isnan(group_id):
                for i_syn in range(self.params.nr_syn_within_strategy):
                    allele_id = allele_id + 1
                    alleles.append((allele_id, strategy, tm_id, group_id, tran))
            else:
                for i_non_syn in range(self.params.nr_non_syn_within_strategy):
                    allele_id = allele_id + 1
                    alleles.append((allele_id, strategy, tm_id, group_id, tran))

        df_alleles = pd.DataFrame(alleles, columns=['allele_id', 'strategy', 'tm_id', 'group_id', 'tran'])
        df_alleles.drop('tran', axis=1, inplace=True)

        for i_c in range(self.params.nr_cells):
            cell_id = i_c + 1

            old_group_id = -1
            new_dtmc_trace = True
            dtmc_list = []
            self.transcripts = []

            for allele_id, strategy, tm_id, group_id, tran in alleles:
                if not math.isnan(group_id):
                    if group_id == old_group_id:
                        new_dtmc_trace = False  # correlation
                        # and dtmc_list stays the same within coordination group
                    else:
                        new_dtmc_trace = True
                        dtmc_list = []
                else:
                    new_dtmc_trace = True
                    dtmc_list = []
                self.tran = tran
                dtmc_list = self.run_transcripts(allele_id, cell_id, tm_id, group_id,
                                                 new_dtmc_trace=new_dtmc_trace, dtmc_list=dtmc_list)
                old_group_id = group_id

            # sample per cell; transcripts -> transcripts_sampled
            # transcripts for all cells are in self.transcripts
            df_cell_transcripts = pd.concat(self.transcripts)

            df_sampled_transcripts = self.sample_transcripts(df_cell_transcripts)

            # TODO: Add both real and sampled transcripts
            df_counts_cell = self.count_transcripts(cell_id, df_cell_transcripts)
            # df_counts_cell = self.count_transcripts(cell_id, df_sampled_transcripts)
            self.cell_counts.append(df_counts_cell)

        df_counts = pd.concat(self.cell_counts)

        df_counts = pd.merge(df_counts, df_alleles, how="left",
                             left_on=['allele_id'],
                             right_on=['allele_id'])

        df_counts.rename(columns={'count_s': 'real_count'}, inplace=True)

        df_counts.fillna("", inplace=True)
        df_counts["fraction"] = df_counts["real_count"] / df_counts["count_all"]
        df_counts["strategy_group"] = df_counts.strategy + "_" + df_counts.group_id.map(str)
        df_counts["allele_label"] = df_counts.allele_id.map(str) + "_" + df_counts.strategy_group

        return df_counts

    # count per cell per allele
    @staticmethod
    def count_transcripts(cell_id, df_cell_transcripts):

        df_label_counts = df_cell_transcripts.groupby(['allele_id', 'trace_id', 'label'])['count_s'].count().reset_index()

        # now we would like to group on allele_id and count and
        df_all_counts = df_cell_transcripts.groupby(['allele_id'])['count_s'].count().reset_index().\
            rename(columns={'count_s': 'count_all'})

        df_counts = pd.merge(df_label_counts, df_all_counts, how='left',
                             left_on=['allele_id'],
                             right_on=['allele_id'])

        df_counts["cell_id"] = cell_id
        return df_counts

    @staticmethod
    def sample_transcripts(df_cell_transcripts):
        efficiency = 0.1  # TODO: parametrize
        len_sample_transcripts = int(efficiency * len(df_cell_transcripts))
        selected_transcripts = df_cell_transcripts.sample(len_sample_transcripts, replace=False)
        return selected_transcripts

    # run returns count matrix + percentage of time of gene being ON or OFF during labeling window(s)
    # returns (per label so we can easily generalize to any number of labels):
    # cell_id ; allele_id ; label; perc_label_on; real_count; real_count_unlabeled
    def run(self) -> pd.DataFrame:

        self.read_strategies()

        # first determine all unique alleles
        alleles = []
        allele_id = 0
        allele_group_id = 0
        for i_cg in range(self.params.nr_coordinated_groups):
            allele_group_id = allele_group_id + 1
            for i_syn in range(self.params.nr_syn_within_strategy):
                allele_id = allele_id + 1
                alleles.append((allele_group_id, allele_id))

        for params in self.strategies:
            self.tran = Transcription(params)

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

        df_counts["fraction"] = df_counts["real_count"] / (df_counts["real_count"] + df_counts["real_count_unlabeled"])
        df_counts["strategy_group"] = df_counts.strategy + "_" + df_counts.allele_group_id.map(str)
        df_counts["allele_label"] = df_counts.allele_id.map(str) + "_" + df_counts.strategy_group

        return df_counts

    # only run transcripts without count
    def run_transcripts(self, allele_id, cell_id, tm_id, group_id, new_dtmc_trace=False, dtmc_list=[]) -> list:

        if new_dtmc_trace:
            self.trace_id = self.trace_id + 1

        df_dtmc, dtmc_list = self.tran.run_bursts(max_minutes=self.params.freeze,
                                                  windows=self.params.windows,
                                                  new_dtmc_trace=new_dtmc_trace,
                                                  dtmc_list=dtmc_list
                                                  )
        # are we still interested in the percentage that the burst was ON?
        # we can derive this from df_dtmc

        # save real transcripts
        df_transcripts = self.tran.df_transcripts
        df_transcripts["cell_id"] = cell_id
        df_transcripts["allele_id"] = allele_id
        df_transcripts["tm_id"] = tm_id
        df_transcripts["group_id"] = group_id
        df_transcripts["trace_id"] = self.trace_id

        self.transcripts.append(df_transcripts)

        return dtmc_list

    # run and count is on allele_id level
    def run_and_count(self, allele_group_id, allele_id, cell_id, new_dtmc_trace=False):

        if new_dtmc_trace:
            self.trace_id = self.trace_id + 1

        # this one runs with complete_trace = True to count in a very inefficient way
        # it needs dfs_labeled_events and df_unlabeled_events
        df_dtmc, dtmc_list = self.tran.run_bursts(max_minutes=self.params.freeze,
                                                  windows=self.params.windows,
                                                  new_dtmc_trace=new_dtmc_trace,
                                                  complete_trace=True)

        # now count on the real transcripts (preprocessed in Transcription tran)
        dfs_labeled_events = self.tran.dfs_labeled_events
        df_unlabeled_events = self.tran.df_unlabeled_events

        # calculate unlabeled counts
        count_un = self.calculate_unlabeled_count(df_unlabeled_events)
        self.counts.append([cell_id, allele_group_id, allele_id, self.trace_id,
                            self.tran.params.name, "", 0, count_un, count_un])

        # calculate labeled counts
        for window in self.params.windows:
            label = window[WINDOW_LABEL]
            count = self.calculate_count(label, dfs_labeled_events)
            perc_burst_time = self.perc_burst_time(df_dtmc, label)
            self.counts.append([cell_id, allele_group_id, allele_id, self.trace_id,
                                self.tran.params.name, label, perc_burst_time, count, count_un])

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

