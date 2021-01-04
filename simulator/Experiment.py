from simulator.StrategyReader import StrategyReader
from simulator.Transcription import *
import pandas as pd
import math
from typing import NamedTuple
import logging

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class ExperimentParams(NamedTuple):
    nr_cells: int
    strategies_file: str
    nr_syn_within_strategy: int
    nr_non_syn_within_strategy: int

    efficiency: float

    windows: list
    fix_time: int


class Experiment:

    params = None

    df_all_transcripts = None
    strategies = []
    tran = None

    trace_id = 0

    transcripts = []
    transcripts_sampled = []
    counts = []
    cell_counts = []

    def __init__(self, params):

        self.params = params

    def read_strategies(self):
        sr = StrategyReader(self.params.strategies_file)
        self.strategies = sr.select_all()

    # run returns count matrix + percentage of time of gene being ON or OFF during labeling window(s)
    # returns (per label so we can easily generalize to any number of labels):
    # allele_id;trace_id;label;real_count;count_all;cell_id;fraction;perc; [extra allele information]
    def run(self) -> pd.DataFrame:

        logging.info("start Experiment.run with efficiency {eff}".format(eff=self.params.efficiency))
        # run is able to
        # - split run and count
        # - insert a sampling step between run and count on cell level
        # - share DTMC traces for different (k_syn, k_d) within a coordination group

        percentages = []
        # to do: extract first part for determining df_alleles
        self.read_strategies()

        logging.info("nr. of strategies read from {file}: {nr}"
                     .format(nr=len(self.strategies), file=self.params.strategies_file))

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

        logging.info("nr. of alleles generated from strategies: {nr_all}".format(nr_all=len(df_alleles)))
        # to do: end of extract

        for i_c in range(self.params.nr_cells):
            cell_id = i_c + 1

            old_group_id = -1
            dtmc_list = []  # list with discrete time markov chain (active <-> inactive) events
            self.transcripts = []  # reinitialize transcripts for all alleles for this cell

            logging.info("start simulation for cell: {cell_id}".format(cell_id=cell_id))

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

                df_dtmc, dtmc_list = self.run_transcripts(allele_id, cell_id, tm_id, group_id,
                                                          new_dtmc_trace=new_dtmc_trace, dtmc_list=dtmc_list)
                # here we have the dtmc_list per allele for one cell
                # from it we can determine the percentage "active period" per window label
                # We join later this information with the labeled Poisson counts
                # this join cannot be done here yet, because the transcripts are going to be sampled on cell level first
                for window in self.params.windows:
                    label = window[WINDOW_LABEL]
                    perc = Experiment.perc_active_state(self.params.windows, df_dtmc, label)

                    # now we have the percentage and we should store it with cell_id, allele_id, label
                    # just add to an array and combine with df_counts at the end of the function
                    percentages.append((cell_id, allele_id, label, perc))

                old_group_id = group_id

            # sample per cell; transcripts -> df_sampled_transcripts
            # transcripts for all alleles for this cell are in self.transcripts
            df_cell_transcripts = pd.concat(self.transcripts)

            df_sampled_transcripts = self.sample_transcripts(df_cell_transcripts)

            # TODO: Add both real and sampled transcripts (now we only add the sampled transcripts)
            # and put this in real_count (which is now a misleading name for sampled counts)

            logging.info("start counting for cell: {cell_id}".format(cell_id=cell_id))
            # df_counts_cell = self.count_transcripts(cell_id, df_cell_transcripts)
            df_counts_cell = self.count_transcripts_per_allele_per_label_for_cell(cell_id, df_sampled_transcripts)

            # now df_counts_cell contains the counts per label per allele for one cell
            # (where there is also a row with an empty label (means unlabeled)

            self.cell_counts.append(df_counts_cell)

        df_counts = pd.concat(self.cell_counts)

        df_counts.rename(columns={'count_s': 'real_count'}, inplace=True)

        df_counts.fillna("", inplace=True)
        # fraction denotes the fraction of cells
        df_counts["fraction"] = df_counts["real_count"] / df_counts["count_all"]
        df_counts.fraction = df_counts.fraction.round(4)

        df_percentages = pd.DataFrame(percentages, columns=['cell_id', 'allele_id', 'label', 'perc'])

        # join df_percentages and join with df_counts
        # motivation outer join: there will be missing rows on both sides
        # df_percentages does not contain the unlabeled (label="") record
        # and df_counts may miss count records for some (cell_id,allele_id,label) combinations
        df_counts = pd.merge(df_counts, df_percentages, how="outer",
                             left_on=['cell_id', 'allele_id', 'label'],
                             right_on=['cell_id', 'allele_id', 'label']).\
            sort_values(['cell_id', 'allele_id', 'label'])

        df_counts = pd.merge(df_counts, df_alleles, how="left",
                             left_on=['allele_id'],
                             right_on=['allele_id'])

        # for displaying purposes
        df_counts["strategy_group"] = df_counts.strategy + "_" + df_counts.group_id.map(str)
        df_counts["allele_label"] = df_counts.allele_id.map(str) + "_" + df_counts.strategy_group

        return df_counts

    @staticmethod
    def count_transcripts_per_allele_per_label_for_cell(cell_id, df_cell_transcripts):

        df_label_counts = df_cell_transcripts.groupby(['allele_id', 'trace_id', 'label'])['count_s'].count().reset_index()

        # now we would like to group on allele_id and count
        df_all_counts = df_cell_transcripts.groupby(['allele_id'])['count_s'].count().reset_index().\
            rename(columns={'count_s': 'count_all'})

        # add allele information
        df_counts = pd.merge(df_label_counts, df_all_counts, how='left',
                             left_on=['allele_id'],
                             right_on=['allele_id'])

        df_counts["cell_id"] = cell_id
        return df_counts

    def sample_transcripts(self, df_cell_transcripts):
        len_sample_transcripts = int(self.params.efficiency * len(df_cell_transcripts))
        selected_transcripts = df_cell_transcripts.sample(len_sample_transcripts, replace=False)
        return selected_transcripts

    # only run transcripts without count
    def run_transcripts(self, allele_id, cell_id, tm_id, group_id, new_dtmc_trace=False, dtmc_list=[]) -> list:

        if new_dtmc_trace:
            self.trace_id = self.trace_id + 1

        df_dtmc, dtmc_list = self.tran.run_bursts(max_minutes=self.params.fix_time,
                                                  windows=self.params.windows,
                                                  new_dtmc_trace=new_dtmc_trace,
                                                  dtmc_list=dtmc_list
                                                  )

        # save real transcripts (to sample later on cell level)
        df_transcripts = self.tran.df_transcripts
        df_transcripts["cell_id"] = cell_id
        df_transcripts["allele_id"] = allele_id
        df_transcripts["tm_id"] = tm_id
        df_transcripts["group_id"] = group_id
        df_transcripts["trace_id"] = self.trace_id

        # add to self.transcripts
        self.transcripts.append(df_transcripts)

        return df_dtmc, dtmc_list

    @classmethod    # this is a class method because it is also used out of context of an instance of Experiment
    def perc_active_state(cls, windows, df_dtmc, label) -> int:

        for window in windows:
            if window[WINDOW_LABEL] == label:
                start = window[WINDOW_START]
                end = window[WINDOW_END]

        interval = end - start

        active_states = df_dtmc[(df_dtmc["state"] == "1") &
                                (df_dtmc["end_time"] >= start) &
                                (df_dtmc["begin_time"] <= end)].copy(deep=True)

        if len(active_states) > 0:
            active_states["begin_time"] = np.maximum(active_states["begin_time"], start)
            active_states["end_time"] = np.minimum(active_states["end_time"], end)
            active_states["state_time"] = active_states["end_time"] - active_states["begin_time"]

            sum_state_time = active_states["state_time"].sum()
            perc = sum_state_time / interval
        else:
            perc = 0

        return perc

