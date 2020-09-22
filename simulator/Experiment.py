from simulator.Transcription import *
import pandas as pd
from typing import NamedTuple

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class ExperimentParams(NamedTuple):
    nr_cells: int
    nr_alleles: int
    windows: list
    freeze: int
    trans_params: object


class Experiment:

    params = None

    trans = None
    df_all_transcripts = None

    def __init__(self, params):

        self.params = params
        self.trans = Transcription(params.trans_params)

    # run returns count matrix + true positives for gene being ON or OFF during labeling window(s)
    # returns:
    # (per label so we can easily generalize to two types of labels):
    # cell_id ; allele_id ; label; perc_label_on; real_count; real_count_unlabeled
    def run(self) -> pd.DataFrame:

        transcripts = []
        counts = []

        # loop cells
        for i_c in range(self.params.nr_cells):
            for i_a in range(self.params.nr_alleles):
                cell_id = i_c + 1
                allele_id = i_a + 1

                df_dtmc, df_events = self.trans.run_bursts(max_minutes=self.params.freeze
                                                                     , windows=self.params.windows)

                df_labeled_events = self.trans.df_labeled_events
                df_unlabeled_events = self.trans.df_unlabeled_events

                df_events["cell_id"] = cell_id
                df_events["allele_id"] = allele_id

                # TODO: we should remember all the individual poisson arrivals and decays on single transcript level
                df_transcripts = self.trans.df_transcripts
                df_transcripts["cell_id"] = cell_id
                df_transcripts["allele_id"] = allele_id

                transcripts.append(df_transcripts)

                # calculate unlabeled counts
                count_un = self.calculate_unlabeled_count(df_unlabeled_events)
                counts.append([cell_id, allele_id, "", 0, count_un, count_un])

                # calculate labeled counts
                for window in self.params.windows:
                    label = window[WINDOW_LABEL]
                    count = self.calculate_count(label, df_labeled_events)
                    perc_burst_time = self.perc_burst_time(df_dtmc, label)
                    counts.append([cell_id, allele_id, label, perc_burst_time, count, count_un])

        df_counts = pd.DataFrame(
            counts, columns=["cell_id", "allele_id", "label", "perc_label_on", "real_count", "real_count_unlabeled"])

        self.df_all_transcripts = pd.concat(transcripts)

        return df_counts

    def calculate_count(self, par_label, df_labeled_events):

        cum_count = 0
        detected = False

        for label, df_labeled_event in df_labeled_events:
            if label == par_label:
                detected = True
                if len(df_labeled_event) > 0:
                    df_before_freeze = df_labeled_event[df_labeled_event.arrival < self.params.freeze]
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

