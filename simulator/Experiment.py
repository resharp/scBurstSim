from simulator.Transcription import *
import pandas as pd


class Experiment:
    nr_cells = 0
    nr_alleles = 0
    windows = []
    freeze = 0
    trans = None

    def __init__(self, nr_cells, nr_alleles, params, windows, freeze):

        self.nr_cells = nr_cells
        self.nr_alleles = nr_alleles
        self.windows = windows
        self.freeze = freeze

        self.trans = Transcription(params)

    # run returns count matrix + TPs for gene being ON or OFF during labeling window(s)
    # so something like (per label so we can easily generalize to two types of labels):
    # cell_id ; gene_id ; label; percentage_on ; real_count
    def run(self) -> pd.DataFrame:

        poisson_arrivals = []
        counts = []

        # loop cells
        for i_c in range(self.nr_cells):
            for i_a in range(self.nr_alleles):
                cell_id = i_c + 1
                allele_id = i_a + 1

                df_dtmc, df_poisson_arrivals = self.trans.run_bursts(max_minutes=self.freeze, windows=self.windows)

                df_labeled_arrivals = self.trans.df_labeled_arrivals
                df_unlabeled_arrivals = self.trans.df_unlabeled_arrivals

                # we remember all the individual molecules to sample later and add technical noise
                poisson_arrivals.append(df_poisson_arrivals)

                for label, count in self.calculate_count(df_labeled_arrivals, df_unlabeled_arrivals):
                    counts.append([cell_id, allele_id, label, "OFF", count])

        df_counts = pd.DataFrame(counts, columns=["cell_id", "allele_id", "label", "state", "real_count"])

        return df_counts

    def calculate_count(self, df_labeled_arrivals, df_unlabeled_arrivals):

        counts = []

        for df_labeled_arrival in df_labeled_arrivals:

            if len(df_labeled_arrival) > 0:
                label = df_labeled_arrival.iloc[0]["label"]

                df_before_freeze = df_labeled_arrival[df_labeled_arrival.arrival < self.freeze]
                if len(df_before_freeze) > 0:
                    cum_count = df_before_freeze.iloc[-1]["cum_count"]
                else:
                    cum_count = 0
                counts.append([label, cum_count])

        return counts

