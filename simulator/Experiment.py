from simulator.Transcription import *
import pandas as pd


# Experiment should include windows and freeze moment (correct name?)
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
    # cell_id ; gene_id ; label; state_80 (80% ON) ; real_count
    def run(self) -> pd.DataFrame:

        poisson_arrivals = []
        counts = []

        # loop cells
        for i_c in range(self.nr_cells):
            for i_a in range(self.nr_alleles):
                df_dtmc, df_poisson_arrivals = self.trans.run_bursts(max_minutes=self.freeze, windows=self.windows)

                # we remember all the individual molecules to sample later and add technical noise
                poisson_arrivals.append(df_poisson_arrivals)

                # calculate counts
                # to do: replace with real code, this is just some scaffolding code
                count = [i_c+1, i_a+1, "", "OFF", 0]
                counts.append(count)

        df_counts = pd.DataFrame(counts, columns=["cell_id", "allele_id", "label", "state", "count"])

        return df_counts
