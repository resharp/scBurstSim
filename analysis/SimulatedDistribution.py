from utils.utils import round_sig


class SimulatedDistribution:

    df_counts = None
    nr_cells = 0

    def __init__(self, df_counts, nr_cells):

        self.df_counts = df_counts

        self.nr_cells = nr_cells

    def create(self, strategy, label=None):

        if label is None:
            df_counts = self.df_counts
            measure = "count_all"
        else:
            df_counts = self.df_counts[self.df_counts.label == label]
            measure = "real_count"

        df_allele_cell_counts = df_counts.groupby(['allele_id', 'strategy', 'cell_id'])[measure].max().reset_index()

        # TODO: improve; filtering on strategy can be done earlier
        df_one_allele_counts = df_allele_cell_counts[df_allele_cell_counts.strategy == strategy]

        df_one_allele_counts = df_one_allele_counts.set_index('cell_id'). \
            reindex(range(1, self.nr_cells + 1)).fillna(0).reset_index()

        df_distribution = df_one_allele_counts.groupby(measure)['cell_id'].count().to_frame().reset_index()
        df_distribution[measure] = df_distribution[measure].astype(int)

        max_count = df_distribution[measure].max()

        df_distribution = df_distribution.set_index(measure).reindex(range(0, max_count + 1)).fillna(
            0).reset_index()

        # cell_id contains the number of cells with the "measure" value
        # (measure may be count_all or real_count (for label))
        nr_cells = int(sum(df_distribution.cell_id))
        df_distribution["chance"] = df_distribution.cell_id / nr_cells

        # total_chance = df_distribution["chance"].sum() # for debugging, should add up to 1
        df_distribution['weighted'] = df_distribution[measure] * df_distribution.cell_id
        sum_weighted = df_distribution['weighted'].sum()
        real_mean = round_sig(sum_weighted / nr_cells, 4)

        return df_distribution, real_mean

