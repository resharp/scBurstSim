import pandas as pd

from simulator.Transcription import TranscriptParams


class StrategyReader:
    filename = ""
    df_strategies = None

    def __init__(self, filename):
        self.filename = filename

    def select_all(self):

        self.read_strategies()

        params_list = [TranscriptParams(k_on=item.k_on, k_off=item.k_off, nr_refractions=2,
                                        tm_id=item.tm_id,
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

            # we sort on transcription matrix (k_on, k_off), coordination group for synchronization
            # and finally the synthesis and decay rate to be able to share DTMC traces
            # between alleles with different (k_syn, k_d) in the same coordination group
            self.df_strategies.sort_values(by=['k_on', 'k_off', 'coord_group', 'k_syn', 'k_d'], inplace=True)

            # add a unique transition matrix id for every unique combination of (k_on, k_off)
            df_tms = self.df_strategies.groupby(['k_on', 'k_off']).max().reset_index()[['k_on', 'k_off']]
            df_tms['count'] = 1
            df_tms['tm_id'] = df_tms['count'].cumsum()
            df_tms.drop('count', axis=1, inplace=True)
            self.df_strategies = pd.merge(self.df_strategies, df_tms, how='left',
                                          left_on=['k_on', 'k_off'],
                                          right_on=['k_on', 'k_off'])

            self.df_strategies["fraction_ON"] = self.df_strategies.k_on / \
                                                ( self.df_strategies.k_on + self.df_strategies.k_off )
            self.df_strategies["fraction_OFF"] = self.df_strategies.k_off / \
                                                ( self.df_strategies.k_on + self.df_strategies.k_off )

    @staticmethod
    def convert_to_params(df_strategy):
        params = TranscriptParams(k_on=df_strategy.k_on.item(), k_off=df_strategy.k_off.item(), nr_refractions=2,
                                  tm_id=df_strategy.tm_id.item(),
                                  k_syn=df_strategy.k_syn.item(), k_d=df_strategy.k_d.item(),
                                  coord_group=df_strategy.coord_group.item(),
                                  name=df_strategy.name.item())
        return params
