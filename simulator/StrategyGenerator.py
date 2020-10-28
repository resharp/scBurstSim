import logging

import numpy as np
import pandas as pd

from simulator.Transcription import TranscriptParams


class StrategyGenerator:

    filename = ""

    range_k_on = []
    range_k_off = []
    range_k_sys = []
    range_k_d = []

    counter = 0

    def __init__(self, range_k_on, range_k_off, range_k_syn, range_k_d, filename):

        self.range_k_on = range_k_on
        self.range_k_off = range_k_off
        self.range_k_sys = range_k_syn
        self.range_k_d = range_k_d

        self.filename = filename

        self.log_parameters()

    def log_parameters(self):
        logging.info("StrategyGenerator started with parameters:")

        logging.info("k_on from {min} to {max} / minute".format(min=self.range_k_on[0], max=self.range_k_on[1]))
        logging.info("k_off from {min} to {max} / minute".format(min=self.range_k_off[0], max=self.range_k_off[1]))

        min_k_on = self.range_k_on[0]
        max_k_on = self.range_k_on[1]
        min_k_off = self.range_k_off[0]
        max_k_off = self.range_k_off[1]

        max_chance_on = round(100 * max_k_on / (max_k_on + min_k_off ), 3)
        min_chance_on = round(100 * min_k_on / ( min_k_on + max_k_off), 3)
        logging.info("percentage ON from {min} to {max} percent".format(min=min_chance_on, max=max_chance_on))

        logging.info("k_syn from {min} to {max} / minute".format(min=self.range_k_sys[0], max=self.range_k_sys[1]))
        k_syn_limits_hours = [limit * 60 for limit in self.range_k_sys]
        logging.info("k_syn from {min} to {max} / hour".format(min=k_syn_limits_hours[0], max=k_syn_limits_hours[1]))

        logging.info("k_d from {min} to {max} / minute".format(min=self.range_k_d[0], max=self.range_k_d[1]))

        half_lives = [round(np.log(2)/limit, 2) for limit in self.range_k_d]

        logging.info("half-lives from {min} to {max} minutes".format(min=half_lives[0], max=half_lives[1]))

    def get_random_parameters(self):

        name, k_on, k_off, group_id, k_syn, k_d = self.get_random()

        tp = TranscriptParams(k_on=k_on, k_off=k_off, k_syn=k_syn, k_d=k_d,
                              nr_refractions=2, name="generated", coord_group=0, tm_id=0)

        return tp

    def get_random(self):

        # we want the ON times shorter than the OFF times: average length ON ~ 1/k_off < 1/k_on
        # => k_off > k_on
        # TODO QUESTION: Is it possible having ON times larger than OFF times?
        k_on = self.sample_value(self.range_k_on)
        min_k_off = max(k_on, self.range_k_off[0])
        k_off = self.sample_value([min_k_off, self.range_k_off[1]])

        # it makes sense that k_syn > k_off always (because there would be no burst)
        min_k_syn = max(k_off, self.range_k_sys[0])
        k_syn = self.sample_value([min_k_syn, self.range_k_sys[1]], exponential=True)

        k_d = self.sample_value(self.range_k_d)

        self.counter = self.counter + 1
        name = "generated_" + str(self.counter)

        return [name, k_on, k_off, np.nan, k_syn, k_d]

    def generate_and_write_strategies(self, nr_of_alleles):

        columns = ['name', 'k_on', 'k_off', 'coord_group', 'k_syn', 'k_d']

        data = []

        for i in range(0, nr_of_alleles):
            data.append(self.get_random())

        df_strategies = pd.DataFrame(data=data, columns=columns)
        df_strategies.to_csv(path_or_buf=self.filename, sep=';', index=False)

    @staticmethod
    def sample_value(limits, exponential=False):
        decimals = 4

        if exponential:
            max_value = limits[1]
            lambda_exp = max_value

            ret_value = max_value + 1
            while ret_value > max_value:
                ret_value = round(np.random.exponential(scale=1 / lambda_exp), decimals)
        else:
            ret_value = round(np.random.uniform(limits[0], limits[1]), decimals)

        return ret_value
