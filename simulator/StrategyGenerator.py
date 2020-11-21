import logging

import numpy as np
import pandas as pd

from simulator.Transcription import TranscriptParams


class StrategyGenerator:

    filename = ""

    range_k_on_exp = []
    range_k_off_exp = []
    range_k_syn_exp = []
    range_k_d_exp = []

    counter = 0

    def __init__(self, range_k_on_exp, range_k_off_exp, range_k_syn_exp, range_k_d_exp, filename):

        self.range_k_on_exp = range_k_on_exp
        self.range_k_off_exp = range_k_off_exp
        self.range_k_syn_exp = range_k_syn_exp
        self.range_k_d_exp = range_k_d_exp

        self.filename = filename

        self.log_parameters()

    def log_parameters(self):

        range_k_on = [10**k/60 for k in self.range_k_on_exp]
        range_k_off = [10**k/60 for k in self.range_k_off_exp]
        range_k_syn = [10**k/60 for k in self.range_k_syn_exp]
        range_k_d = [10**k/60 for k in self.range_k_d_exp]

        logging.info("StrategyGenerator started with parameters:")

        logging.info("k_on from {min} to {max} / minute".format(min=range_k_on[0], max=range_k_on[1]))
        logging.info("k_off from {min} to {max} / minute".format(min=range_k_off[0], max=range_k_off[1]))

        min_k_on = range_k_on[0]
        max_k_on = range_k_on[1]
        min_k_off = range_k_off[0]
        max_k_off = range_k_off[1]

        max_chance_on = round(100 * max_k_on / (max_k_on + min_k_off ), 3)
        min_chance_on = round(100 * min_k_on / ( min_k_on + max_k_off), 3)
        logging.info("percentage ON from {min} to {max} percent".format(min=min_chance_on, max=max_chance_on))

        logging.info("k_syn from {min} to {max} / minute".format(min=range_k_syn[0], max=range_k_syn[1]))
        k_syn_limits_hours = [limit * 60 for limit in range_k_syn]
        logging.info("k_syn from {min} to {max} / hour".format(min=k_syn_limits_hours[0], max=k_syn_limits_hours[1]))

        logging.info("k_d from {min} to {max} / minute".format(min=range_k_d[0], max=range_k_d[1]))

        half_lives = [round(np.log(2)/limit, 2) for limit in range_k_d]

        logging.info("half-lives from {min} to {max} minutes".format(min=half_lives[0], max=half_lives[1]))

    # get_random_parameters is for getting random parameters directly
    # (not by first writing to a csv file
    def get_random_parameters(self):

        name, k_on, k_off, group_id, k_syn, k_d = self.create_random()

        tp = TranscriptParams(k_on=k_on, k_off=k_off, k_syn=k_syn, k_d=k_d,
                              nr_refractions=1, name="generated", coord_group=0, tm_id=0)

        return tp

    def round_sig(self, x, n=4):

        round_to_n = round(x, -int(np.floor(np.log10(x))) + (n - 1))
        # round_to_n = round(x, n)

        return round_to_n

    def create_random(self, k_on_first=False):

        if k_on_first:
            k_on_exp = self.sample_value(self.range_k_on_exp)

            k_on = self.round_sig(10 ** k_on_exp / 60)

            # we want the ON times shorter than the OFF times: average length ON ~ 1/k_off < 1/k_on
            # => k_off > k_on
            min_k_off_exp = max(k_on_exp, self.range_k_off_exp[0])

            k_off_exp = self.sample_value([min_k_off_exp, self.range_k_off_exp[1]])

            k_off = self.round_sig(10 ** k_off_exp / 60)
        else:
            k_off_exp = self.sample_value(self.range_k_off_exp)
            k_off = self.round_sig(10 ** k_off_exp / 60)

            # we want the ON times shorter than the OFF times: average length ON ~ 1/k_off < 1/k_on
            # => k_on < k_off
            max_k_on_exp = min(k_off_exp, self.range_k_on_exp[1])
            k_on_exp = self.sample_value([self.range_k_off_exp[0], max_k_on_exp])

            k_on = self.round_sig(10 ** k_on_exp / 60)

        # it makes sense that k_syn > k_off always (because there would be no burst)
        min_k_syn_exp = max(k_off_exp, self.range_k_syn_exp[0])
        k_syn_exp = self.sample_value([min_k_syn_exp, self.range_k_syn_exp[1]])

        k_syn = self.round_sig(10 ** k_syn_exp / 60)

        k_d_exp = self.sample_value(self.range_k_d_exp)
        k_d = self.round_sig(10 ** k_d_exp / 60)

        self.counter = self.counter + 1
        name = "generated_" + str(self.counter)

        return [name, k_on, k_off, np.nan, k_syn, k_d]

    def generate_and_write_strategies(self, nr_of_alleles, k_on_first=True):

        columns = ['name', 'k_on', 'k_off', 'coord_group', 'k_syn', 'k_d']

        data = []

        for i in range(0, nr_of_alleles):
            data.append(self.create_random(k_on_first))

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
