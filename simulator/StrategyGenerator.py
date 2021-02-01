import logging

import numpy as np
import pandas as pd
import random

from simulator.Transcription import TranscriptParams
from utils.utils import round_sig

from itertools import product


class StrategyGenerator:

    def __init__(self, range_k_on_exp, range_k_off_exp, range_k_syn_exp, range_k_d_exp, filename):

        self.counter = 0
        self.range_k_on_exp = range_k_on_exp
        self.range_k_off_exp = range_k_off_exp
        self.range_k_syn_exp = range_k_syn_exp
        self.range_k_d_exp = range_k_d_exp

        self.filename = filename

        self.log_parameters()

    def log_parameters(self):

        # convert to ranges expressed in rates per minute
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

        k_on, k_off, group_id, k_syn, k_d = self.create_random()

        tp = TranscriptParams(k_on=k_on, k_off=k_off, k_syn=k_syn, k_d=k_d,
                              nr_refractions=1, name="generated", coord_group=0, tm_id=0, tran_type="S")

        return tp

    def create_random(self, k_on_first=False):

        if k_on_first:
            k_on_exp = self.sample_value(self.range_k_on_exp)

            # convert to units per minute
            k_on = round_sig(10 ** k_on_exp / 60)

            # we want the ON times shorter than the OFF times: average length ON ~ 1/k_off < 1/k_on
            # => k_off > k_on
            min_k_off_exp = max(k_on_exp, self.range_k_off_exp[0])

            k_off_exp = self.sample_value([min_k_off_exp, self.range_k_off_exp[1]])

            k_off = round_sig(10 ** k_off_exp / 60)
        else:
            k_off_exp = self.sample_value(self.range_k_off_exp)
            k_off = round_sig(10 ** k_off_exp / 60)

            # we want the ON times shorter than the OFF times: average length ON ~ 1/k_off < 1/k_on
            # => k_on < k_off
            max_k_on_exp = min(k_off_exp, self.range_k_on_exp[1])
            k_on_exp = self.sample_value([self.range_k_off_exp[0], max_k_on_exp])

            k_on = round_sig(10 ** k_on_exp / 60)

        # it makes sense that k_syn > k_off always (because there would be no burst)
        min_k_syn_exp = max(k_off_exp, self.range_k_syn_exp[0])
        k_syn_exp = self.sample_value([min_k_syn_exp, self.range_k_syn_exp[1]])

        k_syn = round_sig(10 ** k_syn_exp / 60)

        k_d_exp = self.sample_value(self.range_k_d_exp)
        k_d = round_sig(10 ** k_d_exp / 60)

        return [k_on, k_off, np.nan, k_syn, k_d]

    def generate_and_write_strategies(self, nr_of_alleles, k_on_first=True):

        columns = ['name', 'k_on', 'k_off', 'coord_group', 'k_syn', 'k_d']

        data = []

        for i in range(0, nr_of_alleles):
            self.counter = self.counter + 1
            name = "generated_" + str(self.counter)

            [k_on, k_off, np.nan, k_syn, k_d] = self.create_random(k_on_first)

            while not self.test_restrictions(k_on, k_off, k_syn, k_d):
                [k_on, k_off, np.nan, k_syn, k_d] = self.create_random(k_on_first)

            data.append([name] + [k_on, k_off, np.nan, k_syn, k_d])

        df_strategies = pd.DataFrame(data=data, columns=columns)
        df_strategies.to_csv(path_or_buf=self.filename, sep=';', index=False)

    # should return True when passing tests
    @staticmethod
    def test_restrictions(k_on, k_off, k_syn, k_d) -> bool:

        # hour production when assuming no decay
        hour_production = k_syn * 60 * k_on/(k_on + k_off)

        # steady state approximation (valid for high k_syn)
        # mean_ssa = (k_syn/k_d) * k_on/(k_on + k_off)
        # steady state approximation if the burst takes long enough to reach it
        mean_ssa = k_syn/k_d

        ret_val = ((hour_production > 1) & (mean_ssa < 100))

        return ret_val

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

    # generate a mixed set of strategies with equal fraction ON and OFF for combination of:
    # - stochastic and fluctuating
    # - coordinated/non-coordinated
    # - length of periods (range between 1h to 24h)
    # - synthesis rates (1 or 0.1 per minute)
    # - half lives (1 to 5 hours)
    # - number of alleles in cluster
    def generate_mixed_set(self):

        columns = ['name', 'k_on', 'k_off', 'coord_group', 'k_syn', 'k_d', 'tran_type']

        data = []

        # 6 increasing periods
        periods_in_hours = [1, 2, 3, 5, 12, 24]

        tran_types = ['S', 'F']
        coord_types = [True, False]

        half_lives_in_hours = [1, 5, 2, 4, 3]
        k_syns_in_minutes = [0.1, 1]

        k_d_in_minutes = [np.log(2) / (60 * h) for h in half_lives_in_hours]

        coord_group_counter = 0

        for period_h in periods_in_hours:
            logging.info("start period_h {}".format(period_h))

            # calculate k_on and k_off
            one_state_m = 60 * period_h / 2
            k_on = round(1 / one_state_m, 4)
            k_off = round(1 / one_state_m, 4)

            categories = product(tran_types, coord_types)
            for tran_type, sync_type in categories:
                logging.info("start category{} {}".format(tran_type, sync_type))
                nrs_in_cluster = [1, 2, 4, 8, 32]

                for nr_in_cluster in nrs_in_cluster:

                    # each cluster is either coordinated or not
                    # determine coord_group
                    if sync_type:
                        coord_group_counter = coord_group_counter + 1
                        coord_group = coord_group_counter
                        coord_label = "C"
                    else:
                        coord_group = np.nan
                        coord_label = "NC"
                    for allele in range(0, nr_in_cluster):
                        self.counter = self.counter + 1
                        name = "a_{}_{}_{}_{}".format(self.counter, tran_type, coord_label, period_h)

                        # pick k_syn and pick k_d
                        k_syn = np.random.choice(k_syns_in_minutes)
                        k_d = round(np.random.choice(k_d_in_minutes), 4)

                        # add generated allele/strategy
                        data.append([name, k_on, k_off, coord_group, k_syn, k_d, tran_type])

        df_strategies = pd.DataFrame(data=data, columns=columns)
        df_strategies.to_csv(path_or_buf=self.filename, sep=';', index=False)
