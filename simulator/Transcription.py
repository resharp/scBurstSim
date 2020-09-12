import numpy as np

WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2


class Transcription:
    l_01 = 0
    l_10 = 0
    k_syn = 0  # k_syn = synthesis rate = transcription rate
    nr_refractions = 1
    # decay is strictly not part of transcription but we include it in the model
    k_d = 0.01  # k_d = decay rate

    def __init__(self, l_01, l_10, k_syn, nr_refractions, k_d):
        self.l_01 = l_01
        self.l_10 = l_10
        self.k_syn = k_syn
        self.k_d = k_d

    def new_poisson_arrivals(self, start_time, interval, windows = []) -> list:
        poisson_list = []

        last_arrival = 0
        label = ""
        while last_arrival < interval:
            arrival = np.random.exponential(scale=1.0, size=None) / self.k_syn
            last_arrival = last_arrival + arrival

            # we immediately determine the decay time once a transcript comes into existence
            # to do: check if this is the right distribution for a death process
            decay = np.random.exponential(scale=1.0, size=None) / self.k_d
            decay_time = start_time + last_arrival + decay

            arrival_time = start_time + last_arrival

            if last_arrival < interval:
                for window in windows:
                    if window[WINDOW_START] < arrival_time < window[WINDOW_END]:
                        label = window[WINDOW_LABEL]
                    else:
                        label = ""
                poisson_list.append([label, arrival_time, 1, decay_time, -1])

        return poisson_list
