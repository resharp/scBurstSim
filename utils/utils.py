
import numpy as np


def round_sig(x, n=4):

    if x < 0:
        x = np.abs(x)
        round_to_n = round(x, -int(np.floor(np.log10(x))) + (n - 1))
        return -round_to_n
    round_to_n = round(x, -int(np.floor(np.log10(x))) + (n - 1))
    return round_to_n
