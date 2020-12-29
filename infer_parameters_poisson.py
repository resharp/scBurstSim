# here we fit a Poisson distribution to dummy data and try to find a global minimum
# by changing inital conditions

from scipy.stats import poisson
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from utils.utils import round_sig


# simplified model
def poisson_model(x, mu):

    return poisson.pmf(x, mu)


def generate_dummy_data(mu):
    len_interval = 3 * mu

    x = [r for r in range(0, len_interval)]
    y = [0, 2, 3, 4, 5, 6, 5, 3, 2, 2, 3, 2, 1, 2, 3, 3, 3, 3, 2, 3, 2, 2, 1, 0, 1]
    len_y = len(y)
    y = y + [0]*(len_interval-len_y)

    return x, y


mu = 20
x, y = generate_dummy_data(mu)

y_norm = [n/sum(y) for n in y]


expected_guesses = [int(mu / 4), int(mu / 2), mu]

for expected in expected_guesses:

    popt, pcov = curve_fit(poisson_model, x, y_norm, expected)

    estimated_mean = round_sig(popt[0], 3)

    plt.title("initial guess={expected}; estimated mean={mean}".format(expected=expected, mean=estimated_mean))
    plt.plot(x, poisson.pmf(x, estimated_mean), label="Poisson")
    plt.plot(x, y_norm, label="real")
    plt.legend()
    plt.show()
    plt.close(1)

expected_guesses = range(0, 20, 2)
cov = []
means = []
for expected in expected_guesses:
    popt, pcov = curve_fit(poisson_model, x, y_norm, expected)
    estimated_mean = round_sig(popt[0], 3)
    cov.append(pcov[0])
    means.append(estimated_mean)

plt.title("find fitting optimum by changing initial guesses")
plt.plot(expected_guesses, cov)
plt.xlabel("initial guess")
plt.ylabel("variance")
plt.show()
plt.close(1)

