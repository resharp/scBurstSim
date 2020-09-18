import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp
from sklearn.linear_model import LogisticRegression
import pandas as pd
import seaborn as sns


def show_distribution_real_counts(df_counts, nr_cells):
    # NB: not all zeroes are shown! (when burst was absent; however, can be derived from unlabeled counts)
    plt.title("Distribution of real labeled mRNA counts over {nr} single cells".format(nr=nr_cells))
    plt.xlabel("count labeled transcripts")
    plt.ylabel("occurence")
    plt.hist(df_counts["real_count"], bins=df_counts["real_count"].max())
    plt.show()


# to do: fix: this logistic regression does not work well on inbalanced ON/OFF dataset
def try_out_logistic_regression(perc, df_counts_label):

    # instantiate the model (using the default parameters)
    logreg = LogisticRegression()

    x_train = np.nan_to_num(df_counts_label[["fraction"]].values)

    df_counts_label.loc[df_counts_label["state_{perc}".format(perc=perc)] == "OFF", "state"] = 0
    df_counts_label.loc[df_counts_label["state_{perc}".format(perc=perc)] == "ON", "state"] = 1

    y_train = df_counts_label[["state"]].values.ravel()

    # fit the model with data
    # to do: split training and validation data set and determine AUC etc
    logreg.fit(x_train, y_train)

    # sanity check: good separation between OFF and ON separation?
    x_test = [[0.0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.7], [0.8], [0.9], [1.0]]
    y_pred = logreg.predict(x_test)

    for i in range(len(x_test)):
        print("{} -> {}".format(x_test[i][0], y_pred[i]))


# goal: visually inspect ON and OFF for a certain percentage ON
def violin_plot_fraction(perc_num, perc, df_counts_eu):

    measure = "state_{perc}".format(perc=perc)
    df_counts_eu[measure] = np.where(df_counts_eu.perc_label_on > perc_num, "ON", "OFF")

    df_counts_eu_on = df_counts_eu[df_counts_eu[measure] == "ON"]["fraction"]
    df_counts_eu_off = df_counts_eu[df_counts_eu[measure] == "OFF"]["fraction"]

    ks_result = ks_2samp(df_counts_eu_on, df_counts_eu_off)

    minus_log_pvalue = -np.log10(ks_result.pvalue).round(2)

    sns.violinplot(data=df_counts_eu, y="fraction", x=measure, inner="stick")
    plt.ylim(0, None)
    plt.title("Labeled EU transcript fraction with OFF <= {perc}% burst or ON > {perc}% burst during window"
              "; -log10(pvalue)={pvalue}".
              format(perc=perc, pvalue=minus_log_pvalue))
    plt.show()

    return df_counts_eu


# goal: see what cut off of percentage ON bests separates the two distributions
def do_kolmogorov_smirnov_tests_for_percentages_on(df_counts_eu):
    ks_results = []

    for i in range(10):
        j = i + 1

        perc_num = j/10
        perc = str(j*10)

        measure = "state_{perc}".format(perc=perc)
        df_counts_eu[measure] = np.where(df_counts_eu.perc_label_on > perc_num, "ON", "OFF")

        df_counts_eu_on = df_counts_eu[df_counts_eu[measure] == "ON"]["fraction"]
        df_counts_eu_off = df_counts_eu[df_counts_eu[measure] == "OFF"]["fraction"]

        if len(df_counts_eu_off > 0) and len(df_counts_eu_on > 0):
            ks_result = ks_2samp(df_counts_eu_on, df_counts_eu_off)

            minus_log_pvalue = -np.log10(ks_result.pvalue).round(2)

            ks_results.append([int(perc), minus_log_pvalue])

    df_ks = pd.DataFrame(ks_results, columns=["perc", "minus_log_pvalue"])
    plt.plot(df_ks.perc, df_ks.minus_log_pvalue)
    plt.show()

def regression_plot(x, y, data):

    # g0 = sns.lmplot(x=x, y=y,
    #                 # maybe
    #                 # hue="noise_level",
    #                 # hue_order=["level_1", "level_2", "level_3"],
    #                 data=data,
    #                 height=5, aspect=1.5)

    sns.jointplot(x=x, y=y, data=data, kind="reg",
                  marginal_kws=dict(bins=20, rug=True))

    plt.xlim(-0.1, 1.1)
    if y == "fraction":
        plt.ylim(-0.1, 1.1)

    # title = "{x} agains {y}".format(x=x, y=y)
    # plt.title(title)

    plt.show()

