from simulator.Experiment import *
from simulator.Transcription import *
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp
from sklearn.linear_model import LogisticRegression

max_minutes = 1440  # 24 hours = 1440 minutes
windows = [[400, 520, 'EU']]  # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = 550  # freeze 30 minutes after end of last window

params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)

nr_cells = 1000
nr_alleles = 1
exp = Experiment(nr_cells, nr_alleles, params, windows, freeze)

df_counts = exp.run()

df_counts["fraction"] = df_counts["real_count"] / (df_counts["real_count"] + df_counts["real_count_unlabeled"])

df_counts_eu = df_counts[df_counts.label == "EU"].copy(deep=True)


def violin_plot_fraction(perc_num, perc):

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


df_counts_eu = violin_plot_fraction(0.5, "50")


def do_kolmogor_smirnov_tests_for_percentages_on():
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


do_kolmogor_smirnov_tests_for_percentages_on()

df_all_arrivals = exp.df_all_arrivals


# this logistic regression does not work well on inbalanced ON/OFF dataset
def experiment_with_logistic_regression(perc):

    # instantiate the model (using the default parameters)
    logreg = LogisticRegression()

    x_train = np.nan_to_num(df_counts_eu[["fraction"]].values)

    df_counts_eu.loc[df_counts_eu["state_{perc}".format(perc=perc)] == "OFF", "state"] = 0
    df_counts_eu.loc[df_counts_eu["state_{perc}".format(perc=perc)] == "ON", "state"] = 1

    y_train = df_counts_eu[["state"]].values.ravel()

    # fit the model with data
    logreg.fit(x_train, y_train)

    # sanity check: good separation between OFF and ON separation?
    x_test = [[0.0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.7], [0.8], [0.9], [1.0]]
    y_pred = logreg.predict(x_test)

    for i in range(len(x_test)):
        print("{} -> {}".format(x_test[i][0], y_pred[i]))


experiment_with_logistic_regression(perc="50")


def show_distribution_real_counts():
    # NB: not all zeroes are shown! (when burst was absent; however, can be derived from unlabeled counts)
    plt.title("Distribution of real labeled mRNA counts over {nr} single cells".format(nr=nr_cells))
    plt.xlabel("count labeled transcripts")
    plt.ylabel("occurence")
    plt.hist(df_counts["real_count"], bins=df_counts["real_count"].max())
    plt.show()

    print("Number of counts: {counts}.".format(counts=len(df_counts)))


show_distribution_real_counts()
