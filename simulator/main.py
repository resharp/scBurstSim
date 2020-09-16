from simulator.Experiment import *
from simulator.Transcription import *
import matplotlib.pyplot as plt
import seaborn as sns

max_minutes = 1440  # 24 hours = 1440 minutes
windows = [[400, 520, 'EU']]  # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = 550  # freeze 30 minutes after end of last window

params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)

nr_cells = 100
nr_alleles = 1
exp = Experiment(nr_cells, nr_alleles, params, windows, freeze)

df_counts = exp.run()
df_counts["state_80"] = np.where((df_counts.label != "") & (df_counts.perc_label_on > 0.8), "ON", "OFF")

df_counts_eu = df_counts[df_counts.label == "EU"]
sns.boxplot(data=df_counts, y="real_count", x="state_80")
sns.swarmplot(data=df_counts, y="real_count", x="state_80", color=".25")
plt.title("EU transcript real counts with OFF < 80% burst or ON > 80% burst")
plt.show()


df_all_arrivals = exp.df_all_arrivals


# NB: not all zeroes are shown! (when burst was absent; however, can be derived from unlabeled counts)
plt.title("Distribution of real labeled mRNA counts over {nr} single cells".format(nr=nr_cells))
plt.xlabel("count labeled transcripts")
plt.ylabel("occurence")
plt.hist(df_counts["real_count"], bins=df_counts["real_count"].max())
plt.show()

print("Number of counts: {counts}.".format(counts=len(df_counts)))
