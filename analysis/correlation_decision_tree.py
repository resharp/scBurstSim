# predict transaction type fluctuating (F)/ stochastic (S) based on:
#   Pearson correlation for different pulse lengths
#   for a certain efficiency of RNA retrieval
#   assuming resolution on allele level
# based on df_corr_all
# train on 70% of data
# irrespective of period

# inspired on
# https://www.datacamp.com/community/tutorials/decision-tree-classification-python
import os
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics


if os.name == 'nt':
    dir_sep = "\\"
    prj_dir = r"D:\26 Battich Oudenaarden transcriptional bursts"
else:
    dir_sep = "/"
    prj_dir = "sc_runs"

gap = 0
label_1 = "EU"
label_2 = "4SU"

eff = 0.2
out_dir = r"{}{}runs_on_server_{}".format(prj_dir, dir_sep, eff)
plot_dir = out_dir + dir_sep + "correlation_labels.plots"

corr_name = "{od}{dir_sep}df_corr_all_G{gap}.csv".format(
    od=plot_dir, dir_sep=dir_sep, gap=gap)

df_corr_all = pd.read_csv(corr_name, sep=';')

# now we want to convert into feature columns
# we may have missing values for pulse lengths for which correlation could not be calculated
# https://stats.stackexchange.com/questions/96025/how-do-decision-tree-learning-algorithms-deal-with-missing-values-under-the-hoo


df_corr_sign = df_corr_all[df_corr_all.p_value < 0.05]

df_corr_sign = df_corr_sign[['strategy', 'tran_type', 'corr', 'len_win']]

data = df_corr_sign.set_index(['strategy', 'tran_type', 'len_win'])

data = data.unstack()
data.columns = ["corr_" + str(x[1]) for x in data.columns.ravel()]
data = data.reset_index()

# quick and dirty: replace np.nans with 0
data = data.fillna(0)

feature_cols = ['corr_15', 'corr_30', 'corr_45', 'corr_60', 'corr_75', 'corr_90',
                'corr_105', 'corr_120', 'corr_135', 'corr_150', 'corr_165', 'corr_180',
                'corr_195']

data["label"] = 0
data.loc[data.tran_type == 'F', "label"] = 1

X = data[feature_cols]
y = data.label

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1)  # 70% training and 30% test

# Create Decision Tree classifier object (uses gini as a default)
clf = DecisionTreeClassifier()

# Train Decision Tree Classifier
clf = clf.fit(X_train, y_train)

# Predict the response for test dataset
y_pred = clf.predict(X_test)


print("Accuracy:", metrics.accuracy_score(y_test, y_pred))

debug = True
