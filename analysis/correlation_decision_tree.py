# predict transaction type fluctuating (F)/ stochastic (S) based on:
#   Pearson correlation for different pulse lengths
#   for a certain efficiency of RNA retrieval
#   assuming resolution on allele level
# based on df_corr_all
# train on 70% of data
# irrespective of period

# inspired on
# https://www.datacamp.com/community/tutorials/decision-tree-classification-python
import io
import os
import pandas as pd
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.model_selection import train_test_split
from sklearn import metrics
from IPython.display import Image
import pydotplus

if os.name == 'nt':
    dir_sep = "\\"
    prj_dir = r"D:\26 Battich Oudenaarden transcriptional bursts"
else:
    dir_sep = "/"
    prj_dir = "sc_runs"

gap = 0
label_1 = "EU"
label_2 = "4SU"

eff = 1

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

# quick and dirty: replace np.nans with Pearson correlation 0
data = data.fillna(0)

# feature_cols = ['corr_15', 'corr_30', 'corr_45', 'corr_60', 'corr_75', 'corr_90',
#                 'corr_105', 'corr_120', 'corr_135', 'corr_150', 'corr_165', 'corr_180',
#                 'corr_195']
# manual selection of features:
feature_cols = ['corr_15', 'corr_30', 'corr_45', 'corr_60', 'corr_90', 'corr_120']

data["label"] = 0
data.loc[data.tran_type == 'F', "label"] = 1

accuracies = []

for i in range(0, 5):
    X = data[feature_cols]
    y = data.label

    # 70% training and 30% test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1)

    # Create Decision Tree classifier object (uses gini as a default)
    # clf = DecisionTreeClassifier()
    clf = DecisionTreeClassifier(criterion="gini", max_depth=4)

    # Train Decision Tree Classifier
    clf = clf.fit(X_train, y_train)

    # Predict the response for test dataset
    y_pred = clf.predict(X_test)

    accuracy = metrics.accuracy_score(y_test, y_pred).round(3)

    print("Accuracy:", accuracy)
    accuracies.append(accuracy)

print("min accuracy for eff {}: {}".format(eff, min(accuracies)))
print("max                  {}: {}".format(eff, max(accuracies)))


# now visualize the last decision tree
plot_dir = r"{}{}runs{}{}".format(prj_dir, dir_sep, dir_sep, "correlation_labels.plots")

# dirty trick to find the Graphviz executables
os.environ['PATH'] = os.environ['PATH']+';' + r"C:\Program Files\Graphviz\bin"

dot_data = io.StringIO()
export_graphviz(clf, out_file=dot_data,
                filled=True, rounded=True,
                special_characters=True, feature_names=feature_cols, class_names=['S', 'F'])
graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
graph.write_png(plot_dir + dir_sep + 'decision_tree.png')
Image(graph.create_png())

