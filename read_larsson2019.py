# parse the gene, k_on, k_off and k_syn and mean expression from larsson 2019
# there are three tabs
#
# C57: 6962 alleles
# CAST: 6935 alleles
# Active X: 268 genes
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

input_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\source\scBurstSim\data"
out_dir = r"D:\26 Battich Oudenaarden transcriptional bursts\runs"

if os.name == 'nt':
    dir_sep = "\\"
else:
    dir_sep = "/"

c57_alleles_file = input_dir + dir_sep + "Larsson 2019 C57 allele.csv"
cast_alleles_file = input_dir + dir_sep + "Larsson 2019 CAST allele.csv"


def extract_float_at_index(params_string, index):

    param = 0
    if not pd.isna(params_string):
        params = params_string.replace("[", "").replace("]", "").split()
        param = float(params[index])

    return param


def read_df(file_name, label):

    df_alleles = pd.read_csv(file_name, sep=';',
                             names=["gene", "kon_koff_ksyn", "frequency", "burst_size", "mean_expression"],
                             skiprows=1)

    df_alleles["label"] = label

    df_alleles["k_on"] = df_alleles.kon_koff_ksyn.apply(extract_float_at_index, args=[0])
    df_alleles["k_off"] = df_alleles.kon_koff_ksyn.apply(extract_float_at_index, args=[1])
    df_alleles["k_syn"] = df_alleles.kon_koff_ksyn.apply(extract_float_at_index, args=[2])

    return df_alleles


df_alleles1 = read_df(c57_alleles_file, "C57")
df_alleles2 = read_df(cast_alleles_file, "CAST")
print(len(df_alleles1))
print(len(df_alleles2))

df_all_alleles = pd.concat([df_alleles1, df_alleles2])
print(len(df_all_alleles))

# df_all_alleles[['gene', 'kon_koff_ksyn', 'k_on', 'k_off', 'k_syn']]

df_subset_alleles = df_all_alleles[df_all_alleles["k_on"] < 10]

plt.hist(df_subset_alleles["k_on"], bins=30, log=True)
plt.title("Distribution of relative k_on (in terms of k_d)")
plt.show()

# Read half-lives from Herzog 2017
half_life_file = input_dir + dir_sep + "Herzog 2017 table 2.csv"

df_half_lives = pd.read_csv(half_life_file, sep=';',
                            usecols=[0, 3, 4, 6, 7, 8, 9, 10],
                            # names=["chr", "start", "end", "name", "length", "strand",
                            #        "half_life", "std_half_life", "k_d_h", "std_k_d_h", "rsquare"],
                            names=["chr", "name", "length",
                                   "half_life", "std_half_life", "k_syn_h", "std_k_syn_h", "rsquare"],
                            skiprows=1)

df_half_lives["k_d_calc_h"] = round(np.log(2)/df_half_lives["half_life"], 4)
df_half_lives["k_d_calc_m"] = df_half_lives["k_d_calc_h"] / 60

df_half_lives["half_life_m"] = df_half_lives["half_life"] * 60

df_merge = pd.merge(df_all_alleles, df_half_lives, how="inner",
                    left_on="gene", right_on="name")

# plt.hist(df_half_lives["k_d_calc_m"], bins=100)
# plt.title("Distribution of k_d_calc_m")
# plt.show()

plt.hist(df_half_lives["half_life_m"], bins=100)
plt.title("Distribution of half_life_m (minutes)")
plt.show()

len_intersection = len(df_merge)/2

print("intersection of two data sets Larsson 2019 and Herzog 2017: {len} genes".format(len=len_intersection))

debug = True


