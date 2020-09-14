from simulator.Experiment import *
from simulator.Transcription import *

max_minutes = 1440 # 24 hours = 1440 minutes
# part of Experiment (Uridine analog windows)
windows = [[400, 520, 'EU']] # e.g. 120 minutes of EU labeling
WINDOW_START = 0; WINDOW_END = 1; WINDOW_LABEL = 2
freeze = 490

params = TranscriptParams(l_01=0.02, l_10=0.02, k_syn=0.16, nr_refractions=1, k_d=0.01)

# Let's do an experiment
nr_cells = 10
nr_alleles = 10
exp = Experiment(nr_cells, nr_alleles, params, windows, freeze)

df_counts = exp.run()

print("Number of counts: {counts}.".format(counts=len(df_counts)))

debug = "True"