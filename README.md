# scBurstSim
Single cell transcriptional bursts/fluctuations and noise simulator

Stochastic modeling of transcription and sampling noise on single cell level. 

# usage 

``` 
usage: main.py [-h] -nc [number of cells] -sf [strategies_file]
               [-g [gap minutes)]] [-w [window length(minutes]]
               [-e [efficiency of RNA retrieval]] [-o [out_dir]]

optional arguments:
  -h, --help            show this help message and exit
  -nc [number of cells], --nr_cells [number of cells]
                        Nr of cells for which to run simulation
  -sf [strategies_file], --strategies_file [strategies_file]
                        Strategies file with burst parameters for alleles
  -g [gap (minutes)], --gap [gap (minutes)]
                        Length of gap (in minutes); default 0
  -w [window length(minutes)], --length_window [window length(minutes)]
                        Length of windows (in minutes); default 60
  -e [efficiency of RNA retrieval], --efficiency [efficiency of RNA retrieval]
                        Efficiency of RNA retrieval on single cell level
                        (default 0.1)
  -o [out_dir], --out_dir [out_dir]
                        Output directory for scBurstSim
```

Example for running from command line
```
python3 "[your source location]/scBurstSim/main.py" -nc 10  -sf "[your path]/strategies.csv" -o "[your output dir]"
```
An example strategies.csv file is in the data folder. See explanation of file below.

A file with a mix of 1128 strategies strategies_mixed.csv is also in the data folder, 
which was used for the analysis of correlation between normalized label counts. 

Output: A resulting file df_counts_W<len_win>_G<len_gap>.csv is put in the output directory, when
no output directory is specified it is put in the directory you run your 
script from.

For now, further parameters, e.g. the (number of) labeling window(s) can be set by adapting main.py.

# Dependencies 

See file requirements.txt for dependency on Python packages.
We used conda environments for handling the dependencies. 

# Strategies

Strategies file with burst parameters for alleles

An example of the strategies file:

```                        
# named strategies for the simulator/Transcription class (; separated!)
name;k_on;k_off;coord_group;k_syn;k_d
# two strategies with the same burst coordination
frequent_coor;0.02;0.02;1;0.16;0.01;S
frequent_high;0.02;0.02;1;0.8;0.01;S
# no coordination: leave coord_group empty
frequent_uncoor;0.02;0.02;;0.16;0.01;S
large_swing;0.005;0.01;;0.16;0.005;S
real_bursty_coor;0.005;0.02;2;0.32;0.02;S
real_bursty_uncoor;0.005;0.02;;0.32;0.02;S
first_example;0.03;0.01;;2;0.02;S
# example fluctuating alleles with fixed period
second_example;0.024;0.03;;1;0.02;F
third_example;0.02;0.02;;2;0.02;S
bimodal;0.01;0.01;;2;0.02;S
powerlaw;0.01;0.04;;2;0.02;S
```
For each strategy a number of alleles will be generated. 
All parameters are per minute.

The first column name is the name of a strategy e.g. for displaying in plots.

k_on and k_off are parameters of the **transition matrix**. k_on is the probability of
switching from a silent (OFF) state to an active (ON) state and k_off is for
switching from active to silent.

k_syn and k_d are the **synthesis (transcription)** and **decay rate** respectively.
k_syn is the rate of the Poisson arrivals during an active state. 
k_d is the decay rate of the transcripts both during the active and the silent state.

If coord_group is not empty and contains a number the alleles are
**synchronized** within each simulated cell, which means that they have the same DTMC
trace (same bursts=same active and silent periods).

Alleles with different k_syn and k_d parameters may be coordinated 
by assigning the same coord_group id to the strategy.   

The last parameter tran_type denotes if there is 
stochastic (S) switching between Markov states 
or deterministic fluctuating (F) behaviour for active/inactive states with fixed period. 

In the above example alleles with the strategy 
frequent_coor and frequent_high are synchronized by a shared coordination group 1.
This is only possible if the transition matrix is identical, so k_on and k_off 
should be identical for coordinated strategies.  


# Examples traces and distributions

## single_allele_example.py 
1. traces: run single_allele_example.py for plotting example traces of specific strategies (see further documentation inline)
2. distributions: quickly create stationary distribution for a single allele example (by sampling snapshots from single trace)

# Analysis

## correlation labels ## 
analysis of correlation between normalized label counts

output: 
- directory correlation_labels.plots containing
- a csv corr_labels_<len_win>.csv with Pearson correlation values
- some plots with correlation values for different categories (tran_type S/F and length of period)


## cluster alleles ##
to do

## analyze parameters

- compare simulated (time-dependent and stationary) distributions 
against theoretical stationary distribution
- examine how quickly the means converge towards the means of the stationary distributions 
(use multiple window lengths)
- make two dimensional scatter and density plot for two labels 
for one strategy (set of dynamical parameters)
- examine correlation between k_d and inferred k_d from means 
from label_1 and label_2 (assuming non-changing parameters)
- extra analysis of chance of being ON in 2nd window give chance ON or OFF in 1st window

## analyze active state
explore the shape of the active state fraction distribution 
(since we miss the theoretical distribution)
- build intuition: create distributions for percentage of active state ("active state fraction distribution")
- explore: try to fit to beta distribution

## infer parameters example
fit to time dependent function of chance of having activity of any length 
during a single labeling window

## infer parameters Poisson
to do


