# scBurstSim
Single cell transcriptional bursts and noise simulator

Stochastic modeling of transcription and sampling noise on single cell level. 

# usage 

``` 
usage: main.py [-h] -nc [number of cells] -sf [strategies_file]
               [-e [efficiency of RNA retrieval]] [-o [out_dir]]

optional arguments:
  -h, --help            show this help message and exit
  -nc [number of cells], --nr_cells [number of cells]
                        Nr of cells for which to run simulation
  -sf [strategies_file], --strategies_file [strategies_file]
                        Strategies file with burst parameters for alleles
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

Output: A resulting file df_counts.csv is put in the output directory, when
no output directory is specified it is put in the directory you run your 
script from.

Also some plots for clustering will be placed here (for now only some simple 
clustering based on the fraction of count of one label as compared to the total counts)  

For now, further parameters, e.g. the labeling window(s) can be set by adapting main.py.

# Dependencies 

See file requirements.txt for dependency on Python packages.
 

# Strategies

Strategies file with burst parameters for alleles

An example of the strategies file:

```                        
# named strategies for the simulator/Transcription class (; separated!)
name;k_on;k_off;coord_group;k_syn;k_d
# two strategies with the same burst coordination
frequent_coor;0.02;0.02;1;0.16;0.01
frequent_uncoor;0.02;0.02;;0.16;0.01
frequent_high;0.02;0.02;1;0.8;0.01
large_swing;0.005;0.01;;0.16;0.005
real_bursty_coor;0.005;0.02;2;0.32;0.02
real_bursty_uncoor;0.005;0.02;;0.32;0.02
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

In the above example alleles with the strategy 
frequent_coor and frequent_high are synchronized by a shared coordination group 1.
This is only possible if the transition matrix is identical, so k_on and k_off 
should be identical for coordinated strategies.  


# Example traces 
run single_allele_example.py for plotting example traces of 
specific strategies (see further documentation inline)

single_allele_example.py uses the Transcription class that is also
used by Experiment.py


