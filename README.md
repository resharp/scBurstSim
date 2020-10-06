# scBurstSim
Single cell transcriptional bursts and noise simulator

``` 
usage: main.py [-h] -nc [number of cells] [-e [efficiency of RNA retrieval]]
               [-o [out_dir]] -sf [strategies_file]

optional arguments:
  -h, --help            show this help message and exit
  -nc [number of cells], --nr_cells [number of cells]
                        Nr of cells for which to run simulation
  -e [efficiency of RNA retrieval], --efficiency [efficiency of RNA retrieval]
                        Efficiency of RNA retrieval on single cell level
                        (default 0.1)
  -o [out_dir], --out_dir [out_dir]
                        Output directory for scBurstSim
  -sf [strategies_file], --strategies_file [strategies_file]
                        Strategies file with burst parameters for alleles
```
                        