#conda deactivate
#conda activate python37
par_input=win_lens.txt
parallel -a $par_input -j 8 bash run_main_for_win_len.sh {}