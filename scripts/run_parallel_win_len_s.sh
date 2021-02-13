# example for running for all window lengths in win_lens.txt with efficiency 50%:
#	run_parallel_win_len_s.sh 0.5
#conda deactivate
#conda activate python37

#win_lens.txt contains all window lengths
par_input=win_lens.txt
efficiency=$1

mutant=$(hostname | tr '.' '\t' | cut -f1)
out_dir=/hosts/linuxhome/${mutant}/tmp/richard/sc_runs_${efficiency}

parallel -a $par_input -j 8 bash run_main_for_win_len.sh {} $efficiency $out_dir