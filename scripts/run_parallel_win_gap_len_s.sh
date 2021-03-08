# example for running for all window lengths in win_lens.txt with efficiency 20%:
#	run_parallel_win_gap_len_s.sh 0.2
# conda deactivate
# conda activate python37

# win_gap_lens.txt contains all window lengths with gaps
par_input=win_gap_lens.txt

efficiency=$1
win_len=60

mutant=$(hostname | tr '.' '\t' | cut -f1)
out_dir=/hosts/linuxhome/${mutant}/tmp/richard/sc_runs_${win_len}_${efficiency}
# out_dir=temp

parallel -a $par_input -j 8 bash run_main_for_win_gap_len.sh {} $efficiency $out_dir
# ls $par_input