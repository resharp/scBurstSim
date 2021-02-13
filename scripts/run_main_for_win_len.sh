# this script can be run alone but is intented to be run from run_parallel_win_len_s.sh
# example for running with window of length 15, 50% efficiency and output folder sc_runs_15
#
# run_main_for_win_len 15 0.5 sc_runs_15
#
window="${1//[$'\t\r\n ']}"
efficiency=$2
out_dir=$3
date_start=$(date)
mutant=$(hostname | tr '.' '\t' | cut -f1)
log_line="${date_start} Running main for window length ${window} on ${mutant}"
echo ${log_line}

source_dir=../..

echo $source_dir
echo "output in" $out_dir 
mkdir $out_dir
ls $source_dir/scBurstSim/main.py
echo -w $window  -o $out_dir

# strategies_mixed_new.csv also contains strategies with transcript half life of 30 minutes
strategy_file=~/sc_runs/strategies_mixed_new.csv

python3 $source_dir/scBurstSim/main.py -nc 500 -sf $strategy_file -g 0 -w $window -o $out_dir -e ${efficiency}
