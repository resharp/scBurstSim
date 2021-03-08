# this script is intented to be run from run_parallel_win_gap_len_s.sh
#
window=$(echo $1|cut -f1 -d ';')
gap=$(echo $1|cut -f2 -d ';')
efficiency=$2
out_dir=$3
date_start=$(date)
mutant=$(hostname | tr '.' '\t' | cut -f1)
log_line="${date_start} Running main for window length ${window} and ${gap} on ${mutant}"
echo ${log_line}

source_dir=../..
#echo $source_dir

echo "output in" $out_dir 
# mkdir $out_dir
# ls $source_dir/scBurstSim/main.py

# echo window ':' $window
# echo gap ':' $gap
# echo -w $window -g $gap -o $out_dir

# strategies_mixed_new.csv also contains strategies with transcript half life of 30 minutes
# strategies_mixed_new_2.csv is a sub set only containing the genes with 120 minute periods
strategy_file=~/sc_runs/strategies_mixed_new_2.csv

python3 $source_dir/scBurstSim/main.py -nc 500 -sf $strategy_file -g $gap -w $window -o $out_dir -e ${efficiency}
