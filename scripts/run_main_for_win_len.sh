window="${1//[$'\t\r\n ']}"
date_start=$(date)
mutant=$(hostname | tr '.' '\t' | cut -f1)
log_line="${date_start} Running main for window length ${window} on ${mutant}"
echo ${log_line}

source_dir=../..

#out_dir=sc_runs_${window}
out_dir=/hosts/linuxhome/mutant1/tmp/richard/sc_runs_${window}

echo $source_dir
echo "output in" $out_dir 
mkdir $out_dir
ls $source_dir/scBurstSim/main.py
echo -w $window  -o $out_dir
python3 $source_dir/scBurstSim/main.py -nc 500 -sf ~/sc_runs/strategies_mixed.csv -g 0 -w $window -o $out_dir -e 1

#sleep $window
