

# download runs for windows in range 15 - 120 to current directory
sshpass -f sshpass.txt scp richard@alive.bio.uu.nl:/hosts/linuxhome/mutant1/tmp/richard/sc_runs_*/df_counts_W*_G0.csv .

sshpass -f sshpass.txt scp richard@alive.bio.uu.nl:/hosts/linuxhome/mutant11/tmp/richard/sc_runs_*/df_counts_W*_G0.csv .
sshpass -f sshpass.txt scp richard@alive.bio.uu.nl:/hosts/linuxhome/mutant36/tmp/richard/sc_runs_*/df_counts_W*_G0.csv .
sshpass -f sshpass.txt scp richard@alive.bio.uu.nl:/hosts/linuxhome/mutant10/tmp/richard/sc_runs_*/df_counts_W*_G0.csv .
# /hosts/linuxhome/mutant1/tmp/richard/sc_runs_1
# /hosts/linuxhome/mutant11/tmp/richard/sc_runs_0.5
# /hosts/linuxhome/mutant36/tmp/richard/sc_runs_0.2
# /hosts/linuxhome/mutant10/tmp/richard/sc_runs_0.05