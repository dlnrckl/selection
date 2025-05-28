
import subprocess
import sys
import os

iteration = sys.argv[1]

#subprocess.call('cd ./' + dom_sel + '/run' + iteration, shell = True)
os.chdir('./Body_mass_index/run_' + iteration)
input_list = subprocess.run('ls *.ms', capture_output=True, text=True, shell = True).stdout.strip().split('\n')
print(input_list)

for inf in input_list:
	sampling_time = inf.split('.')[0].split('_')[-1]
	#epoch = inf.split('.')[0].split('_')[1]
	if sampling_time == "1000":
		sampling_time2 = "Neolithic"
	elif sampling_time == "2000":
		sampling_time2 = "Late_Calcolithic"
	else:
		sampling_time2 = "Modern"

	print('sampling_time: ', sampling_time2)
	subprocess.call('python3 /mnt/NEOGENE1/projects/selection_2023/dilanur/slim/lc_m_dec_slim_scripts/outputs/estimate_summary.py ' + inf + ' ' + iteration + ' ' + sampling_time, shell = True) 


####################################################################
####################################################################
