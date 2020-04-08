import os
import sys
from datetime import datetime

sst_fileroot = '/n/home06/acrawford/mcb_data/sst_files/'
namelist_fileroot = '/n/home06/acrawford/executables/DMRE-V2/'
startup_fileroot = '/n/home06/acrawford/mcb_data/'
sst_base_filename = 'model_1_'
sst_filenum = 512   #default value
case_nameroot = 'v2_'
ref_case = 'mcb_spinup'
ref_date = '0004-08-01'
compset = 'F2000climo'
resolution = 'f09_f09_mg17'
stop_N = 14
stop_option = 'ndays'
run_type = 'branch'

input_args=sys.argv[1:]
if len(input_args) == 2:
    first_case = int(input_args[0])
    last_case  = int(input_args[1])
elif len(input_args) == 0:
    first_case = 1
    last_case = sst_filenum
else:
    raise ValueError('invalid number of arguments - must be either 0 (defaults) or 2 (first and last case)')
assert first_case <= last_case, 'last case cannot precede first case'

for i in range(first_case,last_case+1):
    docn_tag = '{:03d}'.format(i)
    docn = sst_fileroot+sst_base_filename+docn_tag+'.nc'
    
    os.system('cd $CIMEROOT/scripts && ./create_newcase --case {0} --compset {1} --res {2}'.format((case_nameroot+docn_tag),compset,resolution))
    
    os.system('cd $CIMEROOT/scripts/{0} && ./xmlchange STOP_OPTION={1},STOP_N={2},RUN_TYPE={3}'.format((case_nameroot+docn_tag),stop_option,stop_N,run_type))
    os.system('cd $CIMEROOT/scripts/{0} && ./xmlchange RUN_REFCASE={1},RUN_REFDATE={2},SSTICE_DATA_FILENAME={3}'.format(case_nameroot+docn_tag,ref_case,ref_date,docn))

    timestamp = datetime.now().strftime("%H:%M:%S")
    os.system('echo "   Case: {0}   |   Configured at: {1}"'.format(case_nameroot+docn_tag,timestamp))

    os.system('cd $CIMEROOT/scripts/{0} && ./case.setup'.format(case_nameroot+docn_tag))
    os.system('cp {3}{0}/rest/{1}-00000/* $CESMDATAROOT/Run/{2}/run/'.format(ref_case,ref_date,(case_nameroot+docn_tag),startup_fileroot))

    os.system('cp {0}user_nl_cam $CIMEROOT/scripts/{1}/'.format(namelist_fileroot,(case_nameroot+docn_tag)))
    os.system('cd $CIMEROOT/scripts/{0} && ./preview_namelists'.format(case_nameroot+docn_tag))

    timestamp = datetime.now().strftime("%H:%M:%S")
    os.system('echo "   Case: {0}   |   Set up at: {1}"'.format(case_nameroot+docn_tag,timestamp))

    os.system('cd $CIMEROOT/scripts/{0} && ./case.build'.format(case_nameroot+docn_tag))
    
    timestamp = datetime.now().strftime("%H:%M:%S")
    os.system('echo "   Case: {0}   |   Built at: {1}"'.format(case_nameroot+docn_tag,timestamp))
