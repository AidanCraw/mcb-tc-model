import os
import sys
from datetime import datetime

sst_filenum = 512   #default value
case_nameroot = 'v2_'

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

    os.system('cd $CIMEROOT/scripts/{0} && ./case.submit'.format(case_nameroot+docn_tag))

    timestamp = datetime.now().strftime("%H:%M:%S")
    os.system('echo "   Case: {0}   |   Submitted at: {1}"'.format(case_nameroot+docn_tag,timestamp))
