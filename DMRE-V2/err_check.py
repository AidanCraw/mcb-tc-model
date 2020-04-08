import sys
import os.path
from os import path

fileroot = '/n/holyscratch01/keith_lab_seas/acrawford/geomod/data/Run/archive/'
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

print('\nCases with BuildError:\n')

for i in range(first_case,last_case+1):
    docn_tag = '{:03d}'.format(i)
    fullPath = fileroot+case_nameroot+docn_tag
    
    if path.exists('{0}'.format(fullPath)) is False:
        print(case_nameroot+docn_tag)
print('')
