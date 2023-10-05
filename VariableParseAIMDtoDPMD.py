# %%
# # convert multiple systems (w diff atom sizes) outcar from our jdftx to outcar script and prepoutcar sh to the deepmd format
# conda activate deepmd2
#
# import os
# os.chdir('/home/kamron/dpdatajdftx_testing')

# import sys
# sys.path.insert(0, '../dpdatajdftx/')  # to allow testing local python packages 
from dpdata import LabeledSystem,MultiSystems
from glob import glob
import numpy as np
import os
import sys

"""
process multi systems
"""


if len(sys.argv) < 2:
    print("Usage: python VariableParseAIMDtoDPMD.py <int for step stride to take> <checkConvergence optional> conversionInfo")
    exit(1)
    
# add rejection of initial data

step = int(sys.argv[1])
convergence_check = sys.argv[2] if (len(sys.argv) > 2) else True  # default to True
if convergence_check == 'False':
    convergence_check = False

print(f' {step=}')

dirs = ['train','val']
# dirs = ['train']
for curDir in dirs:

    os.chdir(curDir)
    # convFp = open('conversionInfo',"w")
    fs=sorted(glob('./*.jdftxout'))  # remember to change here !!!
    ms=MultiSystems()
    ls=[]
    for f in fs:
        # break
        print(f)
        try:
            # ls=LabeledSystem(f, format='jdftxout',step=1)
            ls=LabeledSystem(f, format='jdftxout',step=step, convergence_check=convergence_check)
            #ls=LabeledSystem(f)
            print(ls)
        except:
            
            print("Fail to read^")
            
        if len(ls)>0:
            ms.append(ls)
    ls.to('vasp/poscar', 'POSCAR', frame_idx=0)
    
    
    ms.to_deepmd_raw('data')
    ms.to_deepmd_npy('data')
    
    os.chdir('..')
    
    # %%
    # Test producing outcar to compare to prior one
    

    # ls['coords'][0][0]
    # ls['forces'][0][0]
    # ls['atom_types']
    
        #create file with last major update noted
    
    
    
    # print(f'shape of forces - frams,atoms,forces', file=convFp)
    # print(np.shape(np.stack(ls['forces'])), file=convFp)  # (201, 54, 3)
    
# f.close()
    
# %%
'''
# Compare to outcar values of coord force energy etc to confirm working
fs=glob('./KFtesting/*.outcar')  # remember to change here !!!
ms=MultiSystems()
ls=[]
for f in fs:
    # break
    # try:
    ls=LabeledSystem(f, format='outcar')
'''
# %%
