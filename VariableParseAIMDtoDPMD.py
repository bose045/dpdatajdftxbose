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
import re
from collections import defaultdict
"""
process multi systems
"""

def extract_params_from_filename(filename):
    # Extract values using regex based on the file naming pattern
    match = re.search(r'sys(\d+)_run(\d+)_(\d+)_(\d+).jdftxout', filename)
    if match:
        return tuple(map(int, match.groups()))
    raise ValueError(f"Invalid filename format: {filename}")

def k_fold_sampling(files, k):
    # Initialize folds
    folds = [[] for _ in range(k)]
    
    for file in files:
        w, x, y, z = extract_params_from_filename(file)  # A function to extract parameters from each individual file
        print(f'w,x,y,z = {w,x,y,z}')
        print('ktot = ',k)
        rem = (z+k) % k
        #print('rem',rem)
        for kk in range (k):
            if rem != kk:
                folds[kk].append(file)
                print('goes to kth fold = ',kk)
    # Print the number of files in each fold
    for idx, fold in enumerate(folds):
        print(f"Number of files in fold {idx}: {len(fold)}")
    return folds

if len(sys.argv) < 2:
    print("Usage: python VariableParseAIMDtoDPMD.py <int for step stride to take> <K = Rejct Kth fold (within 1 to Ktot) optional> <Shortrun True/False optional> <checkConvergence optional> conversionInfo")
    #eg: python /home/pmech/git_repos/dpdatajdftxbose/VariableParseAIMDtoDPMD.py 10 3 True False > log
    exit(1)
    
# add rejection of initial data

step = int(sys.argv[1])
Ktot = int(sys.argv[2]) if (len(sys.argv) > 2) else 1  # default to 1
shortrun = sys.argv[3] if (len(sys.argv) > 3) else False  # default to 1
print('shortrun',shortrun)
convergence_check = sys.argv[4] if (len(sys.argv) > 4) else True  # default to True

if convergence_check == 'False':
    convergence_check = False

# convFp = open('conversionInfo',"w")
# all_files = glob('./*.jdftxout')  # remember to change here !!!   
all_files = [os.path.abspath(file) for file in glob('./*.jdftxout')]
folds = k_fold_sampling(all_files, Ktot)

for K in range(1,Ktot+1):
    print(f' {step=}')
    print(f' {Ktot=}')
    print(f' {K=}') 
    #dirs = ['train','val']
    dirs = ['train']
    for curDir in dirs:
        os.chdir(curDir)
        if shortrun:
            fs =folds[K-1]
        else:
            fs = all_files
            
        #fs=sorted(glob('./*.jdftxout'))  # remember to change here !!!
        print(fs)
        ms=MultiSystems()
        ls=[]
        for f in fs:
            # break
            #print(f)
            try:
                # ls=LabeledSystem(f, format='jdftxout',step=1)
                # ls=LabeledSystem(f, format='jdftxout',step=step, convergence_check=convergence_check)

                ls=LabeledSystem(f, format='jdftxout',step=step, convergence_check=convergence_check, Ktot=Ktot, K=K, shortrun=shortrun)
                #ls=LabeledSystem(f)
                print(ls)
            except:
                
                print("File",f,"failed to read")
                
            if len(ls)>0:
                ms.append(ls)
        #ls.to('vasp/poscar', 'POSCAR', frame_idx=0)
        ms.to_deepmd_raw('K_'+str(K))
        ms.to_deepmd_npy('K_'+str(K))
        
    os.chdir('..')
    print('Current directory is: ', os.getcwd())

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
