# %%
# # convert multiple systems (w diff atom sizes) outcar from our jdftx to outcar script and prepoutcar sh to the deepmd format
# conda activate deepmd2
#
from dpdata import LabeledSystem,MultiSystems
from glob import glob
"""
process multi systems
"""
fs=glob('./KFtesting/*.jdftxout')  # remeber to change here !!!
ms=MultiSystems()
ls=[]
for f in fs:
    # break
    # try:
    ls=LabeledSystem(f, format='jdftxout')



#         #ls=LabeledSystem(f)
#     # except:
#         # print(f)
#     if len(ls)>0:
#         ms.append(ls)

# ms.to_deepmd_raw('data')
# ms.to_deepmd_npy('data')
# %%
