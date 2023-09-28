import numpy as np
import re
import warnings
import sys
import gzip
import glob
from jdftxout import get_frames
import os

# python /home/kamron/dpdatajdftx/dpdata/jdftx/jdftxoutToXYZ.py test.xyz
# python ~/dpdatajdftx/dpdata/jdftx/jdftxoutToXYZ.py 500deg300rot.xyz > conversionInfo


def main(xyzOut):
    pattern = r'*jdftxout'
    with open(xyzOut, 'w') as f:
        for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
            fname = os.path.basename(pathAndFilename)
            print('working on: ', fname)
            unique_atom_names, ions_per_type, atom_types, all_cells, all_coords, all_energies, all_forces = get_frames(fname, step = 10)
            
            # print('Update: working on step ')
            # BUILD xyz output
            AtNumTotal = np.sum(ions_per_type)
            
            # loop over each frame
            for frameNum in range(len(all_energies)):
                f.write(f'{AtNumTotal} \n')
                lattice = ' '.join(['{:.3f}'.format(x) for x in all_cells[frameNum].flatten()])

                # print(lattice)
                # print(all_cells[frameNum])
                # sys.exit(1)
                f.write(f'Lattice="{lattice}" ')
                f.write(f'Properties=species:S:1:pos:R:3:forces:R:3 energy={all_energies[frameNum]} \n')  # add energy!!

                for atIdx in range(AtNumTotal):
                    forces = ' '.join(['{:10.6f}'.format(x) for x in all_forces[frameNum][atIdx].flatten()])
                    coords = ' '.join(['{:10.6f}'.format(x) for x in all_coords[frameNum][atIdx].flatten()])
                    f.write(f'{unique_atom_names[atom_types[atIdx]]} {coords} {forces} \n')

                # f.write('FORCE: {:10.6f} ...  ENERGY: {:16.8f}\n'.format(all_forces[all_forces],energies[i]))
                # f.write('FORCE: {:10.6f} {:10.6f} {:10.6f} ...  ENERGY: {:16.8f}\n'.format(*tuple(all_forces[frameNum]),all_energies[frameNum]))
                # print(all_forces[frameNum][])
                # sys.exit(1)
                
                # f.write('FORCE: {:10.6f} {:10.6f} {:10.6f} ...  ENERGY: \n'.format(*tuple([1,2,3])))
                # loop over each atom
                
                    
                

# Nima Ex 
# 196
# Lattice="-39.110339 -0.016635 1.119012 0.0 -8.676835 0.01733 0.0 0.0 9.284299" 
# Properties=species:S:1:pos:R:3:forces:R:3:magmoms:R:1 energy=-945.97137432 
# stress="0.0029958369992962018 0.0005261467362937069 -0.0031765099639631447 0.0005261467362937069 0.0022300537616233326 -0.0013310018210946047 -0.0031765099639631447 -0.0013310018210946047 0.006073050794575669" 
# magmom=4.06e-05 free_energy=-945.97137432 pbc="T T T"
# H       -3.69440000      -1.04475000       6.12938000       0.35695600      -0.04549100       0.57023800       0.00000000
if __name__ == "__main__":
    xyzOut = sys.argv[1]
    main(xyzOut)