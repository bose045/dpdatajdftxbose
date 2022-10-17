import numpy as np
import re
import warnings
import sys
import gzip


def get_frames(fname, begin = 0, step = 10, ml = False, convergence_check=True, type_idx_zero = True):
    #Units:
    eV = 1./27.2114  # divide by eV to go from jdftx H to eV
    Angstrom = 1/0.5291772
    
    fp = open(fname)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []    
    
    nSteps = 0 #number of processed steps
    nEvery = step #select this many frames
    stepActive = False #Whether to process current data
    latvecActive = False #Whether reading lattice vectors
    stressActive = False #Whether reading stress tensor
    atposActive = False #Whether reading atomic positions
    forcesActive = False #Whether reading forces

    for iLine,line in enumerate(fp):
        if line.find('total atoms') > 0:
            atTotalNumb = int(line.split()[4])  # total atoms
        if line.startswith('IonicDynamics: Step:'):
        # if line.startswith('IonicMinimize: Iter:'):
            tokens = line.split()
            iStep = int(tokens[2])
            stepActive = (iStep % nEvery == 0)
            PE_tot = float(tokens[4])/eV
        #Lattice vectors:
        if latvecActive and iLine<refLine+3:
            iRow = iLine-refLine
            R[iRow] = [ float(tok)/Angstrom for tok in line.split()[1:-1] ]
            if iRow==2:
                latvecActive = False
        if line.startswith('R ='):
            latvecActive = True
            refLine = iLine+1
            R = np.zeros((3,3))
        #Stress tensor:
        if stressActive and iLine<refLine+3:
            iRow = iLine-refLine
            stress[iRow] = [ float(tok)/(eV/Angstrom**3) for tok in line.split()[1:-1] ]
            if iRow==2:
                stressActive = False
        if stepActive and line.startswith('# Stress tensor in'):
            stressActive = True
            refLine = iLine+1
            stress = np.zeros((3,3))
        #Atomic positions:
        if atposActive and iLine<refLine+atTotalNumb:
            iRow = iLine-refLine
            tokens = line.split()
            atom_names.append(tokens[1])
            atpos[iRow] = [ float(tok) for tok in tokens[2:5] ]
            if iRow+1==atTotalNumb:
                atposActive = False
                if coordsType == "cartesian":
                    atpos *= 1./Angstrom
                else:
                    atpos = np.dot(atpos, R.T) #convert to Cartesian (Angstrom)
                atom_names = np.array(atom_names)
        if stepActive and line.startswith('# Ionic positions in '):
            atposActive = True
            refLine = iLine+1
            atpos = np.zeros((atTotalNumb,3))
            atom_names = []
            coordsType = line.split()[4]
        #Forces:
        if forcesActive and iLine<refLine+atTotalNumb:
            iRow = iLine-refLine
            tokens = line.split()
            forces[iRow] = [ float(tok) for tok in tokens[2:5] ]
            if iRow+1==atTotalNumb:
                forcesActive = False
                if coordsType == "Cartesian":
                    forces *= 1./(eV/Angstrom)
                else:
                    forces = np.dot(forces, np.linalg.inv(R)/eV) #convert to Cartesian (eV/Angstrom)
        if stepActive and line.startswith('# Forces in '):
            forcesActive = True
            refLine = iLine+1
            forces = np.zeros((atTotalNumb,3))
            coordsType = line.split()[3]
        #Energy components:
        if stepActive and line.startswith("# Energy components:"):

            #Not actually reading energy components at the moment (just calc/reporting PE)
            #Frame complete: write to OUTCAR
            # NumAt = [np.sum(np.char.count(atNames, 'Na')),
            #     np.sum(np.char.count(atNames, 'Mg')), 
            #     np.sum(np.char.count(atNames, 'Cl')) ]
            # PE = PE_tot - (E_NaPlus * NumAt[0] + E_MgPlusPlus * NumAt[1] + E_ClMinus * NumAt[2])/eV
            energy = PE_tot
            # accumulate items after each step
            all_coords.append(atpos)
            all_cells.append(R)
            all_energies.append(energy)
            all_forces.append(forces)

    fp.close()

    unique_atom_names = np.unique(atom_names)

    ions_per_type = []
    for atom_name in unique_atom_names:
        ions_per_type.append(len(np.where(atom_names == atom_name)[0]))

    atom_types = []
    for idx,ii in enumerate(ions_per_type) :
        for jj in range(ii) :
            if type_idx_zero :
                atom_types.append(idx)
            else :
                atom_types.append(idx+1)

    return unique_atom_names, ions_per_type, atom_types, np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces)

# fname = './KFtesting/CCPBED3md0112978.jdftxout'
# begin = 0
# step = 10 
# ml = False
# type_idx_zero = True