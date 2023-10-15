import numpy as np
import re
import warnings
import sys
import gzip
import glob


def get_frames(fname, begin = 0, step = 10, ml = False, convergence_check=True, type_idx_zero = True):

    #create file with last major update noted
    with open('VersionLatestFix',"w") as f:
        print(f'alphabetized elements to ensure order',file=f)
        print(f'updated version with fix to append to list now using np.copy!',file=f)
        print(f'added convergence check',file=f)
    
        #Units:    
    eV = 1./27.2114  # divide by eV to go from jdftx H to eV
    Angstrom = 1/0.5291772 # divide by Ang to go from jdftx bohr to ang
    
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
    appliedPotential = False #Whether applied potential is used
    converged = False  # check for convergence before allowing data taken
    # set convergence_check to False to ignore convergence check

    for iLine,line in enumerate(fp):
        # print(line)
        # Check for applied potential
        if line.startswith('ionic-gaussian-potential'):
            appliedPotential = True
            # [] store or correct after each iteration
        # Collect applied potential energy
        if line.startswith('EextIonic'):
            tokens = line.split()
            PE_Applied = float(tokens[2])/eV
            
        if line.find('total atoms') > 0:
            atTotalNumb = int(line.split()[4])  # total atoms
            #initialize size for forces since this doesn't work too well now with the hooks given before the two
            forces = np.zeros((atTotalNumb,3))
            extForces = np.zeros((atTotalNumb,3))
        
        # ElecMinimize: None of the convergence criteria satisfied after 30 iterations.
        # or ElecMinimize: Converged (but not with SCF!)
        if line.startswith('ElecMinimize: Converged') or line.startswith('SCF: Converged'):
            converged = True

        if line.startswith('IonicDynamics: Step:'):
        # if line.startswith('IonicMinimize: Iter:'):
            tokens = line.split()
            iStep = int(tokens[2])
            if converged or not convergence_check:
                stepActive = (iStep % nEvery == 0)
                converged = False
            
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

        #External Forces from Applied Potential (comes first):
        if forcesActive and iLine<refLine+atTotalNumb and line.startswith("forceExtIonic"):
            iRow = iLine-refLine
            tokens = line.split()
            extForces[iRow] = [ float(tok) for tok in tokens[2:5] ]
            # if iRow+1==atTotalNumb:
            #     forcesActive = False
            #     if coordsType == "Cartesian":
            #         extForces *= 1./(eV/Angstrom)
            #     else:
            #         extForces = np.dot(extForces, np.linalg.inv(R)/eV) #convert to Cartesian (eV/Angstrom)

        #Forces (comes after external forces):
        if forcesActive and iLine<refLine+atTotalNumb and line.startswith("force "):
            iRow = iLine-refLine
            tokens = line.split()
            forces[iRow] = [ float(tok) for tok in tokens[2:5] ]
            # when done
            if iRow+1==atTotalNumb:
                forcesActive = False
                # for testing
                # print(extForces[0])
                # print(forces[0])
                                
                # subtract out external forces from this step
                forces = forces - extForces
                
                if coordsType == "Cartesian":
                    forces *= 1./(eV/Angstrom)
                else:
                    forces = np.dot(forces, np.linalg.inv(R)/eV) #convert to Cartesian (eV/Angstrom)

        # if stepActive and line.startswith('# Forces in '):
        # need to read this in even if step not active to gather external force even before step is decided or not
        if line.startswith('# Forces in '):
            forcesActive = True
            refLine = iLine+1
            coordsType = line.split()[3]
        
        # if stepActive and line.startswith('Setting wave functions'):
        #     forcesActive = True
        #     refLine = iLine+2
        
        
        #Energy components:
        if stepActive and line.startswith('     Etot ='):
            #was# Energy components:

            #Not actually reading energy components at the moment (just calc/reporting PE)
            #Frame complete: write to OUTCAR
            # NumAt = [np.sum(np.char.count(atNames, 'Na')),
            #     np.sum(np.char.count(atNames, 'Mg')), 
            #     np.sum(np.char.count(atNames, 'Cl')) ]
            # PE = PE_tot - (E_NaPlus * NumAt[0] + E_MgPlusPlus * NumAt[1] + E_ClMinus * NumAt[2])/eV

            #Removed applied potential from forces and energy
            if appliedPotential:
                # print('total | applied')
                # print(PE_tot*eV, PE_Applied*eV)
                energy = PE_tot - PE_Applied
            else:
            # if True:
                energy = PE_tot
            # accumulate items after each step
            all_coords.append(np.copy(atpos))
            all_cells.append(np.copy(R.T))  # need to transpose jdftx convention to map to outcar style
            all_energies.append(np.copy(energy))
            all_forces.append(np.copy(forces))
            # len(all_forces) 
            # if iStep==3: break  # HACK for testing
    fp.close()

    unique_atom_names = np.unique(atom_names) # storted!! H O 
    atom_dict = dict(zip(list(unique_atom_names),np.arange(0,len(unique_atom_names))))

    # print(np.arange(0,len(unique_atom_names)))
    # atom_dict = dict([list(unique_atom_names), np.arange(0,len(unique_atom_names))])
    print(f'{atom_dict=}')
    
    # unique_atom_names_sort = np.unique(atom_names) # sorted but we want this to appear as they do in output
    # uniqueIndexes = np.unique(atom_names, return_index=True)[1]
    # unique_atom_names = [atom_names[index] for index in sorted(uniqueIndexes)]

    ions_per_type = []
    for atom_name in unique_atom_names:
        ions_per_type.append(len(np.where(atom_names == atom_name)[0]))  # H O

    atom_types = []
    # print(atom_dict[atom_names[0]])
    for idx in range(len(atom_names)) :
        if type_idx_zero :
            atom_types.append(atom_dict[atom_names[idx]])
        else :
            atom_types.append(atom_dict[atom_names[idx+1]])
    # print(f'{unique_atom_names=}  {ions_per_type=} {atom_types=}')
    print('atom map fixed')
    print('np.copy fixed version')
    return unique_atom_names, ions_per_type, atom_types, np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces)

# fname = './KFtesting/CCPBED3md0112978.jdftxout'
# begin = 0
# step = 10 
# ml = False
# type_idx_zero = True