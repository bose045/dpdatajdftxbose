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

'''

# E_ClMinus = -15.102402370381217 #reference energy of Cl- ion
# E_NaPlus = -47.435762646473520 #reference energy of Na+ ion
# E_MgPlusPlus = -61.7620540939393123 #reference energy of Mg++ ion PBED2

#Read JDFTx input
writeForSNN = True #True writes outcar for SNN, False for DeePMD
writeOutcar = True #write outcar file or not.  If False then only write RDF
processRDF = True #True will process and save the RDF




def system_info(lines, type_idx_zero = False):
    atom_names = []
    atom_numbs = None
    nelm = None
    for ii in lines:
        ii_word_list=ii.split()
        #atom_names here

        #atom_names here
        if 'ions per type' in ii :
            atom_numbs_ = [int(s) for s in ii.split()[4:]]
            if atom_numbs is None :                
                atom_numbs = atom_numbs_
            else :
                assert (atom_numbs == atom_numbs_), "in consistent numb atoms in OUTCAR"
    assert(nelm is not None), "cannot find maximum steps for each SC iteration"
    assert(atom_numbs is not None), "cannot find ion type info in OUTCAR"
    atom_names = atom_names[:len(atom_numbs)]
    atom_types = []
    for idx,ii in enumerate(atom_numbs):
        for jj in range(ii) :
            if type_idx_zero :
                atom_types.append(idx)
            else :
                atom_types.append(idx+1)
    return atom_names, atom_numbs, np.array(atom_types, dtype = int), nelm


def get_outcar_block(fp, ml = False):
    blk = []
    energy_token = ['free  energy   TOTEN', 'free  energy ML TOTEN']
    ml_index = int(ml)
    for ii in fp :
        if not ii :
            return blk
        blk.append(ii.rstrip('\n'))
        if energy_token[ml_index] in ii:
            return blk
    return blk

# we assume that the force is printed ...
def oldget_frames(fname, begin = 0, step = 1, ml = False, convergence_check=True):
    fp = open(fname)
    blk = get_outcar_block(fp)

    atom_names, atom_numbs, atom_types, nelm = system_info(blk, type_idx_zero = True)
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []    

    cc = 0
    rec_failed = []
    while len(blk) > 0 :
        if cc >= begin and (cc - begin) % step == 0 :
            coord, cell, energy, force, virial, is_converge = analyze_block(blk, ntot, nelm, ml)
            if len(coord) == 0:
                break
            if is_converge or not convergence_check: 
                all_coords.append(coord)
                all_cells.append(cell)
                all_energies.append(energy)
                all_forces.append(force)
                if virial is not None :
                    all_virials.append(virial)
            if not is_converge:
                rec_failed.append(cc+1)

        blk = get_outcar_block(fp, ml)
        cc += 1
    
    if len(rec_failed) > 0 :
        prt = "so they are not collected." if convergence_check else "but they are still collected due to the requirement for ignoring convergence checks."
        warnings.warn(f"The following structures were unconverged: {rec_failed}; "+prt)
        
    if len(all_virials) == 0 :
        all_virials = None
    else :
        all_virials = np.array(all_virials)
    fp.close()
    return atom_names, atom_numbs, atom_types, np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces), all_virials


def analyze_block(lines, ntot, nelm, ml = False):
    coord = []
    cell = []
    energy = None
    force = []
    virial = None
    is_converge = True
    sc_index = 0
    #select different searching tokens based on the ml label
    energy_token = ['free  energy   TOTEN', 'free  energy ML TOTEN']
    energy_index = [4, 5]
    virial_token = ['FORCE on cell =-STRESS in cart. coord.  units', 'ML FORCE']
    virial_index = [14, 4]
    cell_token = ['VOLUME and BASIS', 'ML FORCE']
    cell_index = [5, 12]
    ml_index = int(ml)
    for idx,ii in enumerate(lines):
        #if set ml == True, is_converged will always be True
        if ('Iteration' in ii) and (not ml):
            sc_index = int(ii.split()[3][:-1])
            if sc_index >= nelm:
                is_converge = False
        elif energy_token[ml_index] in ii:
            energy = float(ii.split()[energy_index[ml_index]])
            assert((force is not None) and len(coord) > 0 and len(cell) > 0)
            return coord, cell, energy, force, virial, is_converge
        elif cell_token[ml_index] in ii:
            for dd in range(3) :
                tmp_l = lines[idx+cell_index[ml_index]+dd]
                cell.append([float(ss) 
                             for ss in tmp_l.replace('-',' -').split()[0:3]])
        elif virial_token[ml_index] in ii:
            in_kB_index = virial_index[ml_index]
            while idx+in_kB_index < len(lines) and (not lines[idx+in_kB_index].split()[0:2] == ["in", "kB"]) :
                in_kB_index += 1
            assert(idx+in_kB_index < len(lines)),'ERROR: "in kB" is not found in OUTCAR. Unable to extract virial.'
            tmp_v = [float(ss) for ss in lines[idx+in_kB_index].split()[2:8]]
            virial = np.zeros([3,3])
            virial[0][0] = tmp_v[0]
            virial[1][1] = tmp_v[1]
            virial[2][2] = tmp_v[2]
            virial[0][1] = tmp_v[3]
            virial[1][0] = tmp_v[3]
            virial[1][2] = tmp_v[4]
            virial[2][1] = tmp_v[4]
            virial[0][2] = tmp_v[5]
            virial[2][0] = tmp_v[5]
        elif 'TOTAL-FORCE' in ii and (("ML" in ii) == ml):
            for jj in range(idx+2, idx+2+ntot) :
                tmp_l = lines[jj]
                info = [float(ss) for ss in tmp_l.split()]
                coord.append(info[:3])
                force.append(info[3:6])
    return coord, cell, energy, force, virial, is_converge
'''