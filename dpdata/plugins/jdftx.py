import dpdata.jdftx.jdftxout
import numpy as np
from dpdata.format import Format
from dpdata.utils import sort_atom_names, uniq_atom_names

# rotate the system to lammps convention
@Format.register("jdftxout")
@Format.register("jdftx/jdftxout")
class JDFTXOutFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, convergence_check=True,Ktot=1, K=1, **kwargs):
        data = {}
        ml = kwargs.get("ml", False)  #XX
        data['atom_names'], \
            data['atom_numbs'], \
            data['atom_types'], \
            data['cells'], \
            data['coords'], \
            data['energies'], \
            data['forces'], \
            = dpdata.jdftx.jdftxout.get_frames(file_name, begin=begin, step=step, ml=ml, convergence_check=convergence_check, Ktot=Ktot, K=K)
        # if tmp_virial is not None:
        #     data['virials'] = tmp_virial
        # scale virial to the unit of eV
        # if 'virials' in data:
        #     v_pref = 1 * 1e3 / 1.602176621e6
        #     for ii in range(data['cells'].shape[0]):
        #         vol = np.linalg.det(np.reshape(data['cells'][ii], [3, 3]))
        #         data['virials'][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        return data