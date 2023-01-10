import biotite.structure.io.pdb
import biotite.structure

import numpy as np


class Compare:
    def __init__(self):
        pdbfile = biotite.structure.io.pdb.PDBFile.read('molecules/1l2y.pdb')
        molecule = biotite.structure.io.pdb.get_structure(pdbfile, model=1)

        self.atom_positions: np.ndarray[np.ndarray] = np.array([pos for pos in molecule.coord])  
        
        # print(molecule.res_id, molecule.coord)
        # print(molecule.res_id[0], molecule.coord[0])

        prev_res = molecule.res_id[0]
        self.aa_positions: np.ndarray[np.ndarray] = np.array([molecule.coord[0]])
        for pos, res in zip(molecule.coord, molecule.res_id):
            if res != prev_res:
                self.aa_positions = np.append(self.aa_positions, np.array([pos]), axis=0)
                prev_res = res
        self.aa_positions = np.delete(self.aa_positions, 0, axis=0)

        self.similarity()


    def similarity(self, atoms, only_aa=True) -> float:
        comp_atom_positions = np.array(a.position for a in atoms)
        if only_aa:
            pass

        return 1.0

if __name__ == '__main__':
    c = Compare()