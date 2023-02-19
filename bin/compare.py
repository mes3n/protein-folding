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

    def similarity(self, positions) -> float:
        aligned = self.umeyama(positions, self.aa_positions)

        # LinePlt.plot(self.aa_positions, aligned, coor=True)

        return np.mean(np.linalg.norm(aligned - self.aa_positions, axis=1))

    @staticmethod
    def similarity2(positions1, positions2) -> float:
        aligned = Compare.umeyama(positions1, positions2)

        # LinePlt.plot(self.aa_positions, aligned, coor=True)

        return np.mean(np.linalg.norm(aligned - positions2, axis=1))

    @staticmethod
    def umeyama(P, Q):
        assert P.shape == Q.shape
        n, _ = P.shape

        centeredP = P - P.mean(axis=0)
        centeredQ = Q - Q.mean(axis=0)

        C = np.dot(np.transpose(centeredP), centeredQ) / n

        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        R = np.dot(V, W)

        t = Q.mean(axis=0) - P.mean(axis=0).dot(1.0*R)

        return np.dot(P, R) + t


if __name__ == '__main__':
    c = Compare()
