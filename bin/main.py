import numpy as np
import math

import biotite.structure.io.pdb
import biotite.structure

from plot import Plot
from assemble import AssembleTC


FORCE_CONSTANT = 1

'''
Fold protein based on partial charges of atoms
'''


def rotation_matrix(axis, theta) -> np.ndarray[float, np.dtype]:
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


class Atom:
    def __init__(self, position: np.ndarray[float, np.dtype], charge: float, symbol: str, *args, **kwargs):
        self.position: np.ndarray = position
        self.partial_charge: float = charge
        self.symbol: str = symbol
        

class Protein:
    def __init__(self, file_name='', *args, **kwargs) -> None:
        if len(file_name):
            try:
                pdbfile = biotite.structure.io.pdb.PDBFile.read(
                    file_name)
                self.molecule: biotite.structure.AtomArray = biotite.structure.io.pdb.get_structure(
                    pdbfile, model=True)  # only kept for displaying the molecule
            except FileNotFoundError as ex:
                print(f"Could not load file: {file_name}. Loading from biotite info.")
                self.molecule: biotite.structure.AtomArray = AssembleTC.assemble()
        else:
            self.molecule: biotite.structure.AtomArray = AssembleTC.assemble()

        self.molecule.bonds = biotite.structure.connect_via_residue_names(  # type: ignore
            self.molecule)
        self.charges = biotite.structure.partial_charges(self.molecule)  # type: ignore

        self.aa_start_index: list[int] = [
            j for i, j, _ in self.molecule.bonds.as_array()
            if self.molecule.res_id[i] != self.molecule.res_id[j]
        ]

        self.atoms: np.ndarray[Atom, np.dtype] = np.array([Atom(pos, q, s) for pos, q, s in zip(
            self.molecule.coord, self.charges, self.molecule.element)])

        self.plot: Plot = None

    def calculate_torque_on_segs(self, split: int) -> tuple[np.ndarray[float, np.dtype], np.ndarray[float, np.dtype]]:
        pivot_point = self.atoms[split].position

        net_torque_left: np.ndarray[float, np.dtype] = np.array([0., 0., 0.])
        net_torque_right: np.ndarray[float, np.dtype] = np.array([0., 0., 0.])

        axis_left = self.atoms[0].position - pivot_point
        axis_right = self.atoms[-1].position - pivot_point
        for atom_left in self.atoms[:split]:
            for atom_right in self.atoms[split:]:
                delta_posistion = atom_left.position - atom_right.position
                distance = np.linalg.norm(delta_posistion)
                if distance > 0.1:  # TODO: idk what to put here unit wise (seems good tho)
                    direction = delta_posistion / distance
                    radius_left = axis_left * \
                        (np.linalg.norm(atom_left.position - pivot_point) / np.linalg.norm(axis_left))
                    radius_right = axis_right * \
                        (np.linalg.norm(atom_right.position - pivot_point) / np.linalg.norm(axis_right))

                    force = FORCE_CONSTANT*(atom_left.partial_charge*atom_right.partial_charge)*(distance**-2)

                    net_torque_left += np.cross((direction * force), radius_left)
                    net_torque_right += np.cross((-direction * force), radius_right)

        return net_torque_left, net_torque_right

    def rotate_segments(self, index: int) -> None:
        origin = self.atoms[index].position
        
        axis_left = self.atoms[0].position - origin
        axis_right = self.atoms[-1].position - origin

        torque_left, torque_right = \
            self.calculate_torque_on_segs(index)

        rotation_axis_left = np.cross(torque_left, axis_left)
        rotation_axis_right = np.cross(torque_right, axis_right)

        for atom in self.atoms[:index]:
            atom.position = np.dot(rotation_matrix(
                rotation_axis_left, np.linalg.norm(torque_left)
            ), atom.position - origin) + origin

        for atom in self.atoms[index:]:
            atom.position = np.dot(rotation_matrix(
                rotation_axis_right, np.linalg.norm(torque_right)
            ), atom.position - origin) + origin

    def fold(self, iterations, gui=False) -> None:
        if gui: self.plot = Plot(self.molecule, self.charges)

        # for i in range(iterations):
        #     for index in self.aa_start_index:
        #         self.rotate_segments(index)
        #         print(f'{i}: {index}')
        #         if gui: self.update_plot()

        index = self.aa_start_index[10]
        origin = self.atoms[index].position

        for atom in self.atoms[:index]:
            atom.position = np.dot(rotation_matrix(
                np.array([1., 0., 0.]), np.pi * 0.25
            ), atom.position - origin) + origin

        for atom in self.atoms[index:]:
            atom.position = np.dot(rotation_matrix(
                np.array([1., 0., 0.]), -np.pi * 0.25
            ), atom.position - origin) + origin

        self.create_plot()

    def update_plot(self) -> None:
        self.molecule._coord = np.array([atom.position for atom in self.atoms])
        self.plot.update(self.molecule)

    def create_plot(self) -> None:
        self.molecule._coord = np.array([atom.position for atom in self.atoms])
        Plot.plot(self.molecule, self.charges)


def main() -> None:

    trp_cage: Protein = Protein()
    trp_cage.fold(1)

    return


if __name__ == '__main__':
    main()
