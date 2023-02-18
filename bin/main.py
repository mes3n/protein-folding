import numpy as np
import math

import biotite.structure.io.pdb
import biotite.structure

from plot import Plot
from lineplt import LinePlt
from assemble import AssembleTC
from compare import Compare

import sys

from matplotlib import pyplot as plt

import json


# These are all relative and differences correspond to averge values
STEP                 = 8.0e-4
ELECTROSTATIC_WEIGHT = 1.0e-4
HYDROPHOBIC_WEIGHT   = 1.0e-1
HBOND_WEIGHT         = 1.0e-2


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
    def __init__(self) -> None:

        self.molecule: biotite.structure.AtomArray = AssembleTC.assemble()
        
        self.molecule.bonds = biotite.structure.connect_via_residue_names(  # type: ignore
            self.molecule)
        self.charges = biotite.structure.partial_charges(self.molecule)  # type: ignore


        print(f'{sum(abs(q) for q in self.charges) / len(self.charges)=}')
        # TODO: create atom without H-C H atoms
        
        self.aa_start_index: list[int] = [
            j for i, j, _ in self.molecule.bonds.as_array()
            if self.molecule.res_id[i] != self.molecule.res_id[j]
        ]  # add the latter index of the bond if the bonds are on seperate aminoacids

        self.atoms: np.ndarray[Atom, np.dtype] = np.array([Atom(pos, q, s) for pos, q, s in zip(
            self.molecule.coord, self.charges, self.molecule.element)])  
        # create custom molecule class from biotite molecule

        self.plot: Plot = None  # variable to manage background plot
        self.comp: Compare = Compare()

        self.torques = []
        self.prev_aa_positions = self.aa_positions

    @property
    def aa_positions(self):
        return np.array([self.atoms[i].position for i in self.aa_start_index])

    def calculate_torque_on_segs(self, split: int) -> tuple[np.ndarray[float, np.dtype], np.ndarray[float, np.dtype]]:
        pivot_point = self.atoms[split].position

        net_torque_left: np.ndarray[float, np.dtype] = np.array([0., 0., 0.])
        net_torque_right: np.ndarray[float, np.dtype] = np.array([0., 0., 0.])

        molecule_center = np.mean([a.position for a in self.atoms], axis=0)

        # add comment
        for atom_left in self.atoms[:split]:
            for atom_right in self.atoms[split:]:
                delta_posistion = atom_left.position - atom_right.position  # from right to left
                distance = np.linalg.norm(delta_posistion)  
                if distance > 0.0: 
		
                    # CALCULATE ELECTROSTATIC FORCE (WEAK)
                    dir_right_to_left = delta_posistion / distance

                    electrostatic_force = ELECTROSTATIC_WEIGHT * (atom_left.partial_charge * atom_right.partial_charge) * (distance**-2)
                    if atom_left.partial_charge*atom_right.partial_charge > 0.0:
                        electrostatic_force *= -1.0  # repuslive force if same signs

                    # CALCULATE ATTRICTION FROM HYDROGEN BONDS (STRONG)  # TODO
                    h_bond_force = 0.0
                    if distance < 0.2:
                        if (atom_left.symbol == 'H' and atom_right.symbol in ('O', 'N')) \
                        or (atom_right.symbol == 'H' and atom_left.symbol in ('O', 'N')):
                            h_bond_force = HBOND_WEIGHT
                            print('H_BOND')

                    net_torque_left += np.cross(
                        -dir_right_to_left * (electrostatic_force + h_bond_force)
                    , atom_right.position - pivot_point)
                    net_torque_right += np.cross(
                        dir_right_to_left * (electrostatic_force + h_bond_force)
                    , atom_left.position - pivot_point)

        for atom in self.atoms[:split]:
            # CALCULATE HYDROPHOBIC BEHAVIOR (STRONG)  # TODO
            # currently with magic number 0.2 (max charge is 0.4) to create +force on small charges and -force on large ones
            to_center  = molecule_center - atom.position
            hydrophobicity  = HYDROPHOBIC_WEIGHT * (0.1 - abs(atom.partial_charge)) * to_center

            net_torque_left += np.cross(
                hydrophobicity,
            atom.position - pivot_point)

        for atom in self.atoms[split:]:
            # CALCULATE HYDROPHOBIC BEHAVIOR (STRONG)  # TODO
            # currently with magic number 0.2 (max charge is 0.4) to create +force on small charges and -force on large ones
            to_center  = molecule_center - atom.position
            hydrophobicity  = HYDROPHOBIC_WEIGHT * (0.1 - abs(atom.partial_charge)) * to_center

            net_torque_right += np.cross(
                hydrophobicity,
            atom.position - pivot_point)

        return net_torque_left, net_torque_right

    def rotate_segments(self) -> None:

        print('calculating torque')
        torques = []
        for index in self.aa_start_index:
            print(f'at {index=}')
            torques.append(self.calculate_torque_on_segs(index))

        print('rotating segments')
        for n, index in enumerate(self.aa_start_index):
            torque_left, torque_right = torques[n]

            origin = self.atoms[index].position

            axis_left = self.atoms[0].position - origin
            axis_right = self.atoms[-1].position - origin

            abs_torque_left = np.linalg.norm(torque_left) * STEP
            abs_torque_right = np.linalg.norm(torque_right) * STEP

            # self.torques.append((torque_left, torque_right))

            rotation_axis_left = np.cross(torque_left, axis_left)
            rotation_axis_right = np.cross(torque_right, axis_right)

            rot_mat_left = rotation_matrix(rotation_axis_left, abs_torque_left)
            rot_mat_right = rotation_matrix(rotation_axis_right, abs_torque_right)

            for atom in self.atoms[:index]:
                atom.position = np.dot(rot_mat_left, atom.position - origin) + origin

            for atom in self.atoms[index:]:
                atom.position = np.dot(rot_mat_right, atom.position - origin) + origin

            # TODO: TODO: TODO
            # Make sure bonds cant be bent too much
            if 0 < index < len(self.atoms) - 1:
                v1 = self.atoms[index - 1].position - self.atoms[index].position
                v2 = self.atoms[index + 1].position - self.atoms[index].position
                theta = np.arccos(np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)))
                if theta < np.pi * 0.5:
                    print('Tight angle!', theta)
                    rotation_axis_bond = np.cross(v1, v2)
                    rot_mat_left_reverse = rotation_matrix(rotation_axis_bond, -(np.pi*0.5 - theta)*0.5)
                    rot_mat_right_reverse = rotation_matrix(rotation_axis_bond, (np.pi*0.5 - theta)*0.5)
                    for atom in self.atoms[:index]:
                        atom.position = np.dot(rot_mat_left_reverse, atom.position - origin) + origin
                    for atom in self.atoms[index:]:
                        atom.position = np.dot(rot_mat_right_reverse, atom.position - origin) + origin

            # check if atoms collide
            # if atoms collide: undo rotation
            for a1 in self.atoms[:index]:
                for a2 in self.atoms[index:]:
                    if a1 is not a2:
                        if np.linalg.norm(a1.position - a2.position) < 0.2 :
                            print('Collision!')
                            rot_mat_left_reverse = rotation_matrix(rotation_axis_left, -abs_torque_left)
                            rot_mat_right_reverse = rotation_matrix(rotation_axis_right, -abs_torque_right)
                            for atom in self.atoms[:index]:
                                atom.position = np.dot(rot_mat_left_reverse, atom.position - origin) + origin
                            for atom in self.atoms[index:]:
                                atom.position = np.dot(rot_mat_right_reverse, atom.position - origin) + origin
                        

    def fold(self, iterations, gui=0):
        step_sim = []
        similarity = []

        # if gui != 0:
        #     self.plot = Plot(self.molecule, self.charges)
        
        print(f'{self.comp.similarity(self.aa_positions)=}')
        print(f'{self.comp.similarity2(self.aa_positions, self.prev_aa_positions)=}')

        for iteration in range(iterations):
            print(f'{iteration=}')

            self.rotate_segments()

            step_sim.append(self.comp.similarity2(self.aa_positions, self.prev_aa_positions))
            similarity.append(self.comp.similarity(self.aa_positions))
            print(f'{self.comp.similarity(self.aa_positions)=}')
            print(f'{self.comp.similarity2(self.aa_positions, self.prev_aa_positions)=}')
            self.prev_aa_positions = self.aa_positions

            if gui == 1:
                # self.update_plot()
                LinePlt.plot(self.comp.umeyama(self.aa_positions, self.comp.aa_positions), self.comp.aa_positions, coor=True)  # , save=f'assets/comp{iteration}.png'
            
            if gui == 2:
                LinePlt.plot(self.comp.umeyama(self.aa_positions, self.comp.aa_positions), self.comp.aa_positions, coor=True, save=f'out/comp{iteration}.png')

            print(step_sim)
            print(similarity)

        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.family'] = 'STIXGeneral'
        # plt.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

        plt.plot(step_sim)
        plt.ylim(0, max(step_sim) * 1.2)
        plt.xlabel('Iterations')
        plt.ylabel('Distance')
        plt.xticks([])
        plt.yticks([])
        plt.title('Average distance of each atom to previous iteration')
        plt.savefig('out/step_similarity.png')
        plt.show()
        plt.close()

        plt.plot(similarity)
        plt.ylim(0, max(similarity) * 1.2)
        plt.xlabel('Iterations')
        plt.ylabel('Distance')
        plt.xticks([])
        plt.yticks([])
        plt.title('Average distance of each atom to reference protein')
        plt.savefig('out/ref_similarity.png')
        plt.show()
        plt.close()

        with open('results/raw.json', 'r') as f:
            data = json.load(f)
            data['new'] = {
                'ref_similarity': similarity,
                'step_simularity': step_sim
            }
        with open('results/raw.json', 'w') as f:
            json.dump(data, f)

        return

    def update_plot(self) -> None:
        self.molecule._coord = np.array([atom.position for atom in self.atoms])
        self.plot.update(self.molecule)

    def create_plot(self) -> None:
        LinePlt.plot(self.comp.umeyama(self.aa_positions, self.comp.aa_positions), self.comp.aa_positions, coor=True)

        self.molecule._coord = np.array([atom.position for atom in self.atoms])
        Plot.plot(self.molecule, self.charges)


def main() -> None:

    iterations = 6
    gui = 0

    if len(sys.argv) > 1:
        try:
            iterations = int(sys.argv[1])
        except ValueError:
            sys.exit("Invalid input. Enter number of iterations as an int.")

    if len(sys.argv) > 2:
        try:
            gui = int(sys.argv[2])
            assert 0 <= gui <= 2
        except (ValueError, AssertionError):
            sys.exit("Invalid input. Enter gui mode as an 0 <= int <= 2.")


    trp_cage: Protein = Protein()
    trp_cage.fold(iterations, gui)
    trp_cage.create_plot()

    return


if __name__ == '__main__':
    main()
