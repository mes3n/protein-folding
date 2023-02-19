import json

import matplotlib.pyplot as plt
import numpy as np

import biotite.structure.io.pdb
import biotite.structure

class Compare:
    def __init__(self):
        pdbfile = biotite.structure.io.pdb.PDBFile.read('molecules/1l2y.pdb')
        molecule = biotite.structure.io.pdb.get_structure(pdbfile, model=1)

        self.atom_positions: np.ndarray[np.ndarray] = np.array([pos for pos in molecule.coord])  
        
        # print(molecule.res_id, molecule.coord)
        # print(molecule.res_id[0], molecule.coord[0])

        prev_res = molecule.res_id[0]
        self.positions: np.ndarray[np.ndarray] = np.array([molecule.coord[0]])
        for pos, res in zip(molecule.coord, molecule.res_id):
            if res != prev_res:
                self.positions = np.append(self.positions, np.array([pos]), axis=0)
                prev_res = res
        self.positions = np.delete(self.positions, 0, axis=0)

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

class LinePlt:
    @staticmethod
    def set_axes_equal(ax):
        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        plot_radius = 0.5*max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    @staticmethod
    def plot(*atomss, coor=False, save=''):
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.family'] = 'STIXGeneral'

        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        for atoms, c in zip(atomss, ['g', 'b', 'r', 'y', 'p']):
            if not coor:
                x, y, z = \
                    [a.position[0] for a in atoms], \
                    [a.position[1] for a in atoms], \
                    [a.position[2] for a in atoms]
            else:
                x, y, z = \
                    [c[0] for c in atoms], \
                    [c[1] for c in atoms], \
                    [c[2] for c in atoms]

            ax.scatter(x, y, z, c=c, s=100)
            ax.plot(x, y, z, color=c)

        LinePlt.set_axes_equal(ax)
        
        ax.view_init(90, 0, 0)
        # ax.set_axis_off()

        if save:
            plt.savefig(save, bbox_inches='tight')
        else:
            plt.show()
        
        plt.close(fig)


with open('results/raw.json', 'r') as f:
  # data = json.load(f)
  positions = json.load(f)['normal i=500']['positions']

comp = Compare()

# for i, position in enumerate(positions):
#   LinePlt.plot(position, coor=True, save=f'giffy/pos{i}.png')

LinePlt.plot(Compare.umeyama(np.array(positions[-1]), comp.positions), comp.positions, coor=True)

# data['normal i=500']['positions'] = data['normal i=500']['positions'][7:]

# with open('results/raw.json', 'w') as f:
#   json.dump(data, f, indent=2)