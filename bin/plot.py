import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import biotite.structure.graphics as graphics
import biotite.structure as struc
import biotite.structure.info as info


class Plot:
    # The number of iterations for the PEOE algorithm
    ITERATION_NUMBER = 6
    # The size of the element lables
    ELEMENT_FONT_SIZE = 10
    # The scaling factor of the atom 'balls'
    BALL_SCALE = 1.4
    # The higher this number, the more detailed are the rays
    N_RAY_STEPS = 20
    # The scaling factor of the 'ray' of charged molecules
    RAY_SCALE = 100
    # The transparency value for each 'ray ring'
    RAY_ALPHA = 0.01
    # The color map to use to depict the charge
    CMAP_NAME = "bwr_r"

    def __init__(self, molecule, charges):
        plt.ion()

        pca: PCA = PCA(n_components=3)
        pca.fit(molecule.coord)
        molecule = struc.align_vectors(
            molecule, pca.components_[-1], [0, 0, 1])

        # Balls should be colored by partial charge
        # Later this variable stores values between 0 and 1 for use in color map
        normalized_charges = charges.copy()
        # Show no partial charge for atoms
        # that are not parametrized for the PEOE algorithm
        normalized_charges[np.isnan(normalized_charges)] = 0
        # Norm charge values to highest absolute value
        max_charge = np.max(np.abs(normalized_charges))
        normalized_charges /= max_charge
        # Transform range (-1, 1) to range (0, 1)
        normalized_charges = (normalized_charges + 1) / 2
        # Calculate colors
        self.color_map = plt.get_cmap(Plot.CMAP_NAME)
        self.colors = self.color_map(normalized_charges)

        # Ball size should be proportional to VdW radius of the respective atom
        self.ball_sizes = np.array(
            [info.vdw_radius_single(e) for e in molecule.element]
        ) * Plot.BALL_SCALE

        # Gradient of ray strength
        # The ray size is proportional to the absolute charge value
        ray_full_sizes = self.ball_sizes + np.abs(charges) * Plot.RAY_SCALE
        self.ray_sizes = np.array([
            np.linspace(ray_full_sizes[i], self.ball_sizes[i],
                        Plot.N_RAY_STEPS, endpoint=False)
            for i in range(molecule.array_length())
        ]).T

        # The plotting begins here
        self.fig = plt.figure(figsize=(8.0, 6.0))
        self.ax = self.fig.add_subplot(111, projection="3d")

        # Plot the atoms
        # As 'axes.scatter()' uses sizes in points**2,
        # the VdW-radii as also squared
        graphics.plot_ball_and_stick_model(
            self.ax, molecule, self.colors, ball_size=self.ball_sizes**2, line_width=0.2,
            line_color=self.color_map(0.5), background_color=(.05, .05, .05), zoom=1.5
        )

        # Plots the rays
        for i in range(Plot.N_RAY_STEPS):
            self.ax.scatter(
                *molecule.coord.T, s=self.ray_sizes[i]**2, c=self.colors,
                linewidth=0, alpha=Plot.RAY_ALPHA
            )

        # Plot the colorbar
        color_bar = self.fig.colorbar(ScalarMappable(
            norm=Normalize(vmin=-max_charge, vmax=max_charge),
            cmap=self.color_map
        ))
        color_bar.set_label("Partial charge (e)", color="white")
        color_bar.ax.yaxis.set_tick_params(color="white")
        color_bar.outline.set_edgecolor("white")
        for label in color_bar.ax.get_yticklabels():
            label.set_color("white")

        self.fig.tight_layout()
        plt.show()
        plt.pause(0.04)

    def update(self, molecule):
        pca: PCA = PCA(n_components=3)
        pca.fit(molecule.coord)
        molecule = struc.align_vectors(
            molecule, pca.components_[-1], [0, 0, 1])

        self.ax = self.fig.add_subplot(111, projection="3d")

        graphics.plot_ball_and_stick_model(
            self.ax, molecule, self.colors, ball_size=self.ball_sizes**2, line_width=0.2,
            line_color=self.color_map(0.5), background_color=(.05, .05, .05), zoom=1.5
        )
        for i in range(Plot.N_RAY_STEPS):
            self.ax.scatter(
                *molecule.coord.T, s=self.ray_sizes[i]**2, c=self.colors,
                linewidth=0, alpha=Plot.RAY_ALPHA
            )


        plt.show()
        plt.pause(0.04)

    def pause(self):
        plt.ioff()

    @staticmethod
    def plot(molecule, charges, show_symbols=False):

        pca: PCA = PCA(n_components=3)
        pca.fit(molecule.coord)
        molecule = struc.align_vectors(
            molecule, pca.components_[-1], [0, 0, 1])

        # Balls should be colored by partial charge
        # Later this variable stores values between 0 and 1 for use in color map
        normalized_charges = charges.copy()
        # Show no partial charge for atoms
        # that are not parametrized for the PEOE algorithm
        normalized_charges[np.isnan(normalized_charges)] = 0
        # Norm charge values to highest absolute value
        max_charge = np.max(np.abs(normalized_charges))
        normalized_charges /= max_charge
        # Transform range (-1, 1) to range (0, 1)
        normalized_charges = (normalized_charges + 1) / 2
        # Calculate colors
        color_map = plt.get_cmap(Plot.CMAP_NAME)
        colors = color_map(normalized_charges)

        # Ball size should be proportional to VdW radius of the respective atom
        ball_sizes = np.array(
            [info.vdw_radius_single(e) for e in molecule.element]
        ) * Plot.BALL_SCALE

        # Gradient of ray strength
        # The ray size is proportional to the absolute charge value
        ray_full_sizes = ball_sizes + np.abs(charges) * Plot.RAY_SCALE
        ray_sizes = np.array([
            np.linspace(ray_full_sizes[i], ball_sizes[i],
                        Plot.N_RAY_STEPS, endpoint=False)
            for i in range(molecule.array_length())
        ]).T

        # The plotting begins here
        fig = plt.figure(figsize=(8.0, 6.0))
        ax = fig.add_subplot(111, projection="3d")

        # Plot the atoms
        # As 'axes.scatter()' uses sizes in points**2,
        # the VdW-radii as also squared
        graphics.plot_ball_and_stick_model(
            ax, molecule, colors, ball_size=ball_sizes**2, line_width=0.2,
            line_color=color_map(0.5), background_color=(.05, .05, .05), zoom=1.5
        )

        # Plot the element labels
        if show_symbols:
            for atom in molecule:
                ax.text(
                    *atom.coord, atom.element,
                    fontsize=Plot.ELEMENT_FONT_SIZE, color="black",
                    ha="center", va="center", zorder=100
                )

        # Plots the rays
        for i in range(Plot.N_RAY_STEPS):
            ax.scatter(
                *molecule.coord.T, s=ray_sizes[i]**2, c=colors,
                linewidth=0, alpha=Plot.RAY_ALPHA
            )

        # Plot the colorbar
        color_bar = fig.colorbar(ScalarMappable(
            norm=Normalize(vmin=-max_charge, vmax=max_charge),
            cmap=color_map
        ))
        color_bar.set_label("Partial charge (e)", color="white")
        color_bar.ax.yaxis.set_tick_params(color="white")
        color_bar.outline.set_edgecolor("white")
        for label in color_bar.ax.get_yticklabels():
            label.set_color("white")

        fig.tight_layout()
        plt.show()
