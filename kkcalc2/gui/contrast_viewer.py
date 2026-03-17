"""
A class that extends the asf_viewer class to display the contrast of a material.
"""

from kkcalc2.gui.asf_viewer import asf_viewer
from kkcalc2.models import asf_complex, asp_complex, KK_Datatype
from PyQt6 import QtWidgets
from typing import override
import numpy as np


class contrast_viewer(asf_viewer):
    """
    A class that extends the asf_viewer class to display the contrast of a material.

    Parameters
    ----------
    parent : QtWidgets.QWidget, optional
        The parent widget, by default None
    objs : list[type[asf_complex | asp_complex]], optional
        A list of objects to display, by default None
    """

    def __init__(
        self, parent=None, objs: list[asf_complex | asp_complex] | None = None
    ):
        super().__init__(parent=parent)
        # Force the y_datatype to be Betas.
        self.y_datatype_combo.setCurrentIndex(KK_Datatype.REFRACTIVE_INDEX.value - 1)
        self.y_datatype_combo.setEnabled(False)
        # Set the window title
        self.setWindowTitle("Contrast Viewer")
        # Set the objects
        if objs is not None:
            self.scattering_objects = objs

    @override
    def reset_graph(self) -> None:
        """
        Overrides the `asf_viewer.reset_graph` method,
        to additionally process contrast between scattering objects.
        """
        super().reset_graph()
        # Save references to existing axes
        ax1 = self.ax1
        ax2 = self.ax2
        # self.figure.clear()
        # Setup a new grid
        from matplotlib import gridspec

        grid = gridspec.GridSpec(3, 1, self.figure)
        # Reassign the existing axes to different subplot positions
        spec = ax1.get_subplotspec()
        spec1 = spec.__class__(grid, 0)
        spec2 = spec.__class__(grid, 1)
        ax1.set_subplotspec(spec1)
        ax2.set_subplotspec(spec2)
        # Get the scattering objects
        scat_objs = self.scattering_objects
        # Process contrast between any scattering objects for their contrast.
        self.ax3 = self.figure.add_subplot(313, sharex=ax2)
        self.ax1.set_xlabel("")
        self.ax2.set_xlabel("")
        self.ax3.set_xlabel("Energy (eV)")
        self.ax3.set_ylabel(r"Contrast ($\Delta\delta^2 + \Delta\beta^2$)")

        # Plot the contrasts
        for i, obj in enumerate(scat_objs[:-1]):  # Loop through all but the last object
            for j, obj2 in enumerate(
                scat_objs[i + 1 :]
            ):  # Loop through all objects after the current object
                if (
                    isinstance(obj, (asf_complex, asp_complex))
                    and isinstance(obj2, (asf_complex, asp_complex))
                    and obj.can_calc_refractive
                    and obj2.can_calc_refractive
                ):
                    # Use scattering factors
                    if isinstance(obj, asp_complex) and isinstance(obj2, asp_complex):
                        # Get all energy values
                        energies = np.unique(np.r_[obj.energies, obj2.energies])
                        factors = obj.eval_asf(energies)
                        obj = asf_complex.from_asf(
                            energies, factors, **obj._properties_dict
                        )
                        factors2 = obj2.eval_asf(energies)
                        obj2 = asf_complex.from_asf(
                            energies, factors2, **obj2._properties_dict
                        )
                    elif isinstance(obj, asp_complex):
                        energies = obj2.energies
                        factors = obj.eval_asf(energies)
                        obj = asf_complex.from_asf(
                            energies, factors, **obj._properties_dict
                        )
                    elif isinstance(obj2, asp_complex):
                        energies = obj.energies
                        factors = obj2.eval_asf(energies)
                        obj2 = asf_complex.from_asf(
                            energies, factors, **obj2._properties_dict
                        )

                    # Calculate the contrast
                    contrast_e, contrast_I = obj.contrast(obj2)
                    # Plot the contrast
                    self.ax3.plot(
                        contrast_e, contrast_I, label=f"{obj.name} - {obj2.name}"
                    )

        # Adjust the x-axis based on options
        if self.x_scale_combo.currentIndex() == 0:
            self.ax3.set_xscale("linear")
        else:
            self.ax3.set_xscale("log")

        # Match all axis to be the same
        loc = self.ax3.xaxis.get_minor_locator()
        self.ax2.xaxis.set_minor_locator(loc)
        self.ax1.xaxis.set_minor_locator(loc)
        # Remove the ticks for the top two plots
        # self.ax1.set_xticklabels(["" for _ in self.ax1.get_xticklabels()])
        # self.ax2.set_xticklabels(["" for _ in self.ax2.get_xticklabels()])

        # Add legend
        self.ax3.legend()
        self.figure.tight_layout()
        self.canvas.draw()


if __name__ == "__main__":
    # Generate some example data
    import os
    from kkcalc2 import stoichiometry
    from kkcalc2.models import asf_im

    PS_NAME = "Polystyrene"
    PS_STOICHIOMETRY = "CH"
    ps_stoich = stoichiometry(PS_STOICHIOMETRY)
    fake_stoich = stoichiometry("CHS3")

    # Import Data
    data_dir = os.path.join(os.path.dirname(__file__), "../../examples/data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)
    assert data_PS.shape[1] == 2, "Data file must have two columns"

    # Create the atomic scattering factors from NEXAFS data
    asf_PS = asf_im.from_NEXAFS(
        energies=data_PS[:, 0],
        NEXAFS=data_PS[:, 1],
        name=PS_NAME,
        stoichiometry=ps_stoich,
        density=1.1,
        scale_to_database=True,
    )
    asf_PS_complex = asf_PS.calculate_complex_polynomial(max_iter=5)

    asf_fake = asf_im.from_NEXAFS(
        energies=data_PS[:, 0],
        NEXAFS=data_PS[:, 1],
        name=PS_NAME,
        stoichiometry=fake_stoich,
        density=1.3,
        scale_to_database=True,
    )
    asf_fake_complex = asf_fake.calculate_complex_polynomial(max_iter=5)

    # Create the Application
    app = QtWidgets.QApplication([])
    app.setApplicationName("kkcalc: Contrast Viewer")

    # Create the main window
    window = contrast_viewer(objs=[asf_PS_complex, asf_fake_complex])
    # Set the x-snap to the first object
    window.snap_x(1)

    # Run the application
    window.show()
    app.exec()
