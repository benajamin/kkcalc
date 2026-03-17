"""
GUI interface for the Kramers-Kronig Calculator.

A simple interface that allows the user to input data, assign to real
or imaginary components, and perform a Kramers-Kronig transform on it.
The interface is built using the PyQt6 library.
"""

from PyQt6 import QtWidgets, QtCore, QtGui
import os
import io
import numpy as np
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
import pkgutil

from kkcalc2.gui.asf_viewer import asf_viewer, GraphType
from kkcalc2.gui.asf_modifier import kk_object_modifier
from kkcalc2.gui.asf_loader import kk_object_list
from kkcalc2.models import (
    asf_abstract,
    asp_abstract,
    asp_im,
    asp_re,
    asf_im,
    asf_re,
)
from kkcalc2.stoich import stoichiometry


class kk_gui(QtWidgets.QWidget):
    """
    Widget for the Kramers-Kronig Calculator.

    The widget contains a viewer, a list of objects, and a modifier for the objects.

    Parameters
    ----------
    parent : QtWidgets.QWidget | None
        The parent widget.
    objs : list[asf_abstract | asp_abstract] | None
        A list of initial objects to load into the GUI.
    autohide_modifier : bool
        Whether to autohide the modifier when no object is selected.

    Attributes
    ----------
    obj_list : kk_object_list
        The list and selector of objects.
    obj_modifier : kk_object_modifier
        The modifier (i.e. of density information or stochiometry) for the objects.
    viewer : asf_viewer
        The graphical viewer for the objects.
    """

    def __init__(
        self,
        parent=None,
        objs: list[asf_abstract | asp_abstract] | None = None,
        autohide_modifier: bool = False,
    ) -> None:  # numpydoc ignore=GL08
        super().__init__(parent=parent)
        self.setWindowTitle("Kramers-Kronig Calculator")

        windowIconPath = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "docs",
            "source",
            "_static",
            "logo2.png",
        )
        if os.path.exists(windowIconPath):
            self.setWindowIcon(QtGui.QIcon(windowIconPath))
        else:
            try:
                data = pkgutil.get_data("kkcalc", "logo2.png")
                if data is not None:
                    self.setWindowIcon(QtGui.QIcon(io.BytesIO(data)))
            except FileNotFoundError:  # as e:
                pass

        self._layout = QtWidgets.QHBoxLayout()
        self.setLayout(self._layout)
        self._draggable = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self._layout.addWidget(self._draggable)

        # Setup margins if parent is provided.
        if parent is not None:
            self.setContentsMargins(0, 0, 0, 0)
            self._layout.setContentsMargins(0, 0, 0, 0)

        # Create the viewer and the UI elements
        self.obj_list = obj_list = kk_object_list(objs=objs)
        self.obj_modifier = obj_modifier = kk_object_modifier()
        self.viewer = viewer = asf_viewer()

        # Add elements to the layout.
        # self._layout.addWidget(obj_list)
        # self._layout.addWidget(object_modifier)
        # self._layout.addWidget(viewer)
        self._draggable.addWidget(obj_list)
        self._draggable.addWidget(obj_modifier)
        self._draggable.addWidget(viewer)
        # self.obj_modifier.setFixedWidth(250)
        self._draggable.setStretchFactor(0, 10)
        self._draggable.setStretchFactor(1, 1)
        self._draggable.setStretchFactor(2, 10)

        # Modify element properties
        obj_list.setMinimumWidth(240)
        self._autohide_modifier = autohide_modifier
        if autohide_modifier:
            obj_modifier.hide()
        self._draggable.setHandleWidth(10)

        # Initialize an internal var for plotting extension handles
        self._has_handle = False
        """Tracks if the viewer has a handle for extending the domain."""
        self._handle: SpanSelector | None = None
        """The handle for extending the domain."""

        # Connect the viewer to the object list
        obj_list.viewSelectionChanged.connect(self.on_view_change)
        obj_list.selectedObjectChanged.connect(self.on_object_select)
        obj_modifier.objectModified.connect(self.on_object_modify)
        obj_modifier.objectCreated.connect(self.on_object_create)
        obj_modifier.hasHandle.connect(self.on_has_handle)

        # Setup the viewer
        self.on_view_change()

    def on_has_handle(self, has_handle: bool):
        """
        Catch the signal from the modifier, for when an extended domain is to be created.

        Parameters
        ----------
        has_handle : bool
            Whether to create the handle or not.
        """
        if has_handle:
            # Get the current obj
            obj = self.obj_list.selected_object
            if self._has_handle:
                # Update the existing handle.
                self._handle.extents = obj.energies.min(), obj.energies.max()
            else:
                # Get the current graph style
                graph_type = self.viewer.graph_type
                # Find the appropriate axes
                handle_ax = None
                if obj is None:
                    # No object selected, do nothing.
                    return
                # Only allow database extensions for re/im objects.
                if (
                    graph_type is GraphType.RE_IM_OVERLAY
                    or graph_type is GraphType.RE_IM_SEPARATE
                ):
                    if isinstance(obj, (asf_re, asp_re)):
                        handle_ax = self.viewer.ax1
                    elif isinstance(obj, (asf_im, asp_im)):
                        handle_ax = self.viewer.ax2
                # Create the handle
                if handle_ax is not None:
                    # Define functions to convert the pixel coordinates to data coordinates
                    def pix_to_data(
                        ax: plt.Axes, x: float, y: float
                    ) -> tuple[float, float]:
                        # Convert pixel coordinates to data coordinates
                        return ax.transData.transform([x, y])
                        # return inv.transform((x, y))

                    def handle_update(x, y):
                        min_x, max_x = pix_to_data(handle_ax, x, y)
                        self.on_handle_update(min_x, max_x)
                        return

                    # Create the handle
                    self._handle = SpanSelector(
                        ax=handle_ax,
                        onselect=handle_update,
                        direction="horizontal",
                        useblit=True,
                        interactive=True,
                        drag_from_anywhere=True,
                    )
                    self._handle.extents = (obj.energies.min(), obj.energies.max())
        else:
            if self._has_handle:
                # Remove the handle
                self._handle.visible = False
                self._handle.clear()
                self._handle.ax.artists.remove(self._handle)
                # Lose track of the handle
                self._handle = None
        self.viewer.reset_graph()

    def on_handle_update(self, min_x: float, max_x: float):
        """
        Catch updates to the handle, and update the modifier values.

        Parameters
        ----------
        min_x : float
            The minimum x value of the handle.
        max_x : float
            The maximum x value of the handle.
        """
        # Update the lb and ub values
        self.obj_modifier.merge_dom_lb_edit.setText(f"{min_x:.2f}")
        self.obj_modifier.merge_dom_ub_edit.setText(f"{max_x:.2f}")
        return

    def on_view_change(self):
        """
        Catch a change of toggling viewed objects.
        """
        objs = self.obj_list.checked_objects
        if self._has_handle:
            # Create a temporary asf object to pass to the viewer using the handle
            objs = objs
        # Update the viewer
        self.viewer.scattering_objects = objs

    def on_object_select(self):
        """
        Catch an object selection change, and updates the modifier view.
        """
        selected_obj = self.obj_list.selected_object
        if selected_obj is not None:
            self.obj_modifier.show()
            self.obj_modifier.object = selected_obj
        else:
            if self._autohide_modifier:
                self.obj_modifier.hide()
            else:
                # Clear the modifier QTextEdits
                self.obj_modifier.clear()

    def on_object_modify(self):
        """
        Catch an object update, and updates the table view of the object.
        """
        self.obj_list.update_kk_obj(self.obj_modifier.object)

    def on_object_create(self, new_obj: type[asf_abstract | asp_abstract]):
        """
        Catch signals that generate new objects, and adds them to the object list.

        Parameters
        ----------
        new_obj : type[asf_abstract | asp_abstract]
            The new object to add.
        """
        if new_obj is not None:
            self.obj_list.add_kk_obj(new_obj)


def demo_app():
    """
    A demo application for the kk_gui widget.
    """
    # Create the Application
    app = QtWidgets.QApplication([])
    app.setApplicationName("kkcalc: Kramers-Kronig Calculator")

    PS_NAME = "Polystyrene"
    PS_STOICHIOMETRY = "CH"
    ps_stoich = stoichiometry(PS_STOICHIOMETRY)

    # Import Data
    data_dir = os.path.join(os.path.dirname(__file__), "../data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)
    assert data_PS.shape[1] == 2, "Data file must have two columns"

    # Create the atomic scattering factors from NEXAFS data
    asf_PS = asf_im.from_NEXAFS(
        energies=data_PS[:, 0],
        NEXAFS=data_PS[:, 1],
        name=PS_NAME,
        stoichiometry=ps_stoich,
        scale_to_database=True,
    )

    # db_poly = asp_db(ps_stoich, name = PS_NAME)

    P3MEEET_NAME = "P3MEEET"
    P3MEEET_STOICHIOMETRY = "C11H16O3S"
    P3MEEET_file_O = os.path.normpath(os.path.join(data_dir, "P3MEEET_Oxygen_K.csv"))
    P3MEEET_file_S = os.path.normpath(os.path.join(data_dir, "P3MEEET_Sulfur_K.csv"))
    P3MEEET_data_O = np.genfromtxt(P3MEEET_file_O, skip_header=0, delimiter=",")
    P3MEEET_data_S = np.genfromtxt(P3MEEET_file_S, skip_header=0, delimiter=",")
    PEDOT_NAME = "PEDOT-C6C8"
    PEDOT_STOICHIOMETRY = "C21H36O2S"
    PEDOT_file_O = os.path.normpath(os.path.join(data_dir, "PEDOTC6C8_Oxygen_K.csv"))
    PEDOT_file_S = os.path.normpath(os.path.join(data_dir, "PEDOTC6C8_Sulfur_K.csv"))
    PEDOT_data_O = np.genfromtxt(PEDOT_file_O, skip_header=0, delimiter=",")
    PEDOT_data_S = np.genfromtxt(PEDOT_file_S, skip_header=0, delimiter=",")
    asf_P3MEEET_O = asf_im.from_NEXAFS(
        energies=P3MEEET_data_O[:, 0],
        NEXAFS=P3MEEET_data_O[:, 1],
        name=P3MEEET_NAME + " O",
        stoichiometry=P3MEEET_STOICHIOMETRY,
        scale_to_database=True,
    )
    asf_P3MEEET_S = asf_im.from_NEXAFS(
        energies=P3MEEET_data_S[:, 0],
        NEXAFS=P3MEEET_data_S[:, 1],
        name=P3MEEET_NAME + " S",
        stoichiometry=P3MEEET_STOICHIOMETRY,
        scale_to_database=True,
    )
    asf_PEDOT_O = asf_im.from_NEXAFS(
        energies=PEDOT_data_O[:, 0],
        NEXAFS=PEDOT_data_O[:, 1],
        name=PEDOT_NAME + " O",
        stoichiometry=PEDOT_STOICHIOMETRY,
        scale_to_database=True,
    )
    asf_PEDOT_S = asf_im.from_NEXAFS(
        energies=PEDOT_data_S[:, 0],
        NEXAFS=PEDOT_data_S[:, 1],
        name=PEDOT_NAME + " S",
        stoichiometry=PEDOT_STOICHIOMETRY,
        scale_to_database=True,
    )

    # Create the main window
    window = kk_gui(
        objs=[asf_PS, asf_P3MEEET_O, asf_P3MEEET_S, asf_PEDOT_O, asf_PEDOT_S],
        autohide_modifier=False,
    )

    # Run the application
    window.show()
    app.exec()


if __name__ == "__main__":
    demo_app()
