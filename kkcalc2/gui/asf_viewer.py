"""
File for Atomic scattering viewer GUI, build in PyQt6
"""

from PyQt6 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from enum import Enum
import matplotlib.pyplot as plt
import matplotlib.colors
from typing import Literal
from kkcalc2.models import asf_abstract, asp_abstract, asf_im, asf_re, asf_complex
from kkcalc2.models.factors import KK_Datatype
import numpy.typing as npt
import numpy as np
import warnings
import pandas as pd

# No way to attach these docs to the enum, so we will just define them here
GRAPH_TYPE_DOCS = [
    "Show the real and imaginary parts of the scattering factor separately",
    "Show the real and imaginary parts of the scattering factor together",
    "Show the absolute value and phase of the scattering factor separately",
    "Show the absolute value and phase of the scattering factor together",
]


class GraphType(Enum):
    """Graphing types for the atomic scattering factor viewer"""

    RE_IM_SEPARATE = 0
    """Show the real and imaginary parts of the scattering factor separately"""
    RE_IM_OVERLAY = 1
    """Show the real and imaginary parts of the scattering factor together"""
    ABS_PHASE_SEPARATE = 2
    """Show the absolute value and phase of the scattering factor separately"""
    ABS_PHASE_OVERLAY = 3
    """Show the absolute value and phase of the scattering factor together"""


class asf_viewer(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self._layout = QtWidgets.QVBoxLayout()
        self.setLayout(self._layout)
        # Set the margins to 0 if a parent is provided
        if parent is not None:
            self.setContentsMargins(0, 0, 0, 0)
            self._layout.setContentsMargins(0, 0, 0, 0)

        # Create settings bars
        self.top_settings_bar = QtWidgets.QHBoxLayout()
        self.bottom_settings_bar = QtWidgets.QHBoxLayout()

        # Create an option to switch graphing type
        self.graph_type_label = QtWidgets.QLabel("Repr:")
        self.graph_type_combo = QtWidgets.QComboBox()
        self.graph_type_combo.addItems(
            [x.name.lower().capitalize().replace("_", " ") for x in GraphType]
        )
        # Setup hints
        [
            self.graph_type_combo.setItemData(  # If __docs__ accessible on enums in future use that.
                i, GRAPH_TYPE_DOCS[i], QtCore.Qt.ItemDataRole.ToolTipRole
            )
            for i, x in enumerate(GraphType)
        ]

        # Create an option to switch x scale
        self.x_scale_label = QtWidgets.QLabel("X Scale:")
        self.x_scale_combo = QtWidgets.QComboBox()
        self.x_scale_combo.addItems(["Linear", "Log"])
        self.x_scale_combo.setCurrentIndex(1)  # Default to log scale

        # Create an option to switch the y-datatype.
        self.y_datatype_label = QtWidgets.QLabel("Y Datatype:")
        self.y_datatype_combo = QtWidgets.QComboBox()
        i = 0
        for x in KK_Datatype:
            if x != KK_Datatype.UNDEFINED:
                self.y_datatype_combo.addItem(x.name)
                self.y_datatype_combo.setItemData(
                    i, x.__doc__, QtCore.Qt.ItemDataRole.ToolTipRole
                )
                i += 1

        # self.y_datatype_combo.addItems([x.name for x in KK_Datatype if x != KK_Datatype.UNDEFINED])
        self.y_datatype_combo.setToolTip("The datatype to display on the y-axis")
        self.y_datatype_combo.setCurrentIndex(
            KK_Datatype.ASF.value - 1
        )  # Default to ASF

        # Create a matplotlib widget to show real and imaginary parts
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.navbar = NavigationToolbar(self.canvas, self)

        # Create a combobox to snap the x-axis to a dataset.
        self.snap_x_label = QtWidgets.QLabel("Snap X:")
        self.snap_x_combo = QtWidgets.QComboBox()

        # Create a combobox to normalise the y-axis to a dataset.
        self.norm_y_label = QtWidgets.QLabel("Norm Y:")
        self.norm_y_combo = QtWidgets.QComboBox()

        # Create a button to copy the graph data to the clipboard
        self.copy_button = QtWidgets.QPushButton("Copy Data to Clipboard")

        # Allow hiding of the legend
        self.legend_button = QtWidgets.QCheckBox("Legend")
        self.legend_button.setChecked(True)

        # Assign elements to the settings bar layouts
        self.top_settings_bar.addWidget(self.navbar)
        self.top_settings_bar.addWidget(self.legend_button)
        self.top_settings_bar.addWidget(self.x_scale_label)
        self.top_settings_bar.addWidget(self.x_scale_combo)
        self.top_settings_bar.addWidget(self.y_datatype_label)
        self.top_settings_bar.addWidget(self.y_datatype_combo)
        self.top_settings_bar.addWidget(self.graph_type_label)
        self.top_settings_bar.addWidget(self.graph_type_combo)
        self.bottom_settings_bar.addWidget(self.snap_x_label)
        self.bottom_settings_bar.addWidget(self.snap_x_combo)
        self.bottom_settings_bar.addWidget(self.norm_y_label)
        self.bottom_settings_bar.addWidget(self.norm_y_combo)
        self.bottom_settings_bar.addWidget(self.copy_button)

        # Add elements to the layout.
        self._layout.addLayout(self.top_settings_bar)
        self._layout.addWidget(self.canvas)
        self._layout.addLayout(self.bottom_settings_bar)

        # Setup the plot
        self._graph_style: GraphType = GraphType.RE_IM_SEPARATE
        """The current graphing style"""

        # Setup the scattering object list
        self._scattering_objects: list[type[asp_abstract] | type[asf_abstract]] = []

        # Connect the graph type change to the switch_graph_style function
        def graph_fn(_):
            self.reset_graph()

        self.graph_type_combo.currentIndexChanged.connect(graph_fn)
        self.x_scale_combo.currentIndexChanged.connect(graph_fn)
        self.y_datatype_combo.currentIndexChanged.connect(graph_fn)
        self.norm_y_combo.currentIndexChanged.connect(graph_fn)
        self.snap_x_combo.currentIndexChanged.connect(graph_fn)
        self.legend_button.clicked.connect(self.switch_legend)
        self.copy_button.clicked.connect(self.copy_data)
        self.reset_graph()

    @property
    def scattering_objects(self) -> list[type[asp_abstract] | type[asf_abstract]]:
        """
        The current scattering objects

        Returns
        -------
        list[type[asp_abstract] | type[asf_abstract]]
            The current scattering objects
        """
        return self._scattering_objects

    @scattering_objects.setter
    def scattering_objects(
        self,
        objs: list[
            asp_abstract | asf_abstract | type[asf_abstract] | type[asp_abstract]
        ],
    ):
        """
        The current scattering objects

        Parameters
        ----------
        list[type[asp_abstract] | type[asf_abstract]]
            The current scattering objects
        """
        # Check inputs
        for obj in objs:
            if not isinstance(obj, asf_abstract) and not isinstance(obj, asp_abstract):
                raise ValueError("Invalid scattering object")
        self._scattering_objects = objs

        # Update the snap x combo box with the new objects
        prev_idx, prev_text = (
            self.snap_x_combo.currentIndex(),
            self.snap_x_combo.currentText(),
        )
        self.snap_x_combo.blockSignals(True)
        self.snap_x_combo.clear()
        self.snap_x_combo.addItem(
            ""
        )  # always add a blank option, used to reset to full x-axis.
        self.snap_x_combo.addItems(
            [f"<{type(obj).__name__}>: {obj.name}" for obj in self.scattering_objects]
        )
        # Restore index if previous text matches an item
        if prev_text == self.snap_x_combo.itemText(prev_idx):
            self.snap_x_combo.setCurrentIndex(prev_idx)
        else:
            # Check if the previous text is in the list
            snapx_labels = [
                self.snap_x_combo.itemText(i) for i in range(self.snap_x_combo.count())
            ]
            if prev_text in snapx_labels:
                self.snap_x_combo.setCurrentIndex(snapx_labels.index(prev_text))
            # Otherwise reset is back to null selection.
        self.snap_x_combo.blockSignals(False)

        # # Update the norm y combo box with the new objects
        prev_idx, prev_text = (
            self.norm_y_combo.currentIndex(),
            self.norm_y_combo.currentText(),
        )
        self.norm_y_combo.blockSignals(True)
        self.norm_y_combo.clear()
        self.norm_y_combo.addItem(
            ""
        )  # always add a blank option, used to reset to full y-axis.
        self.norm_y_combo.addItems(
            [f"<{type(obj).__name__}>: {obj.name}" for obj in self.scattering_objects]
        )
        if prev_text == self.norm_y_combo.itemText(prev_idx):
            self.norm_y_combo.setCurrentIndex(prev_idx)
        else:
            # Check if the previous text is in the list
            normy_labels = [
                self.norm_y_combo.itemText(i) for i in range(self.norm_y_combo.count())
            ]
            if prev_text in normy_labels:
                self.norm_y_combo.setCurrentIndex(normy_labels.index(prev_text))
            # Otherwise reset is back to null selection.
        self.norm_y_combo.blockSignals(False)
        # Reset the graph
        self.reset_graph()

    @property
    def graph_type(self) -> GraphType:
        """
        The current graphing type

        Return
        ------
        GraphType
            The current graphing type
        """
        return GraphType(self.graph_type_combo.currentIndex())

    @property
    def x_scale(self) -> Literal["linear", "log"]:
        """
        The current x scale

        Return
        ------
        Literal["linear", "log"]
            The current x scale
        """
        return self.x_scale_combo.currentText().lower()

    @property
    def y_datatype(self) -> KK_Datatype:
        """
        The current y datatype

        Parameters
        ----------
        KK_Datatype
            The current y datatype

        """
        return KK_Datatype(
            self.y_datatype_combo.currentIndex() + 1
        )  # excluded UNDEFINED so need to add 1.

    @property
    def y_normalisation_scales(
        self,
    ) -> tuple[
        tuple[float, float] | None,
        tuple[float, float] | None,
        tuple[float, float] | None,
    ]:
        """
        Returns the normalisation scales for the y-axes.

        If `snap_y_combo` is not set, returns `None`.
        If `snap_y_combo` is set to a complex object, returns the min/max of the absolute value.

        Return
        ------
        re_norm_scale : tuple[float, float] | None
            The normalisation scale (min, max) for the real part of the scattering factor.
        im_norm_scale : tuple[float, float] | None
            The normalisation scale (min, max) for the imaginary part of the scattering factor.
        abs_norm_scale : tuple[float, float] | None
            The normalisation scale (min, max) for the absolute value of the scattering factor.
        """
        norm_idx = self.norm_y_combo.currentIndex()
        y_datatype = self.y_datatype
        NO_NORM = (None, None, None)
        if norm_idx <= 0 or len(self.scattering_objects) == 0:
            return NO_NORM
        else:
            norm_idx -= 1  # Subtract 1 to account for the blank option
            # Get norm obj
            obj = self.scattering_objects[norm_idx]
            # Get the x-domain to snap to
            snap_idx = self.snap_x_combo.currentIndex()
            snap_dom = self.__x_snap_domain(snap_idx)  # May be None
            # Convert asps to asfs
            if isinstance(obj, asp_abstract):
                obj: type[asp_abstract]
                obj: type[asf_abstract] = (
                    obj.to_asf()
                )  # Convert to factors to process max/min
            # Setup the return values
            ret_vals: list[tuple[float, float] | None]
            # Get the appropriate y-data type
            match y_datatype:
                case KK_Datatype.ASF | KK_Datatype.ASF_DASH:
                    ydata: npt.NDArray = obj.factors
                case KK_Datatype.REFRACTIVE_INDEX:
                    if obj.can_calc_refractive:
                        ydata: npt.NDArray = obj.refractive
                    else:
                        # Cannot calculate betas.
                        return NO_NORM
                case KK_Datatype.NEXAFS:
                    ydata: npt.NDArray = obj.NEXAFS
                case _:
                    raise ValueError(f"Invalid y datatype {y_datatype}")
            # Get the valid x-domain
            if snap_dom is not None:
                snap_idx = (obj.energies >= snap_dom[0]) & (obj.energies <= snap_dom[1])
                ydata = ydata[snap_idx]
            # Get the normalisation scales
            if isinstance(obj, asf_re) or isinstance(obj, asf_im):
                mn, mx = ydata.min(), ydata.max()
                if isinstance(obj, asf_re):
                    ret_vals = [(mn, mx), None, None]
                else:
                    ret_vals = [None, (mn, mx), None]
            elif isinstance(obj, asf_complex):
                mn_re, mx_re = ydata.real.min(), ydata.real.max()
                mn_im, mx_im = ydata.imag.min(), ydata.imag.max()
                abs_factors = np.abs(ydata)
                mn_abs, mx_abs = abs_factors.min(), abs_factors.max()
                ret_vals = [(mn_re, mx_re), (mn_im, mx_im), (mn_abs, mx_abs)]
            else:
                warnings.warn(f"Invalid scattering object {obj}")
                return NO_NORM

            # Return the normalisation scales
            return tuple(ret_vals)

    def __x_snap_domain(self, index: float) -> tuple[float, float] | None:
        """
        Returns the x-domain to snap to, if `snap_x_combo` is not set to None.

        Returns
        -------
        tuple[float, float] | None
            The x-domain (min, max) to snap to.
        """
        if index <= 0 or len(self.scattering_objects) == 0:
            return None
        else:
            index -= 1  # Subtract 1 to account for the blank option
            # Get the energies and maximum/minimum of the dataset
            energies: npt.NDArray = self.scattering_objects[index].energies
            mn, mx = energies.min(), energies.max()
            return (mn, mx)

    def snap_x(self, index: int | None = None):
        """
        Snap the x-axis to the dataset at the given index.

        When `index` is `None`, the value of the `snap_x_combo` is used.

        Parameters
        ----------
        index : int
            The index of the dataset to snap to.
        """
        if index is None:
            index = self.snap_x_combo.currentIndex()
        if index == 0:
            # Restore the x-axis to auto
            self.ax1.set_xlim(auto=True)
            self.ax2.set_xlim(auto=True)
            self.ax1.set_xmargin(0.05)  # reset margin to default
            self.ax2.set_xmargin(0.05)  # reset margin to default
        else:
            # Collect the x-domain to snap to
            dom = self.__x_snap_domain(index)
            # Stop if the domain is None
            if dom is None:
                return
            # Otherwise snap the x-axis
            mn, mx = dom
            # Use a margin:
            if self.x_scale == "linear":
                mn, mx = mn - 0.05 * (mx - mn), mx + 0.05 * (mx - mn)
            else:  # log
                mn, mx = mn / 1.05, mx * 1.05
            # Set both axes
            for ax in [self.ax1, self.ax2]:
                ax: plt.Axes
                # Set the x limits
                ax.set_xlim(mn, mx)

                # Get the y-range within the x-range
                lines = [
                    line for line in ax.get_lines() if len(line.get_xdata()) > 0
                ]  # Get all lines with data
                if len(lines) > 0:
                    y_min, y_max = None, None  # initialize
                    # Loop through the lines
                    for line in lines:
                        x, y = line.get_data()
                        mask = (x >= mn) & (x <= mx)
                        if mask.any():
                            y: npt.NDArray = y[mask]
                            # Update the y-range within the x-domain
                            y_min = y.min() if y_min is None else min(y_min, y.min())
                            y_max = y.max() if y_max is None else max(y_max, y.max())
                    # Add margins to the y-min and y-max
                    y_min, y_max = (
                        y_min - 0.05 * (y_max - y_min),
                        y_max + 0.05 * (y_max - y_min),
                    )
                    ax.set_ylim(y_min, y_max)

        # # Update the canvas
        # self.canvas.draw()

    def reset_graph(self):
        """
        Used to update the graph style and x scale attributes.
        """
        # Get graphing attributes
        graph_style = self.graph_type
        x_scale = self.x_scale
        y_datatype = self.y_datatype
        # Clear the existing graph
        for ax in self.figure.axes:
            self.figure.delaxes(ax)
        self.figure.clear()
        # Choose the y labels and title
        match y_datatype:
            case KK_Datatype.ASF:
                title = "Atomic Scattering Factors (f = f1 + i*f2)"
                yl1 = r"$f_1$ (electrons)"  # f1
                yl2 = r"$f_2$ (electrons)"  # f2
            case KK_Datatype.REFRACTIVE_INDEX:
                title = r"Refractive Index ($n = 1 - \delta - i * \beta $)"
                yl1 = r"$1-\delta$ (Dispersion)"  # Delta
                yl2 = r"$\beta$ (Absorption)"  # Beta
            case KK_Datatype.REFRACTIVE:
                title = r"Refractive Components ($n = 1 - \delta - i * \beta $)"
                yl1 = r"$\delta$ (Dispersion)"  # Delta
                yl2 = r"$\beta$ (Absorption)"  # Beta
            case KK_Datatype.NEXAFS | KK_Datatype.XANES | KK_Datatype.PHOTOABSORPTION:
                title = "Absorption intensities (A.U.)"
                yl1 = "Real (A.U.)"
                yl2 = "Imaginary (A.U.)"
            case KK_Datatype.ASF_DASH:
                title = "Atomic Scattering Factors ($f = f^0 + f' + i*f{''}$)"
                yl1 = r"$f'$ (electrons)"  # f'
                yl2 = r"$f{''}$ (electrons)"  # f''
            case _:
                raise ValueError("Invalid y datatype")
        self.figure.suptitle(title)
        # Add the axes, and x labels
        match graph_style:
            case GraphType.RE_IM_SEPARATE | GraphType.ABS_PHASE_SEPARATE:
                self.ax1 = self.figure.add_subplot(211)
                self.ax2 = self.figure.add_subplot(212)
                self.ax1.set_xlabel("Energy (eV)")
                self.ax2.set_xlabel("Energy (eV)")
                self.ax1.set_xscale(x_scale)
                self.ax2.set_xscale(x_scale)
            case GraphType.RE_IM_OVERLAY | GraphType.ABS_PHASE_OVERLAY:
                self.ax1 = self.figure.add_subplot(111)
                self.ax2 = self.ax1.twinx()
                self.ax1.set_xlabel("Energy (eV)")
                self.ax1.set_xscale(x_scale)
        # Add the y labels
        match graph_style:
            case GraphType.RE_IM_SEPARATE | GraphType.RE_IM_OVERLAY:
                self.ax1.set_ylabel(yl1)
                self.ax2.set_ylabel(yl2)
            case GraphType.ABS_PHASE_SEPARATE | GraphType.ABS_PHASE_OVERLAY:
                self.ax1.set_ylabel("Magnitude (A.U.)")
                self.ax2.set_ylabel("Phase (deg)")

        # Get the normalisation scales
        re_norm_scale, im_norm_scale, abs_norm_scale = self.y_normalisation_scales

        # Setup a discrete colormap for the lines
        scat_objs = self.scattering_objects
        if len(scat_objs) < 10:
            cmap = plt.get_cmap("tab10")
        elif len(scat_objs) < 20:
            cmap = plt.get_cmap("tab20")
        else:
            cmap1 = plt.get_cmap("tab20").colors
            cmap2 = plt.get_cmap("tab20b").colors
            cmap3 = plt.get_cmap("tab20c").colors
            cmap = matplotlib.colors.ListedColormap(cmap1 + cmap2 + cmap3)
        cl = len(cmap.colors)
        c_indx = 0  # index for the colormap

        # Plot the data
        for obj in scat_objs:
            # If a polynomial, convert to asf to display.
            if isinstance(obj, asp_abstract):
                obj = obj.to_atomic_scattering_factors()

            # Plot data of an asf object.
            x = obj.energies
            x_snap = self.__x_snap_domain(self.snap_x_combo.currentIndex())
            if x_snap is not None:
                x_dom_idx = (x >= x_snap[0]) & (x <= x_snap[1])
            else:
                x_dom_idx = None

            if isinstance(obj, asf_abstract):
                # Get the appropriate y-data
                if y_datatype == KK_Datatype.ASF:
                    y = obj.factors
                elif y_datatype == KK_Datatype.ASF_DASH:
                    y = obj.factors
                    if isinstance(
                        obj, (asf_re)
                    ):  # Real part of asf_dash is f0 + f', remove f0.
                        y = obj.factors - obj.stoichiometry.relativistic_correction
                elif y_datatype == KK_Datatype.REFRACTIVE_INDEX:
                    if obj.can_calc_refractive:
                        if isinstance(obj, asf_complex):
                            # Use the complex refractive index
                            y = obj.refractive_indexes
                        elif isinstance(obj, asf_re):
                            y = (
                                1 - obj.refractive
                            )  # The refractive index component is 1-delta.
                        elif isinstance(obj, asf_im):
                            y = obj.refractive
                        else:
                            continue  # Skip this object if it cannot calculate refractive index
                    else:
                        # Skip this object if it cannot calculate beta
                        continue
                elif y_datatype == KK_Datatype.REFRACTIVE:
                    if obj.can_calc_refractive:
                        y = obj.refractive
                    else:
                        continue
                elif y_datatype == KK_Datatype.NEXAFS:
                    if isinstance(obj, (asf_im, asf_complex)):
                        y = obj.NEXAFS
                    else:
                        continue  # Skip this object if it cannot calculate NEXAFS
                else:
                    raise ValueError(f"Invalid y datatype {y_datatype}")

                # Calculate the object normalisation scales
                if isinstance(obj, (asf_re, asf_complex)):
                    if re_norm_scale is not None:
                        if x_snap is not None and x_dom_idx is not None:
                            scale_re = (re_norm_scale[1] - re_norm_scale[0]) / (
                                y[x_dom_idx].real.max() - y[x_dom_idx].real.min()
                            )
                        else:
                            scale_re = (re_norm_scale[1] - re_norm_scale[0]) / (
                                y.real.max() - y.real.min()
                            )
                        y_re = (y.real - y.real.min()) * scale_re + re_norm_scale[0]
                    else:
                        y_re = y.real
                if isinstance(obj, (asf_im, asf_complex)):
                    # asf_im object uses real values, but asf_complex uses imaginary values.
                    # Pull out imaginary values from asf_complex
                    y_im = y.imag if isinstance(obj, asf_complex) else y

                    # Now normalise use real parts.
                    if im_norm_scale is not None:
                        if x_snap is not None and x_dom_idx is not None:
                            scale_im = (im_norm_scale[1] - im_norm_scale[0]) / (
                                y_im[x_dom_idx].max() - y_im[x_dom_idx].min()
                            )
                        else:
                            scale_im = (im_norm_scale[1] - im_norm_scale[0]) / (
                                y_im.max() - y_im.min()
                            )
                        y_im = (y_im - y_im.min()) * scale_im + im_norm_scale[0]

                if isinstance(obj, asf_complex):
                    if abs_norm_scale is not None:
                        if x_snap is not None and x_dom_idx is not None:
                            scale_abs = (abs_norm_scale[1] - abs_norm_scale[0]) / (
                                abs(y[x_dom_idx]).max() - abs(y[x_dom_idx]).min()
                            )
                        else:
                            scale_abs = (abs_norm_scale[1] - abs_norm_scale[0]) / (
                                abs(y).max() - abs(y).min()
                            )
                        y_abs = (abs(y) - abs(y).min()) * scale_abs + abs_norm_scale[0]
                    else:
                        y_abs = abs(y)
                    y_phase = np.angle(y, deg=True)

                # Determine the coloring based on if graphs are separate or overlayed.
                if graph_style in [
                    GraphType.RE_IM_SEPARATE,
                    GraphType.ABS_PHASE_SEPARATE,
                ] or not isinstance(obj, asf_complex):
                    c1 = c2 = cmap(c_indx % cl)
                    c_indx += 1
                else:
                    c1 = cmap(c_indx % cl)
                    c2 = cmap((c_indx + 1) % cl)
                    c_indx += 2

                # Cannot plot re/im on amp/phase graphs unless the object is a complex asf.
                if isinstance(obj, asf_re) and graph_style in [
                    GraphType.RE_IM_SEPARATE,
                    GraphType.RE_IM_OVERLAY,
                ]:
                    self.ax1.plot(obj.energies, y_re, label=obj.name, c=c1)
                elif isinstance(obj, asf_im) and graph_style in [
                    GraphType.RE_IM_SEPARATE,
                    GraphType.RE_IM_OVERLAY,
                ]:
                    self.ax2.plot(obj.energies, y_im, label=obj.name, c=c2)
                elif isinstance(obj, asf_complex):
                    if graph_style in [
                        GraphType.ABS_PHASE_SEPARATE,
                        GraphType.ABS_PHASE_OVERLAY,
                    ]:
                        self.ax1.plot(obj.energies, y_abs, label=obj.name, c=c1)
                        self.ax2.plot(obj.energies, y_phase, label=obj.name, c=c2)
                    elif graph_style in [
                        GraphType.RE_IM_SEPARATE,
                        GraphType.RE_IM_OVERLAY,
                    ]:
                        self.ax1.plot(obj.energies, y_re, label=obj.name, c=c1)
                        self.ax2.plot(obj.energies, y_im, label=obj.name, c=c2)
                else:
                    # If the object doesn't satisfy the above conditions, skip it.
                    if graph_style in [
                        GraphType.RE_IM_SEPARATE,
                        GraphType.ABS_PHASE_SEPARATE,
                    ] or not isinstance(obj, asf_complex):
                        # Undo the color index increment
                        c_indx -= 1
                    else:
                        c_indx -= 2
            else:
                raise ValueError("Invalid scattering object")
        if graph_style in [GraphType.RE_IM_SEPARATE, GraphType.ABS_PHASE_SEPARATE]:
            # Add the legend
            if len(self.ax1.get_lines()) > 0:
                self.ax1.legend()
            if len(self.ax2.get_lines()) > 0:
                self.ax2.legend()
        else:
            # Add the combined legend
            lines = self.ax1.get_lines() + self.ax2.get_lines()
            if len(lines) > 0:
                labels = [line.get_label() for line in lines]
                self.ax1.legend(lines, labels)
        # Snap the x-axis if needed
        self.snap_x()
        # Hid the legend if needed
        self.switch_legend()
        # Draw the plot
        self.figure.tight_layout()
        self.canvas.draw()
        return

    def switch_legend(self, state: bool | None = None):
        """
        Switch the visibility of the legend.

        Parameters
        ----------
        state : bool
            The state to set the legend visibility to.
            True to show the legend, False to hide it.
        """
        # Use the button
        state = self.legend_button.isChecked() if state is None else state
        # Set the visibility of the legend
        for ax in [self.ax1, self.ax2]:
            ax: plt.Axes
            legend = ax.get_legend()
            if legend is not None:
                legend.set_visible(state)
        # Redraw the canvas
        self.canvas.draw()

    def copy_data(self) -> pd.DataFrame:
        """
        Copy the graph data to the clipboard, and also return the data as a DataFrame.

        Returns
        -------
        pd.DataFrame
            The data as a DataFrame object.
            Header rows are (1) the sample names and (2) the data x/y designation.
        """
        # Get the current graph style
        graph_style = self.graph_type
        graph_datatype = self.y_datatype

        # Get the data from the displayed lines
        dsets_ax1: list[str, npt.NDArray, npt.NDArray] = []
        dsets_ax2: list[str, npt.NDArray, npt.NDArray] = []
        for line in self.ax1.get_lines():
            line: plt.Line2D
            x, y = line.get_data()
            name = line.get_label()
            dsets_ax1.append((name, x, y))
        for line in self.ax2.get_lines():
            line: plt.Line2D
            x, y = line.get_data()
            name = line.get_label()
            dsets_ax2.append((name, x, y))

        samples = []
        columns = []
        values = []
        match graph_style:
            case GraphType.RE_IM_SEPARATE | GraphType.RE_IM_OVERLAY:
                for name, x, y in dsets_ax1:
                    samples.append(name)
                    samples.append(name)
                    columns.append("Energy")
                    columns.append(f"{graph_datatype.name} (Re)")
                    values.append(x)
                    values.append(y)
                for name, x, y in dsets_ax2:
                    samples.append(name)
                    samples.append(name)
                    columns.append("Energy")
                    columns.append(f"{graph_datatype.name} (Im)")
                    values.append(x)
                    values.append(y)

        index = pd.MultiIndex.from_tuples(
            list(zip(samples, columns)), names=["Sample", "Data"]
        )
        df = pd.DataFrame(values, index=index).transpose()
        df.to_clipboard()
        return df


if __name__ == "__main__":
    # Create the application
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    app.setApplicationName("Atomic Scattering Factor (ASF) Viewer")
    ui = asf_viewer()

    # Import kkcalc functions
    from kkcalc2.asf_database import asp_db_extended
    from kkcalc2.models.common import atomic_scattering
    from kkcalc2.stoich import stoichiometry

    # Import some example data
    import numpy as np
    import os

    PS = "C10H8"
    PS_stoich = stoichiometry(PS)

    # Data
    data_path = os.path.join(
        os.path.dirname(__file__), "../../examples/data/PS_004_-dc.txt"
    )
    data = np.genfromtxt(data_path, skip_header=4)

    # Database extensions and representations
    polystyrene_data = asf_im.from_NEXAFS(
        data[:, 0], data[:, 1], stoichiometry=PS_stoich, name="Polystyrene Data"
    )
    polystyrene_asp_im = asp_db_extended(
        data_asf=polystyrene_data, database=PS_stoich, name="PS Ext."
    )
    polystyrene_asp_im2 = asp_db_extended(
        data_asf=polystyrene_data,
        database=PS_stoich,
        merge_domain=[280, 320],
        name="PS Ext. Domain.",
    )
    polystyrene_asf_re = polystyrene_asp_im.kk_transform()
    polystyrene_asf_re2 = polystyrene_asp_im2.kk_transform()
    polystyrene_complex = polystyrene_asp_im.calculate_complex_factors(
        name="Polystyrene Complex"
    )
    PS_DENSITY = 1.0  # g/cm^3
    polystyrene_complex_density = polystyrene_asp_im.calculate_complex_factors(
        name="Polystyrene Complex with Beta", density=PS_DENSITY
    )

    # List
    objs: list[type[atomic_scattering]] = [
        polystyrene_data,
        polystyrene_asp_im,
        polystyrene_asp_im2,
        polystyrene_asf_re,
        polystyrene_asf_re2,
        polystyrene_complex,
        polystyrene_complex_density,
    ]

    # Provide the data to the viewer.
    ui.scattering_objects = objs

    ui.show()
    sys.exit(app.exec())
