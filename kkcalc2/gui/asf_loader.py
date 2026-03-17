"""
Object loader and lister for objects that implement asf_abstract and asp_abstract classes.
Allows the loading of raw data and duplication objects.
"""

from PyQt6 import QtWidgets, QtCore
from kkcalc2.models import (
    asf_abstract,
    asf,
    asf_re,
    asf_im,
    asf_complex,
    asp_abstract,
    asp,
    asp_re,
    asp_im,
    asp_complex,
    asp_db_im_extended,
    asp_db_re_extended,
    asp_db_complex_extended,
    KK_Datatype,
)
from kkcalc2.gui.dialogs import (
    factor_complexity_dialog,
    import_data_dialog,
    factor_dtype_dialog,
)
from kkcalc2.gui.contrast_viewer import contrast_viewer
import warnings


class kk_object_list(QtWidgets.QWidget):
    """
    A widget for loading and listing objects that implement the `asf_abstract` and `asp_abstract` classes.

    Each object is represented by a row in a table,
    with columns for the object name, stoichiometry, class type, and
    a visibility checkbox.
    """

    viewSelectionChanged = QtCore.pyqtSignal()
    """A signal emitted when the 'viewing' selection changes, by toggling visibility checkboxes."""

    selectedObjectChanged = QtCore.pyqtSignal()
    """A signal emitted when the selected row (object) changes."""

    def __init__(
        self, parent=None, objs: list[type[asf_abstract | asp_abstract]] | None = None
    ):
        super().__init__(parent=parent)
        self.setWindowTitle("kkcalc Object Loader")
        self._layout = QtWidgets.QVBoxLayout()
        self.setLayout(self._layout)

        # Setup margins if parent is provided.
        if parent is not None:
            self.setContentsMargins(0, 0, 0, 0)
            self._layout.setContentsMargins(0, 0, 0, 0)

        # Create the load data buttons
        self.import_data_btn = QtWidgets.QPushButton("Import Data")

        # Create the table
        self.table = QtWidgets.QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(["Name", "Stoich.", "Type", "Vis."])
        for i in range(4):
            self.table.setColumnWidth(i, 50)
        self.table.setMinimumWidth(200)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )

        # Change the table properties
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )

        # Create duplicate, delete and contrast buttons
        hlayout = QtWidgets.QHBoxLayout()
        self.duplicate_btn = QtWidgets.QPushButton("Duplicate")
        self.duplicate_btn.setToolTip("Duplicate the selected object.")
        self.delete_btn = QtWidgets.QPushButton("Delete")
        self.delete_btn.setToolTip("Delete the selected object.")
        self.contrast_btn = QtWidgets.QPushButton("Calc. Contrast")
        self.contrast_btn.setToolTip(
            "Calculate the contrast between the two selected objects.\nSquare of differences between dispersive and absorptive components."
        )
        # Create extend button for multiple edge calculations
        self.extend_btn = QtWidgets.QPushButton("Extend across multiple edges")
        self.extend_btn.setToolTip("Extend the selected object across multiple edges.")
        self.extend_btn.setEnabled(False)

        # Create a binding for contrast viewers
        self.contrast_viewers: list[contrast_viewer] = []

        # Assign elements to the layout
        self._layout.addWidget(self.import_data_btn)
        self._layout.addWidget(self.table)
        self._layout.addLayout(hlayout)
        hlayout.addWidget(self.duplicate_btn)
        hlayout.addWidget(self.delete_btn)
        hlayout.addWidget(self.contrast_btn)
        self._layout.addWidget(self.extend_btn)

        # Setup the object sets
        self._visible_rows: set[QtWidgets.QTableWidgetItem] = (
            set()
        )  # initialize the set
        """A set of visible rows numbers, corresponding to `_objs` keys."""
        self._objs: dict[int, type[asf_abstract | asp_abstract]] = {}
        """A mapping of table row to object."""

        # Add objects to the table
        if objs is not None:
            self.add_kk_objs(objs)

        # Setup widget properties
        self.duplicate_btn.setEnabled(False)
        self.delete_btn.setEnabled(False)
        self.contrast_btn.setEnabled(False)

        # Setup connections
        self.table.itemClicked.connect(self.itemViewClicked)
        self.import_data_btn.clicked.connect(self.import_data)
        self.table.itemSelectionChanged.connect(self.on_row_selection_change)
        self.duplicate_btn.clicked.connect(self.duplicate)
        self.delete_btn.clicked.connect(self.delete)
        self.contrast_btn.clicked.connect(self.calc_contrast)
        self.extend_btn.clicked.connect(self.extend_multiple)

    def update_kk_obj(self, obj: type[asf_abstract | asp_abstract]) -> None:
        """
        Updates the row matching the object in the table.
        """
        for row, obj_ in self._objs.items():
            if obj_ == obj:
                items = [
                    obj.name,  # Name
                    str(obj.stoichiometry),  # Stoich
                    obj.__class__.__name__,
                ]  # Class type
                for i, item in enumerate(items):
                    self.table.item(row, i).setText(item)
                    if i == 2:
                        self.table.item(row, i).setToolTip(obj.__class__.__doc__)
                    else:
                        self.table.item(row, i).setToolTip(item)
                return

    def add_kk_obj(self, obj: type[asf_abstract | asp_abstract]) -> None:
        """
        Adds a new object to the table.

        Generates a 4 column table row for the object, and also stores
        the object in the internal mapping.

        Parameters
        ----------
        obj : type[asf_abstract | asp_abstract]
            The object to add to the table.
        """
        # Check if the object is an asf or asp base object; needs to be designated as real or imag.
        if obj.__class__ is asf or obj.__class__ is asp:
            # Create dialog to convert to real or imag
            from kkcalc2.gui.dialogs import factor_complexity_dialog

            dialog = factor_complexity_dialog()
            dialog.setWindowTitle(
                "Select Complexity for `"
                + obj.__class__.__name__
                + ":"
                + obj.name
                + "`"
            )
            dialog.show()
            if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                if obj.__class__ is asf:
                    obj = (
                        asf_re.from_asf(obj)
                        if dialog.complexity == dialog.EnumComplexity.REAL
                        else asf_im.from_asf(obj)
                    )
                elif obj.__class__ is asp:
                    obj = (
                        asp_re.from_asp(obj)
                        if dialog.complexity == dialog.EnumComplexity.REAL
                        else asp_im.from_asp(obj)
                    )
                else:
                    return
            else:
                return

        # Add the entry to the table
        rows = self.table.rowCount()
        self.table.setRowCount(rows + 1)
        # Check row doesn't already exist in mappings
        if rows in self._objs or rows in self._visible_rows:
            raise ValueError(f"Row {rows} already exists in the object mapping.")
        # Add the object to the mapping
        self._visible_rows.add(rows)
        self._objs[rows] = obj
        # Add the object data to the table
        obj_name = QtWidgets.QTableWidgetItem(obj.name)
        obj_name.setToolTip(obj.name)
        obj_stoich = QtWidgets.QTableWidgetItem(str(obj.stoichiometry))
        obj_stoich.setToolTip(str(obj.stoichiometry))
        obj_class = QtWidgets.QTableWidgetItem(obj.__class__.__name__)
        obj_class.setToolTip(obj.__class__.__doc__)
        self.table.setItem(rows, 0, obj_name)
        self.table.setItem(rows, 1, obj_stoich)
        self.table.setItem(rows, 2, obj_class)
        checkbox = QtWidgets.QTableWidgetItem()
        checkbox.setFlags(
            QtCore.Qt.ItemFlag.ItemIsUserCheckable | QtCore.Qt.ItemFlag.ItemIsEnabled
        )
        checkbox.setCheckState(QtCore.Qt.CheckState.Checked)
        self.table.setItem(rows, 3, checkbox)

        # Autoscale the table column widths
        # self.table.resizeColumnsToContents()

        # Emit a signal
        self.viewSelectionChanged.emit()
        return

    def add_kk_objs(self, objs: list[type[asf_abstract | asp_abstract]]) -> None:
        """
        Adds multiple objects to the table.

        Parameters
        ----------
        objs : list[type[asf_abstract  |  asp_abstract]]
            A list of objects to add to the table.

        See Also
        --------
        add_kk_obj : Adds a single object to the table.
        """
        for obj in objs:
            self.add_kk_obj(obj)
        # self.table.update()
        return

    def itemViewClicked(self, item: QtWidgets.QTableWidgetItem):
        # Check if the item belongs to the checkbox column
        if item.column() == 3:
            # Check if the item is checked
            if (
                item.checkState() == QtCore.Qt.CheckState.Checked
                and item.row() not in self._visible_rows
            ):
                self._visible_rows.add(item.row())
            elif item.row() in self._visible_rows:
                self._visible_rows.remove(item.row())
            self.viewSelectionChanged.emit()
        return

    def on_row_selection_change(self):
        # Alert the selected object has changed.
        self.selectedObjectChanged.emit()
        # Enable or disable the duplicate and delete buttons
        if self.selected_object is not None:
            self.duplicate_btn.setEnabled(True)
            self.delete_btn.setEnabled(True)
        else:
            self.duplicate_btn.setEnabled(False)
            self.delete_btn.setEnabled(False)
        # Enable or disable the contrast button if two rows are selected.
        selection = self.table.selectedItems()
        rows = set([item.row() for item in selection])
        if len(rows) == 2 and all(
            [isinstance(self._objs[row], (asf_complex, asp_complex)) for row in rows]
        ):
            self.contrast_btn.setEnabled(True)
        else:
            self.contrast_btn.setEnabled(False)
        # Enable or disable the extend button if two rows are selected.
        objs = [self._objs[row] for row in rows]
        materials = [obj.stoichiometry for obj in objs]
        if (
            len(materials) > 1  # Check if there are multiple rows selected
            and all(
                [material == materials[0] for material in materials]
            )  # Check if all materials are the same
            and all(
                [obj.__class__ is objs[0].__class__ for obj in objs]
            )  # Check if all objects are the same class
            and all([not self._objs[row].is_extended for row in rows])
        ):  # Check if all objects are not already extended:
            self.extend_btn.setEnabled(True)
            # self.extend_btn.setHidden(False)
        else:
            self.extend_btn.setEnabled(False)
            # self.extend_btn.setHidden()

    @property
    def checked_objects(self) -> list[type[asf_abstract | asp_abstract]]:
        return [self._objs[row] for row in self._visible_rows]

    @property
    def selected_object(self) -> type[asf_abstract | asp_abstract] | None:
        selected = self.table.selectedItems()
        if len(selected) == 0:
            return None
        return self._objs[selected[0].row()]

    def import_data(self):
        # Collect the raw data
        window_data = import_data_dialog()
        window_data.show()
        if window_data.exec():
            data_e, data_y = window_data.selected_data
            # Collect the complexity
            window_complexity = factor_complexity_dialog(name=window_data.load_filename)
            window_complexity.show()
            if window_complexity.exec():
                complexity = window_complexity.complexity
                if complexity == window_complexity.EnumComplexity.REAL:
                    obj = asf_re(
                        energies=data_e,
                        factors=data_y,
                        origin_dtype=KK_Datatype.ASF,
                        name=window_data.load_filename,
                    )
                    warnings.warn(
                        "Real data loaded as asf_re object, implementation required for form."
                    )
                    self.add_kk_obj(obj)
                elif complexity == window_complexity.EnumComplexity.IMAGINARY:
                    # Collect the datatype
                    window_dtype = factor_dtype_dialog(name=window_data.load_filename)
                    window_dtype.show()
                    if window_dtype.exec():
                        dtype = window_dtype.datatype
                        match dtype:
                            case KK_Datatype.ASF:
                                obj = asf_im(
                                    energies=data_e,
                                    factors=data_y,
                                    origin_dtype=KK_Datatype.ASF,
                                    name=window_data.load_filename,
                                )
                            case KK_Datatype.REFRACTIVE:
                                obj = asf_im.from_refractive(
                                    energies=data_e,
                                    refractive=data_y,
                                    name=window_data.load_filename,
                                )
                            case KK_Datatype.REFRACTIVE_INDEX:
                                obj = asf_im.from_refractive_index(
                                    energies=data_e,
                                    refractive_index=data_y,
                                    name=window_data.load_filename,
                                )
                            case KK_Datatype.NEXAFS:
                                obj = asf_im.from_NEXAFS(
                                    energies=data_e,
                                    NEXAFS=data_y,
                                    name=window_data.load_filename,
                                )
                            case _:
                                raise ValueError("Invalid datatype selected.")
                        self.add_kk_obj(obj)
        return

    def duplicate(self) -> None:
        """Duplicates the currently selected object in the table."""
        selected = self.selected_object
        if selected is None:
            return
        # Duplicate the object
        copy = selected.copy()
        # Get all names in the table
        names = [self.table.item(i, 0).text() for i in range(self.table.rowCount())]
        copy_name = copy.name + "_copy"
        copy_name_i = copy_name
        i = 0
        while copy_name_i in names:
            copy_name_i = copy_name + str(i)
            i += 1
        # Set the name of the copy
        copy.name = copy_name_i
        self.add_kk_obj(copy)

        # Emit view change signal if selected object is visible
        if self.table.selectedItems()[0].row() in self._visible_rows:
            self.viewSelectionChanged.emit()

    def delete(self) -> None:
        """Deletes the currently selected object from the table."""
        selected = self.selected_object
        if selected is None:
            return
        # Delete the row and the object from the mapping
        row = self.table.selectedItems()[0].row()
        self._objs.pop(row)
        was_visible: bool = False
        if row in self._visible_rows:
            self._visible_rows.remove(row)
            was_visible = True

        # Update the table
        self.table.removeRow(row)

        # Adjust the row numbers
        if len(self._objs) > row:
            # Shift each row backwards
            for i in range(row + 1, len(self._objs) + 1):
                obj = self._objs.pop(i)
                if (i - 1) not in self._objs:
                    self._objs[i - 1] = obj
                else:
                    raise ValueError(
                        f"Row {i - 1} already exists in the object mapping."
                    )
                if i in self._visible_rows:
                    self._visible_rows.remove(i)
                    self._visible_rows.add(i - 1)

        # Select the row prior
        if row > 0:
            self.table.selectRow(row - 1)
        elif len(self._objs) > 0:
            self.table.selectRow(0)
        else:
            self.table.clearSelection()

        # Signal a change
        if was_visible:
            self.viewSelectionChanged.emit()

    def calc_contrast(self):
        """Calculates the contrast between multiple selected objects if they are complex."""
        selection = self.table.selectedItems()
        rows = set([item.row() for item in selection])
        # Collect the complex objects
        objs = [
            self._objs[row]
            for row in rows
            if isinstance(self._objs[row], (asf_complex, asp_complex))
        ]

        # Create a contrast viewer
        viewer = contrast_viewer(objs=objs)
        viewer.show()
        self.contrast_viewers.append(viewer)  # Prevent garbage collection
        return

    def extend_multiple(self):
        """Extends the selected objects across multiple edges."""
        selection = self.table.selectedItems()
        rows = set([item.row() for item in selection])
        # Collect the rows
        objs = [self._objs[row] for row in rows]
        # Convert any polynomial objects to factors
        objs_asf: list[asf_abstract] = []
        for obj in objs:
            if isinstance(obj, asf_abstract):
                objs_asf.append(obj)
            elif isinstance(obj, asp_abstract):
                objs_asf.extend(obj.to_asf())
            else:
                raise ValueError(
                    "Invalid object selected to extend across multiple edges."
                )

        # Extend the objects
        if isinstance(objs[0], asf_im) and all(
            [isinstance(obj, asf_im) for obj in objs]
        ):
            new_obj = asp_db_im_extended(
                objs_asf, objs[0].stoichiometry, None, name=objs[0].name + "_extended"
            )
        elif isinstance(objs[0], asf_re) and all(
            [isinstance(obj, (asf_re)) for obj in objs]
        ):
            new_obj = asp_db_re_extended(
                objs_asf, objs[0].stoichiometry, None, name=objs[0].name + "_extended"
            )
        elif isinstance(objs[0], asf_complex) and all(
            [isinstance(obj, asf_complex) for obj in objs]
        ):
            new_obj = asp_db_complex_extended(
                objs_asf, objs[0].stoichiometry, None, name=objs[0].name + "_extended"
            )
        else:
            raise ValueError("Invalid object selected to extend across multiple edges.")
        self.add_kk_obj(new_obj)
        return
