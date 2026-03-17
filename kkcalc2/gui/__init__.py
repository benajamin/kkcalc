"""
The GUI module for kkcalc.

Provides a
"""

from kkcalc2.gui.kk_gui import kk_gui


def kkcalc_app():
    """
    An example application using the kk_gui module.
    """
    import PyQt6.QtWidgets as QtWidgets

    # Create the Application
    app = QtWidgets.QApplication([])
    app.setApplicationName("kkcalc: Kramers-Kronig Calculator")

    # Generate some example data
    from kkcalc2.models import asp_db_im
    from kkcalc2 import stoichiometry

    PS_NAME = "Polystyrene"
    PS_STOICHIOMETRY = "CH"
    ps_stoich = stoichiometry(PS_STOICHIOMETRY)
    db_poly = asp_db_im(ps_stoich, name=PS_NAME)

    # Create the main window
    window = kk_gui(objs=[db_poly], autohide_modifier=True)
    window.show()

    # Run the application
    app.exec()


if __name__ == "__main__":
    kkcalc_app()
