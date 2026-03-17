"""
The Kramers Kronig module.

This base modules contains bindings to essential classes and functions for the calculation
of Kramer-Kronig transforms. In particular, the module provides the following classes:
- `stoichiometry`: A class for the calculation of stoichiometry in a chemical compound.
- `kk_transforms`: A set of functions for the calculation of Kramers-Kronig transforms.
- `conversions`: A set of functions to convert between different data types
(e.g. atomic scattering factors to absorption/dispersion coefficients).
- `KK_Datatype`: An enumeration class for the data types used in `factors` object.
- `factors`: A set of classes to wrap and add methods to experimental data.
- `polynomials`: A set of classes for the calculation of the Kramer-Kronig transforms.
"""

import importlib.metadata
from kkcalc2.stoich import stoichiometry
from kkcalc2 import transforms
from kkcalc2.models import (
    conversions,
    polynomials,
    factors,
    KK_Datatype,
    PROPERTIES_DICT,
    PROPERTIES_DICT_NO_STOICH,
)
from kkcalc2 import models


def get_installed_packages() -> tuple[list[str], list[str]]:
    """
    List all installed packages with their versions.

    Returns
    -------
    list[tuple[str, str]]
        A list of tuples containing package names and their versions.
    """
    distributions = importlib.metadata.distributions()
    installed_packages = []
    versions = []
    for dist in distributions:
        installed_packages.append(dist.metadata["Name"])
        versions.append(dist.version)
    return installed_packages, versions


name = None
installed_packages, _ = get_installed_packages()
installed_packages = [pkg.lower() for pkg in installed_packages]
try:
    # Import the GUI module if appropriate packages are available:
    req = importlib.metadata.requires("kkcalc2")
    if req is not None:
        for value in req:
            value = value.replace("'", '"')
            if 'extra == "gui"' in value and ";" in value:
                name_version = value.split(";")[0]
                name = name_version
                for delim in ["~=", ">=", "==", "<=", "!=", ">", "<"]:
                    if delim in name:
                        name = name.split(delim)[0]
                name = name.strip()
                # Check that the module is available
                if name not in installed_packages:
                    raise ImportError(
                        f"kkcalc2 initialisation: Required package '{name}' is not installed."
                    )
                # Or check that the module can be imported...
                # module = __import__(name)

        # Attempted import on GUI module.
        from kkcalc2.gui import kk_gui
    else:
        print("kkcalc2 initialisation: No requirements loaded.")

except ImportError as e:
    if name is not None:
        print(
            f"kkcalc2 initialisation: GUI module import failed, requires module:\t{name}"
        )
    else:
        print("kkcalc2 initialisation: GUI module import failed.", e)

# Cleanup extra names
locs = locals()
if "value" in locs:
    del value
if "name" in locs:
    del name
if "delim" in locs:
    del delim
if "name_version" in locs:
    del name_version
if "installed_packages" in locs:
    del installed_packages

# Define the version of the package:
__version__ = importlib.metadata.version("kkcalc2")

# Traversable items
if "kk_gui" in locals():
    __all__ = [
        "stoichiometry",
        "transforms",
        "conversions",
        "polynomials",
        "KK_Datatype",
        "factors",
        "PROPERTIES_DICT",
        "PROPERTIES_DICT_NO_STOICH",
        "models",
        "__version__",
        "kk_gui",
    ]
else:
    __all__ = [
        "stoichiometry",
        "transforms",
        "conversions",
        "polynomials",
        "KK_Datatype",
        "factors",
        "PROPERTIES_DICT",
        "PROPERTIES_DICT_NO_STOICH",
        "models",
        "__version__",
    ]
