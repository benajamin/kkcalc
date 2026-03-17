"""
The loading module for the atomic scattering factor database.

Loads the database from a json file created by the db_generator_script.py script.
By default, the database file is expected to be in the same directory as this module.
Can load a compressed gzip version of the database for smaller distribution size.
"""

import os
import json
import numpy as np
import numpy.typing as npt
from typing import TypedDict
import gzip
import pkgutil
import io


class ASFElement(TypedDict):
    """
    Atomic Scattering Factor data for a single element.

    Attributes
    ----------
    E : npt.NDArray
        Photon energies corresponding to the scattering factor data points.
    Im : npt.NDArray
        Imaginary part piecewise polynomial coefficients of the scattering factor.
    Re : npt.NDArray
        Real part coefficients of the scattering factor.
    name : str
        The full name of the element.
    symbol : str
        The chemical symbol of the element.
    mass : float
        The atomic mass of the element.
    """

    E: npt.NDArray
    Im: npt.NDArray
    Re: npt.NDArray
    name: str
    symbol: str
    mass: float


def load_asf_database() -> dict[int, ASFElement]:
    """
    Load the atomic scattering factor database from a json file.

    The database has been previously created by the db_generator_script.py script.

    Returns
    -------
    dict[int, ASFElement]
            A dictionary of elements atomic number keys, with each value
            consisting of a dictionary of values (see `ASFElement`).
    """
    # For package distribution, use pkgutil to load the data file instead of file paths
    json_database = None
    try:
        gzip_json_data = pkgutil.get_data("kkcalc2.asf_database", "ASF.json.gz")
        if gzip_json_data is None:
            # Try to load the uncompressed version
            json_data = pkgutil.get_data("kkcalc2.asf_database", "ASF.json")
            if json_data is None:
                raise FileNotFoundError("ASF database file not found in package.")
            else:
                json_database = json.load(io.BytesIO(json_data))
        else:
            json_database = json.load(gzip.open(io.BytesIO(gzip_json_data)))
    except Exception as e:
        print("Failed to load ASF database via `pkgutil`.", e)

    if json_database is None:
        print("Trying file path loading...")
        path_json = os.path.join(os.path.dirname(__file__), "ASF.json")
        path_gzip_json = os.path.join(os.path.dirname(__file__), "ASF.json.gz")
        # Load all information. This inclues E, Im, and Re but also name and atomic masses.
        if os.path.exists(path_gzip_json):
            with gzip.open(path_gzip_json, "rt") as f:
                json_database = json.load(f)
        elif os.path.exists(path_json):
            with open(path_json, "r") as f:
                json_database = json.load(f)
        else:
            raise FileNotFoundError("ASF database file not found.")

    # Convert lists to numpy arrays and convert dictionary keys to integers
    asf_database = {}
    for Z in json_database.keys():
        try:
            intZ = int(Z)
            # Use the same values but with integer keys
            asf_database[intZ] = json_database[Z]
            asf_database[intZ]["E"] = np.array(json_database[Z]["E"])
            asf_database[intZ]["Im"] = np.array(json_database[Z]["Im"])
            asf_database[intZ]["Re"] = np.array(json_database[Z]["Re"])
        except ValueError:
            continue

    return asf_database


if __name__ == "__main__":
    # Test loading the database
    db = load_asf_database()
    print(f"Loaded ASF database with {len(db)} elements.")
    for Z in sorted(db.keys())[:5]:
        element = db[Z]
        print(f"Element {element['name']} (Z={Z}): {len(element['E'])} data points.")
