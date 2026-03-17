"""
This module accumulates optical database data.

This file is not included in the distributed package, but is provided
for reference.

Data is taken from  different sources and packaged into `ASF.json` to be
used by and distributed with the Kramers-Kronig Calculator software package.

Workflow to accomodate:
1. Read data from .nff files and BL files.
2. Combine Henke and BL data sets.
3. Convert to useful format for internal use.
4. Write to json file for distribution.
5. Load data for use by KKCalc.
6. Combine data for different elements as selected by user.
6.a) figure out energy values (abscissa).
6.b) add coefficients/intensity values in selected proportions
7. Provide combined data in required formats.
7.a) list of tuples for plotting.
7.b) list of energy ranges, each with corresponding list of polynomial coefficients, (i.e. piecewise polynomial format) for PP KK calculation.

Items 1-4 not usually performed by users. Items 5-7 must be integrated into KKCalc program.
"""

import os
import os.path
import scipy.interpolate
import json
import numpy.typing as npt
import numpy as np
import gzip
import pkgutil
import io

BASEDIR = os.path.dirname(os.path.realpath(__file__))
classical_electron_radius = 2.81794029957951365441605230194258e-15  # meters
Plancks_constant = 4.1356673310e-15  # eV*seconds
speed_of_light = 2.99792458e8  # meters per second
Avogadro_constant = 6.02214129e23

try:
    data_resource = pkgutil.get_data(
        "kkcalc2",
        "asf_database/db_data/elements.dat",
    )
    with io.BytesIO(data_resource) as f:
        Elements_DATA = [line.decode("utf-8").strip("\r\n").split() for line in f]
except IOError:
    # Try a relative path if pkgutil fails (e.g. when running the script directly)
    Elements_DATA = [
        line.strip("\r\n").split()
        for line in open(os.path.join(BASEDIR, "db_data", "elements.dat"))
    ]
Database = dict()


#################################################################################################################
def LoadData(filename: str) -> npt.NDArray:
    """
    A loader for the Henke .nff ascii data files.

    Parameters
    ----------
    filename : str
        The path to the .nff file to load.

    Returns
    -------
    npt.NDArray
        The loaded floating data as a numpy array.
    """
    data = []
    if os.path.isfile(filename):
        for line in open(filename):
            try:
                data.append([float(f) for f in line.split()])
            except ValueError:
                pass
        data = np.array(data)
    else:
        print("Error:", filename, "is not a valid file name.")
    if len(data) == 0:
        print("Error: no data found in", filename)
    return np.array(data)


def parse_BL_file(briggs_file: str | io.BytesIO) -> dict[int, npt.NDArray]:
    """
    Parse a Biggs and Lighthill (BL) file.

    Parameters
    ----------
    briggs_file : str
        The path to the BL file to parse.

    Returns
    -------
    dict[int, npt.NDArray]
        A dictionary containing the parsed data.
        The keys are element atomic numbers, and the values are numpy arrays of coefficients.

    Raises
    ------
    FileNotFoundError
        If the specified BL file does not exist.
    """
    continue_norm = True  # Normalise the Biggs and Lighthill data as the published scattering factors do, rather than as Henke et al says.
    BLfile = {}
    # Read the file
    if isinstance(briggs_file, str):
        if os.path.exists(briggs_file) is False:
            raise FileNotFoundError(
                f"Biggs and Lighthill file not found: {briggs_file}"
            )
        with open(briggs_file, "r") as f:
            lines = f.readlines()
    elif isinstance(briggs_file, io.BytesIO):
        briggs_file.seek(0)  # Ensure we're at the start of the BytesIO stream
        lines = briggs_file.readlines()
        lines = [line.decode("utf-8") for line in lines]
    else:
        raise TypeError("briggs_file must be a string path or a BytesIO object.")
    # Process the lines
    for line in lines:
        try:
            values = [float(f) for f in line.split()]
            if values[3] > 10:
                Norm_value = 0  # will calculate actual normalisation value later
                if (
                    not continue_norm
                    and values[2] > 10
                    and values[2] not in [20, 100, 500, 100000]
                ):
                    Norm_value = 1
                elif (
                    not continue_norm
                    and values[0] == 42
                    and values[2] > 10
                    and values[2] not in [100, 500, 100000]
                ):  # Mo needs special handling
                    # print "Mo seen at", values[0], values[2]
                    Norm_value = 1
                values.append(Norm_value)
                if values[2] not in [0.01, 0.1, 0.8, 4, 20, 100, 500, 100000] or (
                    values[0] == 42 and values[2] == 20
                ):
                    values.append(1)  # this is an absorption edge!
                else:
                    values.append(0)  # this is not an absorption edge
                BLfile[int(values[0])].append(values)
        except ValueError:
            pass
        except IndexError:
            pass
        except KeyError:
            BLfile[int(values[0])] = [values]
    for elem, coeffs in list(BLfile.items()):
        BLfile[elem] = np.array(coeffs)[:, 2:]
    return BLfile


def BL_to_ASF(E: npt.ArrayLike, coeffs: npt.NDArray, Atomic_mass: float) -> npt.NDArray:
    r"""
    The conversion factor from Biggs and Lighthill coefficients to Henke scattering factors.

    Biggs and Lighthill offers photoelectric cross-section (PECS) with the sum of

    ..math::
        PECS = \sum_{n=1}^{4} A_n * E^{-n}

    where n is the reciprocal order (n=1-4), with E in keV and PECS in cm^2/g.
    The Henke scattering factors are related by

    ..math::
        f2 = PECS*E/(2*r0*h*c),

    where E is the energy (eV), PECS cm^2/atom.

    Parameters
    ----------
    E : npt.ArrayLike
        The energies to calculate the Henke scattering factors at.
        Can be singular or an array of energies.
    coeffs : npt.NDArray
        The polynomial coefficients corresponding to the energies.
        Should have a shape of at least (4,), but another dimension
        can match the shape of E for vectorised calculations.
    Atomic_mass : float
        The atomic mass of the element being calculated (in atomic mass units).

    Returns
    -------
    npt.NDArray
        The f2 Henke scattering factors.
    """
    # If E is not a singular value, convert to array for vectorised calculation
    if not isinstance(E, (int, float)):
        E = np.asarray(E)
    return (
        (
            coeffs[0]
            + coeffs[1] / (E * 0.001)
            + coeffs[2] / ((E * 0.001) ** 2)
            + coeffs[3] / ((E * 0.001) ** 3)
        )
        * Atomic_mass
        / (
            2
            * Avogadro_constant
            * classical_electron_radius
            * Plancks_constant
            * speed_of_light
        )
        * 0.1
    )


def Coeffs_to_ASF(E: npt.ArrayLike, coeffs: npt.NDArray) -> npt.NDArray:
    r"""
    Calculate Henke scattering factors from polynomial coefficients.

    Uses the linear n=1, n=0 coefficients from Henke data, and the
    n = -1, -2, -3 coefficients from Biggs and Lighthill data.

    ..math::
        f2 = \sum_{n=0}^{4} B_n * E^{1-n}

    E in eV and PECS in cm^2/atom

    Parameters
    ----------
    E : npt.ArrayLike
        The energies to calculate the Henke scattering factors at.
        Can be singular or an array of energies.
    coeffs : npt.NDArray
        The polynomial coefficients corresponding to the energies.
        Should have a shape of at least (5,), with coefficients
        ordered from n=1 to n=-3, but another dimension can match
        the shape of E for vectorised calculations.

    Returns
    -------
    npt.NDArray
        The f2 Henke scattering factors.
    """
    if not isinstance(E, (int, float)):
        E = np.asarray(E)
    return (
        coeffs[0] * E
        + coeffs[1]
        + coeffs[2] / E
        + coeffs[3] / (E**2)
        + coeffs[4] / (E**3)
    )


###########################################################################################################

data = pkgutil.get_data(
    "kkcalc2",
    "asf_database/db_data/original_biggs_file.dat",
)
if data is not None:
    file = io.BytesIO(data)
else:
    # Try a relative path if pkgutil fails (e.g. when running the script directly)
    file = os.path.join(BASEDIR, "db_data", "original_biggs_file.dat")
    if not os.path.exists(file):
        raise FileNotFoundError("Biggs and Lighthill file not found in package data.")

BL_data = parse_BL_file(briggs_file=file)

# for z, symbol, name, atomic_mass, Henke_file in [Elements_DATA[0]]:
for z, symbol, name, atomic_mass, Henke_file in Elements_DATA:
    # print(z, symbol, name, atomic_mass, Henke_file)
    # Get basic metadata
    Element_Database = dict()
    Element_Database["mass"] = float(atomic_mass)
    Element_Database["name"] = name
    Element_Database["symbol"] = symbol

    # Get basic data
    # print("Load nff data from:", os.path.join(BASEDIR, 'data', Henke_file))
    data_path = os.path.join(BASEDIR, "db_data", Henke_file)
    asf_RawData = LoadData(data_path)
    if min(asf_RawData[1:-1, 0] - asf_RawData[0:-2, 0]) < 0:
        print(
            "Warning! Energies in ",
            Henke_file,
            "are not in ascending order! (Sorting now..)",
        )
        asf_RawData.sort()
    # print BL_data[int(z)]

    # Convert and normalise BL data
    # get normalisation values
    ASF_norm = scipy.interpolate.splev(
        10000,
        scipy.interpolate.splrep(asf_RawData[:, 0], asf_RawData[:, 2], k=1),
        der=0,
    )
    BL_norm = BL_to_ASF(10000, BL_data[int(z)][0][3:7], float(atomic_mass))
    # print "Norms:", ASF_norm, BL_norm, BL_norm/ASF_norm

    temp_E = []
    BL_coefficients = []
    for line in BL_data[int(z)]:
        if float(line[1]) >= 30:
            temp_E.append(float(line[0]))
            BL_coefficients.append(
                line[2:7]
                / BL_norm
                * ASF_norm
                * [0, 1, 1000, 1000000, 1000000000]
                * float(atomic_mass)
                / (
                    2
                    * Avogadro_constant
                    * classical_electron_radius
                    * Plancks_constant
                    * speed_of_light
                )
                * 0.1
            )
    # store for use in calculation
    C = np.array(BL_coefficients)
    # (insert 30000.1 here to use linear section from 30000 to 30000.2 to ensure continuity between data sets)
    X = np.array([30.0001] + temp_E[1:]) * 1000

    # Express asf data in PP
    M = (asf_RawData[1:, 2] - asf_RawData[0:-1, 2]) / (
        asf_RawData[1:, 0] - asf_RawData[0:-1, 0]
    )
    B = asf_RawData[0:-1, 2] - M * asf_RawData[0:-1, 0]
    E = asf_RawData[:, 0]
    # asf_RawData (i.e. E, Re, Im) matches dimensions at this stage. Energies only go to 30000 eV, so we need to extend them. Briggs Lighthill adds 4 points.

    Full_coeffs = np.zeros((len(asf_RawData[:, 0]) - 1, 5))
    Full_coeffs[:, 0] = M
    Full_coeffs[:, 1] = B
    # Append B&L data and make sure it is continuous
    E = E[0:-1]
    for i in range(len(X) - 1):
        Y1 = Coeffs_to_ASF(X[i] - 0.1, Full_coeffs[-1, :])
        Y2 = Coeffs_to_ASF(X[i] + 0.1, C[i, :])
        M = (Y2 - Y1) / 0.2
        B = Y1 - M * (X[i] - 0.1)
        E = np.append(E, [X[i] - 0.1, X[i] + 0.1])
        Full_coeffs = np.append(Full_coeffs, [[M, B, 0, 0, 0]], axis=0)
        Full_coeffs = np.append(Full_coeffs, [C[i, :]], axis=0)
    E = np.append(E, X[-1])

    # Store -9999. Re values as np.nan instead
    asf_RawData[asf_RawData[:, 1] == -9999.0, 1] = np.nan

    # convert np arrays to nested lists to enable json serialisation with the default converter.
    Element_Database["E"] = E.tolist()
    Element_Database["Im"] = Full_coeffs.tolist()
    Element_Database["Re"] = asf_RawData[:, 1].tolist()
    Database[int(z)] = Element_Database

output_path = os.path.join(BASEDIR, "ASF.json")
with open(output_path, "w") as f:
    json.dump(Database, f, indent=None)  # Indent: 1 for readability, None for compact.

# Compress the database file
with open(output_path, "rb") as f_in:
    with gzip.open(output_path + ".gz", "wb") as f_out:
        f_out.writelines(f_in)
