"""
Conversions functions between atomic scattering factors (ASF) and various other formats.

The `conversions` module provides functions for data conversion between
atomic scattering factors (ASF) and various formats including:
- photoabsorption data, NEXAFS data, XANES data,
- refractive component values (index of refraction),
- atomic scattering polynomial (ASP) coefficients.
"""

# Standard library imports
from typing import overload
import warnings

# External imports
import numpy as np
import numpy.typing as npt
import scipy.constants as sc
from scipy.constants import (
    Avogadro as N_A,
    speed_of_light as c,
    Planck as h,  # 6.62e-34
    elementary_charge as e,
    pi,
    epsilon_0,
    electron_mass as m_e,
)

# Internal imports
from kkcalc2.stoich import stoichiometry as kk_stoichiometry

E_RADIUS: float
"""Classical electron radius in meters. ~2.818e-15"""
try:
    E_RADIUS = sc.value("classical electon radius")
except KeyError:
    E_RADIUS: float = 1 / (4 * pi * epsilon_0) * e**2 / (m_e * c**2)


@staticmethod
def energy_to_wavelength(
    energies: npt.NDArray[np.float64 | np.int_] | float | int,
) -> npt.NDArray[np.float64] | float:
    """
    Convert photon energies in eV to wavelengths in Angstroms.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.

    Returns
    -------
    npt.NDArray | float
        Wavelengths in Angstroms.
    """
    return h * c / (energies * e) * 1e10


@overload
def refractive_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    refractive_component: npt.NDArray[np.float64 | np.int_],
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def refractive_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    refractive_component: npt.NDArray[np.complex128],
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> npt.NDArray[np.complex128]: ...  # numpydoc ignore=GL08


@overload
def refractive_to_ASF(
    energies: float | int,
    refractive_component: float | int,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> float: ...  # numpydoc ignore=GL08


@overload
def refractive_to_ASF(
    energies: float | int,
    refractive_component: complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> complex: ...  # numpydoc ignore=GL08


@overload
def refractive_to_ASF(
    energies: npt.NDArray | float | int,
    refractive_component: npt.NDArray | float | int | complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> npt.NDArray | float | complex: ...  # numpydoc ignore=GL08


def refractive_to_ASF(
    energies: npt.NDArray | int | float,
    refractive_component: (
        npt.NDArray[np.int_ | np.float64 | np.complex128] | int | float | complex
    ),
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
    reverse: bool = False,
) -> npt.NDArray[np.floating | np.complexfloating] | int | float | complex:
    r"""
    Convert refractive component values (index of refraction) to atomic scattering factors (ASF).

    .. math::
        n(E) &= 1 - \delta + i\beta
                &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-f_1 + i f_2\right)
                &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-(f^0 + f^') + i f^{''}\right)

    This function applies the scaling above scaing factor, where $n_a$ is the number density,
    $r_e$ is the classical electron radius, $\lambda$ is the wavelength.

    The refractive value is either
    - the imaginary part of the index of refraction, representing absorption.
    - the real part of the index of refraction, representing dispersion.

    Requires some form of material density information to convert to ASF.
    As per positional argument order, the function will use the first available density information.
    This can either be:
    - `number_density` in atoms per millilitre (cm^3), or
    - `density` in grams per millilitre (cm^3), and
        - `formula_mass` (molecular mass), or
        - `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.
    refractive_component : array_like
        Refractive ($\delta$, $\beta$) part of the index of refraction.
    number_density : float, optional
        Material density in atoms per millilitre (cm^3).
    density : float, optional
        Material density in grams per millilitre (cm^3).
    formula_mass : float, optional
        Atomic mass sum of the materials chemical formula (molecular mass).
        Equivalent to providing a `stoichiometry`.
    stoichiometry : stoichiometry | str
        Description of the combination of elements composing the material.
    reverse : bool
        Flag to indicate the reverse conversion.

    Returns
    -------
    npt.NDArray[np.float64 | np.complex128] | float | complex
        Atomic scattering factors.
    """

    # Get number density from material density information.
    if number_density:
        pass
    elif density:
        fm = None
        if formula_mass:
            fm = formula_mass
        elif stoichiometry:
            stoichiometry = (
                kk_stoichiometry(stoichiometry)
                if isinstance(stoichiometry, str)
                else stoichiometry
            )
            if isinstance(stoichiometry, kk_stoichiometry):
                fm = stoichiometry.formula_mass
            else:
                raise ValueError("Invalid stoichiometry provided.")
        else:
            raise ValueError(
                "Material `formula_mass` or `stoichiometry` required with `density` required to convert Beta to ASF."
            )
        number_density = density * N_A / fm
    else:
        raise ValueError(
            "No material density information provided to convert Beta to ASF."
        )

    # Prefactor is e_radius * lambda^2 * n / (2 * pi)
    prefactor = (
        1e-6  # Convert from m^3 to cm^3
        * 2
        * pi
        * energies**2  # lambda = hc/E
        / (number_density * E_RADIUS * (h / e * c) ** 2)
    )
    if not reverse:
        # Generate factors from refraction data.
        factors = prefactor * refractive_component
        return factors
    else:
        factors = refractive_component  # relabel refractive_component as factors.
        # Generate refraction data from factors data.
        refractive_component_reverse = factors / prefactor
        return refractive_component_reverse


@overload
def ASF_to_refractive(
    energies: npt.NDArray[np.float64 | np.int_],
    factors: npt.NDArray[np.float64 | np.int_],
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def ASF_to_refractive(
    energies: npt.NDArray[np.float64 | np.int_],
    factors: npt.NDArray[np.complex128],
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray[np.complex128]: ...  # numpydoc ignore=GL08


@overload
def ASF_to_refractive(
    energies: float | int,
    factors: float | int,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> float: ...  # numpydoc ignore=GL08


@overload
def ASF_to_refractive(
    energies: float | int,
    factors: complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> complex: ...  # numpydoc ignore=GL08


@overload
def ASF_to_refractive(
    energies: npt.NDArray | float | int,
    factors: npt.NDArray | float | int | complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray | float | complex: ...  # numpydoc ignore=GL08


def ASF_to_refractive(
    energies: npt.NDArray[np.float64 | np.int_] | float | int,
    factors: (
        npt.NDArray[np.float64 | np.int_ | np.complex128] | float | int | complex
    ),
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray[np.float64 | np.complex128] | float | complex:
    r"""
    Convert atomic scattering factors (ASF) to refractive component values (index of refraction).

    Uses `refractive_to_ASF` with the `reverse` flag set to `True`.

    .. math::
        n(E) &= 1 - \delta + i\beta
                &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-f_1 + i f_2\right)
                &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-(f^0 + f^') + i f^{''}\right)

    This function applies the scaling above scaing factor, where $n_a$ is the number density,
    $r_e$ is the classical electron radius, $\lambda$ is the wavelength.

    The refractive value is either
    - the imaginary part of the index of refraction, representing absorption.
    - the real part of the index of refraction, representing dispersion.

    Requires some form of material density information to convert from ASF.
    As per positional argument order, the function will use the first available density information.
    This can either be:
    - `number_density` in atoms per millilitre (cm^3), or
    - `density` in grams per millilitre (cm^3), and
        - `formula_mass` (molecular mass), or
        - `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.
    factors : array_like
        Atomic scattering factors.
    number_density : float, optional
        Material density in atoms per millilitre (cm^3).
    density : float, optional
        Material density in grams per millilitre (cm^3).
    formula_mass : float, optional
        Atomic mass sum of the materials chemical formula (molecular mass).
        Equivalent to providing a `stoichiometry`.
    stoichiometry : stoichiometry | str
        Description of the combination of elements composing the material.

    Returns
    -------
    npt.NDArray[np.float64 | np.complex128] | float | complex
        The refractive component value(s) (index of refraction).
    """
    return refractive_to_ASF(
        energies=energies,
        refractive_component=factors,
        number_density=number_density,
        density=density,
        formula_mass=formula_mass,
        stoichiometry=stoichiometry,
        reverse=True,
    )


@overload
def NEXAFS_to_ASF(
    energies: npt.NDArray, NEXAFS: npt.NDArray, reverse: bool = False
) -> npt.NDArray: ...  # numpydoc ignore=GL08


@overload
def NEXAFS_to_ASF(
    energies: int | float, NEXAFS: int | float, reverse: bool = False
) -> float: ...  # numpydoc ignore=GL08


@overload
def NEXAFS_to_ASF(
    energies: npt.NDArray | int | float,
    NEXAFS: npt.NDArray | float,
    reverse: bool = False,
) -> npt.NDArray | float: ...  # numpydoc ignore=GL08


def NEXAFS_to_ASF(
    energies: npt.NDArray | float | int,
    NEXAFS: npt.NDArray | float | int,
    reverse: bool = False,
) -> npt.NDArray | float:
    r"""
    Convert NEXAFS photoabsorption data to atomic scattering factors (ASF).

    This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
    $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

    .. math::
        \mu_a = 2 r_e \lambda f_2

    We relabel the cross section as NEXAFS, the scattering factor $f_2$ (or $f''$) as $ASF_i$ and replace the
    wavelength $\lambda$ for the photon energy $E_i$:

    .. math::
        ASF_i = \frac{1}{2 r_e}\frac{e c}{h} E_i \text{NEXAFS}_i

    Where $r_e$ is the classical electron radius (in meters), $e$ is the elementary charge (in Coulombs),
    $c$ is the speed of light (in m/s), and $h$ is Planck's constant (in J.s) and $E_i$ is the photon energy (in eV).
    The prefactor is approximately 1/(6.9876e-21 J.s^2/C).

    Parameters
    ----------
    energies : array_like
        The photon energies in eV.
    NEXAFS : array_like
        The NEXAFS photoabsorption data.
    reverse : bool, optional
        Flag to indicate the reverse conversion from atomic scattering factors to NEXAFS data.

    Returns
    -------
    npt.NDArray | float
        Scaled data.
    """
    # Convert plank constant from J.s to eV.s via elementary charge.
    prefactor = 2 * E_RADIUS * h / e * c  # ~6.9876e-21 eV s^2

    if not reverse:
        # Convert from NEXAFS to ASF.
        factors = NEXAFS * energies / prefactor
        return factors
    else:
        factors = NEXAFS  # relabel NEXAFS as factors.
        # Convert from ASF to NEXAFS.
        nexafs_reverse = prefactor * factors / energies
        return nexafs_reverse


@overload
def ASF_to_NEXAFS(
    energies: npt.NDArray[np.float64 | np.int_],
    factors: npt.NDArray[np.float64 | np.int_ | np.complex128],
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def ASF_to_NEXAFS(
    energies: int | float, factors: int | float | complex
) -> float: ...  # numpydoc ignore=GL08


@overload
def ASF_to_NEXAFS(
    energies: npt.NDArray | int | float,
    factors: npt.NDArray | int | float | complex,
) -> npt.NDArray | float: ...  # numpydoc ignore=GL08


def ASF_to_NEXAFS(
    energies: npt.NDArray | float | int,
    factors: npt.NDArray | float | int | complex,
) -> npt.NDArray[np.float64] | float:
    """
    Convert atomic scattering factors (ASF) to NEXAFS photoabsorption data.

    Uses `NEXAFS_to_ASF` with the `reverse` flag set to `True`.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.
    factors : array_like
        Atomic scattering factors.

    Returns
    -------
    np.ndarray
        NEXAFS photoabsorption data.
    """
    # Use the absorption (imaginary) part of the atomic scattering factors for NEXAFS / XANES.
    if isinstance(factors, complex) or (
        isinstance(factors, np.ndarray) and factors.dtype is np.complex128
    ):
        factors = factors.imag
    return NEXAFS_to_ASF(energies, factors, reverse=True)


@overload
def refractive_component_to_NEXAFS(
    energies: npt.NDArray[np.float64 | np.int_],
    refractive_component: npt.NDArray[np.float64 | np.int_ | np.complex128],
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def refractive_component_to_NEXAFS(
    energies: float | int,
    refractive_component: float | int | complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> float: ...  # numpydoc ignore=GL08


@overload
def refractive_component_to_NEXAFS(
    energies: npt.NDArray | float | int,
    refractive_component: npt.NDArray | float | int | complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray | float: ...  # numpydoc ignore=GL08


def refractive_component_to_NEXAFS(
    energies: npt.NDArray | float | int,
    refractive_component: npt.NDArray | float | int | complex,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray[np.float64] | float:
    r"""
    Convert $\beta$ values (index of refraction) to NEXAFS photoabsorption data.

    Uses `refractive_component_to_asf` and `ASF_to_NEXAFS` to perform the conversion.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.
    refractive_component : array_like
        Imaginary part of the index of refraction.
    number_density : float, optional
        Material density in atoms per millilitre (cm^3).
    density : float
        Material density in grams per millilitre (cm^3).
    formula_mass : float
        Atomic mass sum of the materials chemical formula (molecular mass).
        Equivalent to providing a `stoichiometry`.
    stoichiometry : stoichiometry | str
        Description of the combination of elements composing the material.

    Returns
    -------
    npt.NDArray
        The NEXAFS photoabsorption equivalent data.
    """
    factors = refractive_to_ASF(
        energies,
        refractive_component,
        number_density,
        density,
        formula_mass,
        stoichiometry,
    )
    # Reduce factors to imaginary component for conversion to NEXAFS / XANES.
    if isinstance(factors, complex) or (
        isinstance(factors, np.ndarray) and factors.dtype is np.complex128
    ):
        factors = factors.imag
    return ASF_to_NEXAFS(energies, factors)


@overload
def NEXAFS_to_refractive_component(
    energies: npt.NDArray,
    NEXAFS: npt.NDArray,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray: ...  # numpydoc ignore=GL08


@overload
def NEXAFS_to_refractive_component(
    energies: float | int,
    NEXAFS: float | int,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> float: ...  # numpydoc ignore=GL08


@overload
def NEXAFS_to_refractive_component(
    energies: npt.NDArray | float | int,
    NEXAFS: npt.NDArray | float | int,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray | float: ...  # numpydoc ignore=GL08


def NEXAFS_to_refractive_component(
    energies: npt.NDArray | float | int,
    NEXAFS: npt.NDArray | float | int,
    number_density: float | None = None,
    density: float | None = None,
    formula_mass: float | None = None,
    stoichiometry: kk_stoichiometry | str | None = None,
) -> npt.NDArray | float:
    r"""
    Convert NEXAFS photoabsorption data to $\beta$ values (index of refractive scale).

    Uses `NEXAFS_to_ASF` and `ASF_to_refractive_component` to perform the conversion.

    Parameters
    ----------
    energies : array_like
        Photon energies in eV.
    NEXAFS : array_like
        NEXAFS photoabsorption data.
    number_density : float, optional
        Material density in atoms per millilitre (cm^3).
    density : float
        Material density in grams per millilitre (cm^3).
    formula_mass : float
        Atomic mass sum of the materials chemical formula (molecular mass).
        Equivalent to providing a `stoichiometry`.
    stoichiometry : stoichiometry | str
        Description of the combination of elements composing the material.

    Returns
    -------
    npt.NDArray
        The refractive component value(s) ($\beta$).
    """
    factors = NEXAFS_to_ASF(energies, NEXAFS)
    return ASF_to_refractive(  # type: ignore
        energies=energies,
        factors=factors,
        number_density=number_density,
        density=density,
        formula_mass=formula_mass,
        stoichiometry=stoichiometry,
    )


@overload
def ASF_to_ASP(
    energies: npt.NDArray[np.float64 | np.int_],
    factors: npt.NDArray[np.float64 | np.int_],
    N: int = 5,
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def ASF_to_ASP(
    energies: npt.NDArray[np.float64 | np.int_],
    factors: npt.NDArray[np.complex128],
    N: int = 5,
) -> npt.NDArray[np.complex128]: ...  # numpydoc ignore=GL08


def ASF_to_ASP(
    energies: npt.ArrayLike | npt.NDArray[np.float64 | np.int_],
    factors: npt.ArrayLike | npt.NDArray[np.float64 | np.int_ | np.complex128],
    N: int = 5,
) -> npt.NDArray[np.float64 | np.complex128] | float | complex:
    """
    Convert atomic scattering factors (ASF) to atomic scattering polynomial (ASP) coefficients.

    Calculates `N` polynomial coefficients for the spans between `factors` defined at `energies`.
    Currently only the first two coefficients are calculated (linear interpolation).

    Parameters
    ----------
    energies : array_like
        An array of `N` photon energies in eV.
    factors : array_like
        An array of `N` atomic scattering factors.
    N : int
        The number of coefficients to calculate. Default is 5.

    Returns
    -------
    npt.NDArray
        An array of `N-1` Atomic scattering polynomial coefficients.
    """
    energies = np.asarray(energies)
    factors = np.asarray(factors)
    if isinstance(energies, np.ndarray):
        # Ensure no duplicate energies and ordered.
        diffs = np.diff(energies)
        monotonic = np.all(diffs > 0)
        if not monotonic:
            raise ValueError(
                "Energies must be ordered and unique."
                + f" Negative differences: {np.where(diffs <= 0)}."
            )
        # Calculate the coefficients: setup array.
        coefs = np.zeros((len(energies) - 1, N))
        # Calculate coefficient #0
        coefs[:, 0] = (factors[1:] - factors[:-1]) / (energies[1:] - energies[:-1])
        # Calculate coefficient #1
        coefs[:, 1] = factors[:-1] - coefs[:, 0] * energies[:-1]
    return coefs


@overload
def ASP_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    coefs: npt.NDArray[np.float64 | np.int_],
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def ASP_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    coefs: npt.NDArray[np.complex128],
) -> npt.NDArray[np.complex128]: ...  # numpydoc ignore=GL08


@overload
def ASP_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    coefs: npt.NDArray[np.float64 | np.int_],
    orders: npt.NDArray[np.integer] | None,
) -> npt.NDArray[np.float64]: ...  # numpydoc ignore=GL08


@overload
def ASP_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    coefs: npt.NDArray[np.complex128],
    orders: npt.NDArray[np.integer] | None,
) -> npt.NDArray[np.complex128]: ...  # numpydoc ignore=GL08


def ASP_to_ASF(
    energies: npt.NDArray[np.float64 | np.int_],
    coefs: npt.NDArray[np.float64 | np.int_ | np.complex128],
    orders: npt.NDArray[np.integer] | None = None,
) -> npt.NDArray[np.float64 | np.complex128]:
    """
    Convert the atomic scattering polynomial (ASP) coefficients to atomic scattering factors (ASF).

    Parameters
    ----------
    energies : array_like
        An array of `N` or `N+1` photon energies in eV.
        using the starting energy of each interval.
    coefs : array_like
        An array with dimension (`N`, `M`), with `N` sets of `M` atomic
        scattering polynomial coefficients.
    orders : npt.NDArray[np.int_] | None, optional
        An array of `M` integers defining the polynomial orders for each coefficient set.
        Each integer corresponds to the power of the energy term multipled by
        the corresponding coefficient in the  polynomial, before summation to factors.
        If None (default), assumes coefficients in sequential order decreasing from 1.
        i.e: [1, 0, -1, -2, ...] etc.

    Returns
    -------
    npt.NDArray
        An array of `N` or `N+1` atomic scattering factors, matching the input `energies` length.
        If `energies` has length `N+1`, the last ASF value will be calculated using the last ASP coefficient.
    """
    # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
    energies = np.asarray(energies, dtype=float)
    coefs = np.asarray(coefs)
    # Check dimensions:
    if energies.ndim == 0:
        # Boost to 1D
        energies = np.array([energies])
    if coefs.ndim == 1:
        warnings.warn("Single coefficient set provided, boosting to 2D.")
        # Boost to 2D
        coefs = np.array([coefs])

    # Check shapes:
    if (coefs.shape[0] != energies.shape[0] - 1) and (
        coefs.shape[0] != energies.shape[0]
    ):
        raise ValueError(
            f"Number of coefficients sets ({len(coefs)}) "
            + f"does not match the number of energies ({len(energies) - 1} or {len(energies)})."
        )

    # Create an array of energy powers for each coefficient.
    if orders is not None:
        if orders.ndim == 1 and orders.shape[0] == coefs.shape[1]:
            powers = np.c_[*[energies**i for i in orders]]
        else:
            raise ValueError(
                f"Number of orders ({orders.shape[0]}) "
                + f"does not match the number of coefficients ({coefs.shape[1]})."
            )
    else:
        powers = np.c_[*[energies ** (1 - i) for i in range(5)]]
    # Do energies match the number of coefficient sets?
    if energies.shape[0] == coefs.shape[0] + 1:
        # Duplicate the final polynomial to define the final boundary.
        coefs = np.r_[coefs, coefs[-1:, :]]  # Duplicate the last row.
    return np.sum(coefs * powers, axis=1)
