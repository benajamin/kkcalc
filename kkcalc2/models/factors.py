"""
'Atomic scattering factor' data models.

Defines the types of data that can be used, and conversion between.
"""

# In polynomials.py, the equivalent import is only done via type checking or in functions, to prevent recursion.
from kkcalc2.models.polynomials import (
    asp as asp_type,
    asp_abstract,
    asp_im,
    asp_re,
    asp_complex,
)

## ..
from kkcalc2.models.common import (
    atomic_scattering_abstract,
    atomic_scattering,
    PROPERTIES_DICT,
    PROPERTIES_DICT_NO_STOICH,
)
from kkcalc2 import conversions
from kkcalc2.stoich import (
    stoichiometry as kk_stoichiometry,
)  # To prevent overlap use with the `stoichiometry` argument.
from kkcalc2.transforms import DEF_ITER, DEF_TOL

import numpy as np
import numpy.typing as npt
import abc
import warnings
from enum import Enum
from typing import Self, Iterator, override, overload, Unpack, TypedDict

try:
    import pandas as pd

    has_pandas = True
except ImportError:
    has_pandas = False

KK_DATATYPE_DOCS: dict[str, str] = {
    "UNDEFINED": """For undefined data types.""",
    "NEXAFS": """Near edge X-ray absorption fine structure (NEXAFS).""",
    "XANES": """X-ray absorption near edge structure (XANES).""",
    "PHOTOABSORPTION": """Photoabsorption.""",
    "REFRACTIVE": r"""Refractive components, with dispersive :math:`\delta` and absorptive :math:`\beta` components.

    Not to be confused with the index of refraction (KK_Datatype.REFRACTIVE_INDEX).

    .. math::
        n(E) = \delta(E) + i\beta(E)
    """,
    "REFRACTIVE_INDEX": r"""
    Index of refraction, with dispersive :math:`\delta` and absorptive :math:`\beta` components.

    .. math::
        n(E) = 1 - \delta(E) - i\beta(E)
    """,
    "ASF": r"""
    Atomic scattering factors, real :math:`f_1` & imaginary :math:`f_2` components.
    Both scattering strengths are relative to the Thompson scattering of a free electron.
    Calculated for a set of :math:`m` elements, with number density :math:`N_m`,
    wavelength :math:`\lambda`, photon energy :math:`E` and classical electron radius :math:`r_0`.

    .. math::
        n(E) = 1 - \frac{r_0}{2\pi}\lambda^2\sum_m \left(f_{1m}(E) + i f_{2m}(E)\right)

    See Also:
    kkcalc.stoich.relativistic_correction_eq : The sum of relativistic corrections for a elemental composition.
    """,
    "ASF_DASH": r"""
    Atomic scattering factors, real :math:`f^{0}`, :math:`f^{'}` and imaginary :math:`f^{''}` components.
    Here the relativistic correction is :math:`f^{0}`, and :math:`f^{'}` is the energy dependent real part.
    Both scattering strengths are relative to the Thompson scattering of a free electron.

    See Also
    --------
    kkcalc.stoich.relativistic_correction_eq : The relativistic correction for a elemental composition.
    """,
}


class KK_Datatype(Enum):
    """
    Enum for the type of data to be used in the Kramers-Kronig calculation.
    """

    UNDEFINED = 0
    """For undefined data types."""
    NEXAFS = 1  # AKA Photoabsorption, XANES.
    """Near edge X-ray absorption fine structure (NEXAFS)."""
    XANES = 1  # AKA Photoabsorption, NEXAFS.
    """X-ray absorption near edge structure (XANES)."""
    PHOTOABSORPTION = 1  # AKA NEXAFS, XANES.
    """Photoabsorption."""
    REFRACTIVE = 2
    r"""Refractive components, with dispersive :math:`\delta` and absorptive :math:`\beta` components.

    Not to be confused with the index of refraction (KK_Datatype.REFRACTIVE_INDEX).

    .. math::
        n(E) = \delta(E) + i\beta(E)
    """
    REFRACTIVE_INDEX = 3  # Index of refraction
    r"""
    Index of refraction, with dispersive :math:`\delta` and absorptive :math:`\beta` components.

    .. math::
        n(E) = 1 - \delta(E) - i\beta(E)
    """
    ASF = 4  # Atomic scattering factors as per the original KK Calc; f1 & f2,
    r"""
    Atomic scattering factors, real :math:`f_1` & imaginary :math:`f_2` components.
    Both scattering strengths are relative to the Thompson scattering of a free electron.
    Calculated for a set of :math:`m` elements, with number density :math:`N_m`,
    wavelength :math:`\lambda`, photon energy :math:`E` and classical electron radius :math:`r_0`.

    .. math::
        n(E) = 1 - \frac{r_0}{2\pi}\lambda^2\sum_m \left(f_{1m}(E) + i f_{2m}(E)\right)

    See Also:
    kkcalc.stoich.relativistic_correction_eq : The sum of relativistic corrections for a elemental composition.
    """
    ASF_DASH = 5  # Atomic scattering factors f0, f' & f'',
    r"""
    Atomic scattering factors, real :math:`f^{0}`, :math:`f^{'}` and imaginary :math:`f^{''}` components.
    Here the relativistic correction is :math:`f^{0}`, and :math:`f^{'}` is the energy dependent real part.
    Both scattering strengths are relative to the Thompson scattering of a free electron.

    See Also
    --------
    kkcalc.stoich.relativistic_correction_eq : The relativistic correction for a elemental composition.
    """


#
for i, dtype in enumerate(KK_Datatype):
    name = dtype.name.upper()
    if name in KK_DATATYPE_DOCS:
        dtype.__doc__ = KK_DATATYPE_DOCS[name]


class KK_ASF_DICT(TypedDict, total=False):
    """Initialisation kwargs of the atomic scattering factor classes."""

    origin_dtype: KK_Datatype
    """Original data type of the atomic scattering factors."""
    origin_data: npt.NDArray
    """Original data of the atomic scattering factors."""


class asf_abstract(atomic_scattering_abstract, metaclass=abc.ABCMeta):
    """
    Abstract implementation of the atomic scattering factor class.

    Provides the interface for the atomic scattering factor classes, by
    defining the required properties and method signatures.

    Parameters
    ----------
    energies : npt.NDArray
        Beam energies in eV.
    factors : npt.NDArray
        Atomic scattering factors.
    **kwargs : Unpack[KK_ASF_DICT]
        Additional keyword arguments for the `atomic_scattering` base class,
        and the `KK_ASF_DICT` dictionary (e.g. `origin_dtype`, `origin_data`).
    """

    @abc.abstractmethod
    def __init__(
        self,
        energies: npt.ArrayLike,
        factors: npt.ArrayLike,
        **kwargs: Unpack[
            KK_ASF_DICT
        ],  # TODO: Merge with PROPERTIES_DICT if possible, and implement KK_ASF_DICT across other subclasses...
    ) -> None:  # numpydoc ignore=GL08
        raise NotImplementedError("Cannot instantiate abstract class.")

    @property
    @abc.abstractmethod
    def energies(self) -> npt.NDArray:  # numpydoc ignore=PR02
        """
        Property for the energies of the atomic scattering factors.

        Parameters
        ----------
        energies : npt.NDArray
            Energies in eV.

        Returns
        -------
        npt.NDArray
            Energies in eV.
        """
        raise NotImplementedError("Getter requires implementation.")

    @energies.setter
    @abc.abstractmethod
    def energies(self, energies: npt.ArrayLike) -> None:  # numpydoc ignore=GL08
        raise NotImplementedError("Setter requires implementation.")

    @property
    @abc.abstractmethod
    def factors(self) -> npt.NDArray:  # numpydoc ignore=PR02
        """
        Property for the atomic scattering factors.

        Parameters
        ----------
        factors : array_like
            Atomic scattering factors.

        Returns
        -------
        np.ndarray
            Atomic scattering factors.
        """
        raise NotImplementedError("Getter requires implementation.")

    @factors.setter
    @abc.abstractmethod
    def factors(self, factors: npt.ArrayLike) -> None:  # numpydoc ignore=GL08
        raise NotImplementedError("Setter requires implementation.")

    @property
    def data(self) -> tuple[npt.NDArray, npt.NDArray]:  # numpydoc ignore=PR02
        """
        Property for the atomic scattering factor data (energies and amplitudes).

        Parameters
        ----------
        data : tuple[npt.NDArray, npt.NDArray]
            Tuple of energies (eV) and atomic scattering factors.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Tuple of energies (eV) and atomic scattering factors.
        """
        return self.energies, self.factors

    @data.setter
    def data(
        self, data: tuple[npt.NDArray, npt.NDArray]
    ) -> None:  # numpydoc ignore=GL08
        if (
            not isinstance(data, tuple)
            or len(data) != 2
            or len(data[0]) != len(data[1])
        ):
            raise ValueError("Data must be a tuple of two equal length arrays.")
        self.energies, self.factors = np.asarray(data[0]), np.asarray(data[1])

    @property
    def refractive(self) -> np.ndarray:
        r"""
        Scale the atomic scattering factors to $\delta$ and $\beta$ index-of-refraction values.

        This method converts the atomic scattering factors to refractive values by using the
        `conversions.ASF_to_refractive` method. The factors are scaled to 'index of refraction' `n(E)`
        components, and can consist of only the real ($delta$) or imaginary ($beta$) part, or combined complex value.
        Note, this is NOT equivalent to the 'index of refraction' of a material.

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-f_1 + i f_2\right)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-(f^0 + f^') + i f^{''}\right)

        Where
        - $\\beta$ is the imaginary part of the index of refraction, representing absorption.
        - $1-\\delta$ is the real part of the index of refraction, where $\\delta$ represents dispersion.
        - $f_1$ is the real part of the atomic scattering factor and includes the elemental electron density.
        - $f_2$ is the imaginary part of the atomic scattering factor.

        Requires some form of material density information to convert to ASF.
        This can either be:
        - `number_density` in atoms per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and
            - `formula_mass` (molecular mass), or
            - `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Each of these parameters can be provided as an attribute for the object.

        Returns
        -------
        np.ndarray
            A component of the refractive index, either real or imaginary, or complex,
            depending on the factor data.

        Raises
        ------
        ValueError
            If the object attributes do not provide enough information to calculate the number density.

        See Also
        --------
        kkcalc2.models.factors.asf_abstract.to_refractive : Method
            Converts the atomic scattering factors to refractive values, with density parameter arguments.
        """
        if self.number_density is not None:
            return conversions.ASF_to_refractive(
                energies=self.energies,
                factors=self.factors,
                number_density=self.number_density,
            )
        elif self.density is not None:
            if self.formula_mass is not None:
                return conversions.ASF_to_refractive(
                    energies=self.energies,
                    factors=self.factors,
                    density=self.density,
                    formula_mass=self.formula_mass,
                )
            elif self.stoichiometry is not None:
                return conversions.ASF_to_refractive(
                    energies=self.energies,
                    factors=self.factors,
                    density=self.density,
                    stoichiometry=self.stoichiometry,
                )
        raise ValueError(
            "Material density information is required to convert to 'index of refraction' values."
        )

    def to_refractive(
        self,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
    ) -> np.ndarray:
        """
        Same as `refractive` property, but allows specification of density information.

        Parameters
        ----------
        number_density : float
            Number density of the material in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a stoichiometry.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.

        Returns
        -------
        np.ndarray
            'index of refraction' values corresponding to `energies` property.

        See Also
        --------
        kkcalc2.models.factors.asf_abstract.refractive : Property
            Converts the atomic scattering factors to refractive values.
        """
        if (
            number_density is None
            and density is None
            and formula_mass is None
            and stoichiometry is None
        ):
            # Attempt to use the object's density information to convert to Beta values.
            number_density = self.number_density
            density = self.density
            formula_mass = self.formula_mass
            stoichiometry = self.stoichiometry
        else:
            if stoichiometry is not None and isinstance(stoichiometry, str):
                stoichiometry = kk_stoichiometry(stoichiometry)

        # Attempt to use available density information to convert to Beta values.
        if number_density is not None:
            return conversions.ASF_to_refractive(
                energies=self.energies,
                factors=self.factors,
                number_density=number_density,
            )
        elif density is not None:
            if formula_mass is not None:
                return conversions.ASF_to_refractive(
                    energies=self.energies,
                    factors=self.factors,
                    density=density,
                    formula_mass=formula_mass,
                )
            elif stoichiometry is not None:
                return conversions.ASF_to_refractive(
                    energies=self.energies,
                    factors=self.factors,
                    density=density,
                    stoichiometry=stoichiometry,
                )
        raise ValueError(
            "Material density information is required to convert to 'index of refraction' values."
        )

    @overload
    @classmethod
    def from_refractive(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        number_density: float,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        density: float,
        formula_mass: float,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        density: float,
        stoichiometry: kk_stoichiometry | str,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @classmethod
    @abc.abstractmethod
    def from_refractive(
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Load refractive values to as atomic scattering factors.

        Uses some form of material density information to convert to atomic scattering factors.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3). By default, None.
        density : float, optional
            Material density in grams per millilitre (cm^3). By default, None.
        formula_mass : float, optional
            Atomic mass sum of the materials chemical formula (molecular mass). By default, None.
        stoichiometry : stoichiometry | str, optional
            Description of the combination of elements composing the material. By default, None.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale. By default, False.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` base class.

        Returns
        -------
        np.ndarray
            Refractive values corresponding to the `energies` property.
        """
        raise NotImplementedError("Must be implemented in subclass.")

    @overload
    @classmethod
    def from_refractive_index(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        number_density: float,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive_index(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        density: float,
        formula_mass: float,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive_index(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        density: float,
        stoichiometry: kk_stoichiometry | str,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @overload
    @classmethod
    def from_refractive_index(  # numpydoc ignore=GL08
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self: ...

    @classmethod
    @abc.abstractmethod
    def from_refractive_index(
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Abstract method to convert refractive index values.

        Convert real or imaginary refractive index values (1 - $\delta$ or $\beta$) to atomic scattering factors (ASF).

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive_index : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf
            Atomic scattering factors equivalent representation.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        raise NotImplementedError("Must be implemented in subclass.")

    @staticmethod
    def atomic_scattering_factors_to_coefs(
        energies: npt.NDArray, factors: npt.NDArray, N: int = 5
    ) -> npt.NDArray:
        """
        Convert atomic scattering factors (ASF) to atomic scattering polynomial (ASP) coefficients.

        Alias for `conversions.ASF_to_ASP` to calculate the atomic scattering polynomial coefficients from
        atomic scattering `factors` defined at `energies`.

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
        return conversions.ASF_to_ASP(energies, factors)

    @property
    def atomic_scattering_polynomial(self) -> npt.NDArray:
        """
        Convert atomic scattering factors to atomic scattering polynomial coefficients.

        Uses `energies` and `factors` with length `N` to calculate the atomic scattering polynomial coefficients.
        To convert to an `asp` object, use the `to_atomic_scattering_polynomial` method.

        Returns
        -------
        npt.NDArray
            An array with dimension (`N-1`, 5) of atomic scattering polynomial coefficients.
        """
        return self.atomic_scattering_factors_to_coefs(self.energies, self.factors)

    # @property
    # @doc_copy(atomic_scattering_polynomial)
    # def asp(self) -> npt.NDArray:
    #     """
    #     Alias for `atomic_scattering_polynomial`.
    #     """
    #     return self.atomic_scattering_polynomial
    asp = atomic_scattering_polynomial  # Alias for atomic scattering polynomial coefficients.

    @abc.abstractmethod
    def to_atomic_scattering_polynomial(self) -> asp_abstract:
        """
        Convert factor representation to polynomial representation.

        Uses the `energies` and `factors` attributes (with length `N`) of the object
        to create an atomic scattering polynomial object of coefficients with length (N-1).

        For an array of polynomial coefficients, use the `atomic_scattering_polynomial` property.

        Returns
        -------
        asp
            Atomic scattering polynomial object.
        """
        pass

    # @abc.abstractmethod
    # @doc_copy(to_atomic_scattering_polynomial)
    # def to_ASP(self) -> asp_abstract:
    #     """
    #     Alias for `to_atomic_scattering_polynomial`.
    #     """
    #     return self.to_atomic_scattering_polynomial()
    to_ASP = to_atomic_scattering_polynomial  # Alias for atomic scattering polynomial conversion.

    def dataframe(self) -> "pd.DataFrame":
        """
        Generate a Pandas representation of the factors list, useful for display.

        Returns
        -------
        pd.DataFrame
            A Pandas DataFrame with the energies and factors.

        Raises
        ------
        ImportError
            If Pandas is not installed.
        """
        if not has_pandas:
            raise ImportError("Pandas is required for this method.")
        return pd.DataFrame(  # type: ignore - `has_pandas` check above.
            np.c_[self.energies, self.factors], columns=["Energy (eV)", "ASF"]
        )

    def __str__(self, **kwargs) -> str:
        """
        Create a string representation of the factors list.

        Pandas is used if availalble to create a table representation.
        Rows displayed are the first and last 5 if more than 10 rows.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `pd.dataFrame.to_string` method.

        Returns
        -------
        str
            A string representation of the factors.
        """
        # Create a default max_rows if not provided.
        if has_pandas:
            if "max_rows" not in kwargs:
                kwargs["max_rows"] = 6
            return self.dataframe().to_string(**kwargs)
        else:
            return (
                f"{self.__class__.__name__} object with {len(self.energies)} energies."
            )

    def __getitem__(self, key: int | slice) -> Self:
        """
        A subset of the atomic scattering factors at the specified index/slice.

        Parameters
        ----------
        key : int | slice
            Index or slice of the atomic scattering factors.

        Returns
        -------
        type[asf_abstract]
            Atomic scattering factors object at the specified index.
        """
        if not isinstance(key, (int, slice)):
            raise TypeError("Index must be an integer or slice.")
        common_kwargs = self._properties_dict
        return self.__class__(
            energies=self.energies[key], factors=self.factors[key], **common_kwargs
        )

    def __iter__(self) -> Iterator[tuple[float, float | complex]]:
        """
        Provide each energy and factor of the energy-dependent scattering amplitude.

        Yields
        ------
        energy : float
            The energy at which the scattering factor is defined.
        factor : float
            The atomic scattering factor.
        """
        for i in range(len(self.energies)):
            yield (self.energies[i], self.factors[i])

    @abc.abstractmethod
    def copy(self, **kwargs) -> "asf_abstract":
        """
        Generate a copy of the `asf` object.

        Parameters
        ----------
        **kwargs
            Any keyword arguments for the constructors to update the copy properties.

        Returns
        -------
        asf_abstract
            A new `asf` object with the same atomic scattering factors,
            and properties, but unique memory allocation.
        """
        pass


class asf(asf_abstract, atomic_scattering):
    """
    Generic class for handling atomic scattering factors.

    Provides static methods to convert from `KK_Datatype` to scattering factors.

    Parameters
    ----------
    energies : np.ndarray
        Energies in eV.
    factors : np.ndarray
        Atomic scattering factors.
        If data y-data is instead betas or NEXAFS, use the respective `from_<name>` method.
    origin_dtype : KK_Datatype, optional
        Original data type of the atomic scattering factors.
        If not provided, the original data is assumed to be the same as the input data.
    origin_data : np.ndarray, optional
        Original data of the atomic scattering factors.
        If not provided, the original data is assumed to be the same as the input data.
    scale_to_database : bool, optional
        Whether to scale the atomic scattering factors to the Henke database scale.
    **kwargs :
        Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` such as:
        - `number_density` : float
        - `density` : float
        - `stoich` : stoichiometry
        - `formula_mass` : float
        - `name` : str

    See Also
    --------
    kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
    """

    def __init__(
        self,
        energies: npt.NDArray,
        factors: npt.NDArray,
        origin_dtype: KK_Datatype | None = None,
        origin_data: npt.NDArray | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Initialise the atomic scattering base class
        atomic_scattering.__init__(self, **kwargs)

        # Initialize hidden attributes
        self._energies: npt.NDArray
        self._factors: npt.NDArray | None
        self._energies = energies = np.asarray(energies)
        self._factors = factors = np.asarray(factors)

        # Check energies are sorted
        if not np.all(np.diff(energies) > 0):
            warnings.warn("Dataset energies are not sorted, sorting.")
            idxs = np.argsort(energies)
            self.energies = energies = energies[idxs]
            self.factors = factors = np.asarray(factors)[idxs]

        if origin_dtype is None:
            # If no original data type is provided, assume the input data is the original data.
            self._origin_dtype = KK_Datatype.ASF
        else:
            self._origin_dtype = origin_dtype
        if origin_data is not None:
            # Store a copy of original data.
            origin_data = np.asarray(origin_data)
            origin_data = origin_data.copy()
            if len(origin_data.shape) != 2:
                raise ValueError("Original data must contain two columns.")
            self._origin_data = origin_data
        elif origin_dtype is None:
            # If no original data is provided, assume the input data is the original data.
            self._origin_data = np.c_[energies, factors]  # already creates copies

        if scale_to_database and self.__class__ == asf:
            # Do not allow the base class to scale, as no complexity designation.
            warnings.warn(
                f"Scaling to database is only available for real and imaginary components, not for {self}. Turning off scaling."
            )
            scale_to_database = False

        if scale_to_database:
            self.scale_to_database()

    @property
    def energies(self) -> np.ndarray:  # numpydoc ignore=PR02
        """
        Property for the atomic scattering factors.

        Parameters
        ----------
        energies : np.ndarray
            Energies in eV.

        Returns
        -------
        np.ndarray
            Energies in eV.
        """
        return self._energies

    @energies.setter
    def energies(self, energies: npt.ArrayLike) -> None:  # numpydoc ignore=GL08
        self._energies = np.asarray(energies)
        if self.factors is not None and len(self._energies) != len(self.factors):
            warnings.warn(
                "Length of energies does not match the length of factors. Factors have been discarded."
            )
            self._factors = None

    @property
    def factors(self) -> np.ndarray:  # numpydoc ignore=PR02
        """
        Property for the atomic scattering factors.

        Parameters
        ----------
        factors : array_like
            Atomic scattering factors.

        Returns
        -------
        np.ndarray
            Atomic scattering factors.

        Raises
        ------
        ValueError
            If the factors have been reset due to a change in energies.
        """
        if self._factors is None:
            raise ValueError(
                "Factors have been reset due to a change in energies, and require setting."
            )
        return self._factors

    @factors.setter
    def factors(self, factors: npt.ArrayLike) -> None:  # numpydoc ignore=GL08
        factors = np.asarray(factors)
        if len(factors) != len(self.energies):
            raise ValueError("Length of factors does not match the length of energies.")
        self._factors = factors

    @property
    def origin_dtype(self) -> KK_Datatype:
        """
        The original data type of the atomic scattering factors.

        Returns
        -------
        KK_Datatype
            Enumerate of the original data type.
        """
        return self._origin_dtype

    @property
    def origin_data(self) -> np.ndarray | None:
        """
        The original data provided for the atomic scattering factors.

        Returns
        -------
        np.ndarray | None
            Original data of the atomic scattering factors, matching
            the format described by the `origin_dtype` attribute.
            Returns a copy.
        """
        if self.origin_dtype == KK_Datatype.UNDEFINED:
            return None
        return self._origin_data.copy()

    def scale_to_database(
        self,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
    ) -> None:
        """
        Scale the data to the Henke database reference.

        If the object contains a stoichiometry and has an identifiable complexity (imag, real),
        this method scales the atomic scattering factors to the database scale.

        Origin data is unmodified, but factors are modified.

        Parameters
        ----------
        merge_domain : tuple[float, float], optional
            Energy range (eV) over which to merge the data to the database.
        fix_distortions : bool, optional
            Whether to attempt to fix distortions in the data before scaling.
            By default, False.

        Raises
        ------
        ValueError
            If the object does not contain a stoichiometry.
            If the object is not a real or imaginary designated complexity.
        """
        if self.stoichiometry is not None:
            if isinstance(self, asf_re):
                from kkcalc2.models.db_models import asp_db_re

                self.factors = asp_db_re.scale_data(
                    self.energies,
                    self.factors,
                    self.stoichiometry,
                    merge_domain,
                    fix_distortions,
                )
                return
            elif isinstance(self, asf_im):
                from kkcalc2.models.db_models import asp_db_im

                self.factors = asp_db_im.scale_data(
                    self.energies,
                    self.factors,
                    self.stoichiometry,
                    merge_domain,
                    fix_distortions,
                )
                return
            raise ValueError(
                f"Scaling to database is only available for real and imaginary components. {self} is {self.__class__}."
            )
        raise ValueError("Scaling to database requires a stoichiometry.")

    def to_atomic_scattering_polynomial(
        self, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> asp_type:
        """
        Convert factors representation to polynomial representation.

        Uses the `energies` and `factors` attributes (with length `N`) of the object
        to create an atomic scattering polynomial object of coefficients with length (N-1).

        For an array of polynomial coefficients, use the `atomic_scattering_polynomial` property.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asp
            Atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Get the existing properties, update with kwargs
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        # Return asp object
        return asp_type(
            energies=self.energies[:-1],
            coefs=self.atomic_scattering_polynomial,
            **common_kwargs,
        )

    # @doc_copy(to_atomic_scattering_polynomial)
    # def to_ASP(self, **kwargs) -> asp_type:
    #     """
    #     Alias for `to_atomic_scattering_polynomial`.
    #     """
    #     return self.to_atomic_scattering_polynomial(**kwargs)

    to_ASP = to_atomic_scattering_polynomial  # Alias for atomic scattering polynomial conversion.

    @override
    @classmethod
    def from_refractive(
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert refractive values ($\delta$ or $\beta$) to atomic scattering factors (ASF).

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf
            Atomic scattering factors equivalent representation.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Convert energy and beta data to numpy arrays.
        energies = np.asarray(energies)
        refractive = np.asarray(refractive)
        # Perform conversion
        factors = conversions.refractive_to_ASF(
            energies, refractive, number_density, density, formula_mass, stoichiometry
        )
        # Accumulate keyword arguments
        kwargs.update(
            {
                "number_density": number_density,
                "density": density,
                "formula_mass": formula_mass,
                "stoichiometry": stoichiometry,
            }
        )
        # Return asf instance
        return cls(
            energies=energies,
            factors=factors,
            origin_dtype=KK_Datatype.REFRACTIVE_INDEX,
            origin_data=np.c_[energies, refractive],
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @override
    @classmethod
    @abc.abstractmethod
    def from_refractive_index(
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Unimplemented method to convert refractive index values.

        Always raises a `NotImplementedError` as the conversion to atomic scattering factors is
        ambiguous between $\delta$ & $\beta$ values.

        Instead use `from_refractive` to scale $\delta$ or $\beta$ values in the same way,
        or use the `asf_re`, `asf_im` or `asf_complex` classes to decompose the refractive index.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive_index : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Raises
        ------
        NotImplementedError
            Always is raised as the conversion is not implemented for this class.
        """
        raise NotImplementedError(
            r"Refractive index conversion is not implemented for this class,\
                as the conversion is ambiguous between $\delta$ & $\beta$ values."
        )

    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Generate a copy of the `asp` object.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Any keyword arguments for the constructors to update the copy properties.

        Returns
        -------
        type[asp]
            A new `asp` object with the same polynomial coefficients,
            and properties, but unique memory allocation.
        """
        # Copy the object properties
        common_kwargs = self._properties_dict
        for key in common_kwargs:
            if hasattr(common_kwargs[key], "copy"):
                common_kwargs[key] = common_kwargs[key].copy()
        # Update with kwargs
        common_kwargs.update(kwargs)
        # Create a new object
        return self.__class__(
            energies=self.energies.copy(),
            factors=self.factors.copy(),
            origin_dtype=self.origin_dtype,
            origin_data=self.origin_data.copy(),
            **common_kwargs,
        )


class asf_re(asf):
    """
    Identical to `asf`, but reserved for real component factors.

    Parameters
    ----------
    energies : npt.NDArray
        Photon energies in eV.
    factors : npt.NDArray
        Atomic scattering factors.
        If data y-data is instead deltas, use the respective `from_<name>` classmethod.
    origin_dtype : KK_Datatype, optional
        Original data type of the atomic scattering factors.
        If not provided, the original data is assumed to be the same as the input data.
    origin_data : npt.NDArray, optional
        Original data of the atomic scattering factors.
    scale_to_database : bool, optional
        Whether to scale the scattering factors to the Henke Database background.
        By default False.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments for the `atomic_scattering` object.
    """

    def __init__(
        self,
        energies: npt.NDArray,
        factors: npt.NDArray,
        origin_dtype: KK_Datatype | None = None,
        origin_data: npt.NDArray | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Initialise the atomic scattering base class
        asf.__init__(
            self,
            energies=energies,
            factors=factors,
            origin_dtype=origin_dtype,
            origin_data=origin_data,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @classmethod
    def from_asf(cls: type[Self], asf: asf, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Designate an `asf` object into real dispersive factors representation.

        Parameters
        ----------
        asf : asf
            Atomic scattering factors object.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_re
            Atomic scattering factors object designated as a real component.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Get the existing properties, update with kwargs
        common_kwargs = asf._properties_dict
        common_kwargs.update(kwargs)
        # Return asf object
        return cls(
            energies=asf.energies,
            factors=asf.factors,
            origin_dtype=asf.origin_dtype,
            origin_data=asf.origin_data,
            **common_kwargs,
        )

    def to_atomic_scattering_polynomial(
        self, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> asp_re:
        """
        Convert real factors to a real polynomial representation.

        Uses the `energies` and `factors` attributes (with length `N`) of the object
        to create an atomic scattering polynomial object of coefficients with length (N-1).

        For an array of polynomial coefficients, use the `atomic_scattering_polynomial` property.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asp
            Atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Get the existing properties, update with kwargs
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        # Return asp object
        return asp_re(
            energies=self.energies,
            coefs=self.atomic_scattering_polynomial,
            **common_kwargs,
        )

    # @doc_copy(to_atomic_scattering_polynomial)
    # def to_ASP(self, **kwargs) -> asp_re:
    #     """
    #     Alias for `to_atomic_scattering_polynomial`.
    #     """
    #     return self.to_atomic_scattering_polynomial(**kwargs)
    to_ASP = to_atomic_scattering_polynomial  # alias for atomic scattering polynomial conversion.

    @property
    def deltas(self) -> npt.NDArray[np.floating]:
        r"""
        Calculate the real dispersive refraction values ($\delta$) from atomic scattering factors.

        This is equivalent to the `refractive` property.

        Returns
        -------
        npt.NDArray[np.floating]
            The real part of the refractive index.
        """
        return self.refractive

    @classmethod
    def from_deltas(
        cls: type[Self],
        energies: npt.NDArray,
        dispersion: npt.NDArray,
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asf_re":
        r"""
        Convert real refractive dispersion values ($\delta$) to atomic scattering factors (ASF).

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        dispersion : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_re
            Real atomic scattering factors representation.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        return asf_re.from_refractive(
            energies=energies,
            refractive=dispersion,
            number_density=number_density,
            density=density,
            formula_mass=formula_mass,
            stoichiometry=stoichiometry,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @override
    @classmethod
    def from_refractive_index(
        cls: type[Self],
        energies: npt.NDArray[np.floating],
        refractive_index: npt.NDArray[np.floating],
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert real refractive index values (n = 1 - $\delta$) to atomic scattering factors (ASF).

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive_index : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf
            Atomic scattering factors equivalent representation.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Convert the real refractive index value to $\delta$ then use `from_deltas`.
        deltas = 1 - refractive_index
        # Convert to atomic scattering factors
        return asf_re.from_deltas(
            energies=energies,
            dispersion=deltas,
            number_density=number_density,
            density=density,
            formula_mass=formula_mass,
            stoichiometry=stoichiometry,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    def kk_transform_inv(
        self,
        target_energies: npt.NDArray[np.floating] | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
    ) -> "asf_im":
        """
        Inverse Kramers-Kronig transform for the real part of the atomic scattering factors.

        Converts `asf_re` to `asp_re` and uses the `KK_PP` method to calculate the inverse transform.

        Parameters
        ----------
        target_energies : npt.NDArray
            The energies at which to calculate the inverse transform.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.

        Returns
        -------
        asf_im
            Atomic scattering factors object with imaginary components.
        """
        # Convert
        asp_re = self.to_atomic_scattering_polynomial()
        return asp_re.kk_transform_inv(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            max_iter=max_iter,
            tolerance=tolerance,
        )

    def calculate_complex_polynomial(
        self,
        target_energies: npt.NDArray[np.floating] | npt.ArrayLike | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
        **kwargs,
    ) -> "asp_complex":
        """
        Create a complex atomic scattering polynomial representation.

        Inverse transforms (`kkcalc.kk_transforms.KK_PP_inv`) the imaginary part of the atomic scattering
        factors to real factors, and then uses both to form a complex polynomial representation.

        Parameters
        ----------
        target_energies : npt.NDArray | npt.ArrayLike | None, optional
            The energies at which to calculate the spectrum.
            By default None, uses the object energies.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.
        **kwargs
            Additional keyword arguments for the `asp_complex` and `atomic_scattering` classes.

        Returns
        -------
        asp_complex
            An atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp_complex : Complex atomic scattering polynomial class.
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Return asp object
        re_asp = self.to_atomic_scattering_polynomial()
        return re_asp.calculate_complex_polynomial(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
            **kwargs,
        )

    def calculate_complex_factors(
        self,
        target_energies: npt.NDArray[np.floating] | npt.ArrayLike | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT_NO_STOICH],
    ) -> "asf_complex":
        """
        Create a complex atomic scattering factor representation.

        Inverse transforms (`kkcalc.kk_transforms.KK_PP_inv`) the real part of the atomic scattering
        factors to imaginary factors, and then uses both to form a complex representation.

        Parameters
        ----------
        target_energies : npt.NDArray | npt.ArrayLike | None, optional
            The energies at which to calculate the spectrum.
            By default None, uses the object energies.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
            By default True.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.
        **kwargs : Unpack[PROPERTIES_DICT_NO_STOICH]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asf_complex
            A complex atomic scattering factor object.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        re: asp_re = self.to_atomic_scattering_polynomial()
        im: asf_im = re.kk_transform_inv(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        re_extended: asp_re = re.extend_energies(im.energies)
        # Return asf object
        common_kwargs = {}
        common_kwargs.update(self._properties_dict)
        common_kwargs.update(kwargs)
        return asf_complex(
            re=re_extended.to_asf(), im=im, **common_kwargs
        )  # asf_complex already pulls properties from components.


class asf_im(asf):
    """
    Identical to `asf`, but reserved for imaginary component factors.

    Parameters
    ----------
    energies : np.ndarray
        Energies in eV.
    factors : np.ndarray
        If data y-data is instead betas or NEXAFS, use the respective `from_<name>` method.
    origin_dtype : KK_Datatype, optional
        Original data type of the atomic scattering factors.
        If not provided, the original data is assumed to be the same as the input data.
    origin_data : np.ndarray, optional
        Original data of the atomic scattering factors.
    scale_to_database : bool, optional
        Whether to scale the data to the background from the stoichiometry.
        By default False.
    **kwargs : Unpack[PROPERTIES_DICT]
        Keyword arguments for the `kkcalc2.models.atomic_scattering` base class.
    """

    def __init__(
        self,
        energies: npt.NDArray[np.floating],
        factors: npt.NDArray[np.floating],
        origin_dtype: KK_Datatype | None = None,
        origin_data: npt.NDArray[np.floating] | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Initialise the atomic scattering base class
        asf.__init__(
            self,
            energies=energies,
            factors=factors,
            origin_dtype=origin_dtype,
            origin_data=origin_data,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @classmethod
    def from_asf(
        cls: type[Self], asf: asf, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> "asf_im":
        """
        Convert an `asf` object to an `asf_im` object.

        Parameters
        ----------
        asf : asf
            Atomic scattering factors object.
        **kwargs : Unpack[PROPERTIES_DICT]
            Keyword arguments for the `atomic_scattering` base class.

        Returns
        -------
        asf_im
            The imaginary atomic scattering factors representation.
        """
        common_kwargs = asf._properties_dict
        common_kwargs.update(kwargs)
        return cls(
            energies=asf.energies,
            factors=asf.factors,
            origin_dtype=asf.origin_dtype,
            origin_data=asf.origin_data,
            **common_kwargs,
        )

    def to_atomic_scattering_polynomial(
        self, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> asp_im:
        """
        Convert the factors representation to an atomic scattering polynomial representation.

        Uses the `energies` and `factors` attributes (with length `N`) of the object
        to create an atomic scattering polynomial object of coefficients with length (N-1).

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asp
            Atomic scattering polynomial object.
        asf_im.atomic_scattering_polynomial
            Calculate the polynomial coefficients for the representation.
        """
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        return asp_im(
            energies=self.energies,
            coefs=self.atomic_scattering_polynomial,
            **common_kwargs,
        )

    # @doc_copy(to_atomic_scattering_polynomial)
    # def to_ASP(self) -> asp_im:
    #     """
    #     Alias for `to_atomic_scattering_polynomial`.
    #     """
    #     return self.to_atomic_scattering_polynomial()
    to_ASP = (
        to_atomic_scattering_polynomial  # Alias for to_atomic_scattering_polynomial
    )

    @property
    def NEXAFS(self) -> np.ndarray:
        r"""
        Convert atomic scattering factors to NEXAFS representation.

        This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
        $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

        .. math::
            \mu_a = 2 r_e \lambda f_2

        Returns
        -------
        np.ndarray
            NEXAFS photoabsorption values corresponding to the `energies` property.

        See Also
        --------
        kkcalc2.models.conversions.ASF_to_NEXAFS : Converts atomic scattering factors to NEXAFS/XANES/Photoabsorption data.
        """
        # TODO: Add documentation about what the NEXAFS scaling is...
        return conversions.ASF_to_NEXAFS(self.energies, self.factors)

    def to_NEXAFS(self) -> tuple[np.ndarray, np.ndarray]:
        r"""
        A tuple of energies and NEXAFS photoabsorption values.

        This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
        $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

        .. math::
            \mu_a = 2 r_e \lambda f_2

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Tuple of energies (eV) and NEXAFS photoabsorption values.

        See Also
        --------
        kkcalc2.models.conversions.ASF_to_NEXAFS : Converts atomic scattering factors to NEXAFS/XANES/Photoabsorption data.
        """
        return self.energies, self.NEXAFS

    @classmethod
    def from_NEXAFS(
        cls: type[Self],
        energies: npt.NDArray[np.floating],
        NEXAFS: npt.NDArray[np.floating],
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert NEXAFS data to imaginary absorption atomic scattering factors (ASF).

        This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
        $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

        .. math::
            \mu_a = 2 r_e \lambda f_2

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        NEXAFS : array_like
            NEXAFS/XANES/photoabsorption data.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_im
            Atomic scattering factors object.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        kkcalc2.models.conversions.ASF_to_NEXAFS : Converts atomic scattering factors to NEXAFS/XANES/Photoabsorption data.
        """
        # TODO: update docs about what the factors in the conversion are about...
        return cls(
            energies=energies,
            factors=conversions.NEXAFS_to_ASF(energies, NEXAFS),
            origin_dtype=KK_Datatype.NEXAFS,
            origin_data=np.c_[energies, NEXAFS],
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @property
    def betas(self) -> npt.NDArray[np.floating]:
        """
        Calculate the imaginary refractive absorption values ($\beta$) from the atomic scattering factors.

        This is equivalent to the `refractive` property.

        Returns
        -------
        npt.NDArray[np.floating]
            The imaginary part of the refractive index.
        """
        # Use the refractive calculation.
        return self.refractive

    @classmethod
    def from_betas(
        cls: type[Self],
        energies: npt.NDArray[np.floating],
        absorption: npt.NDArray[np.floating],
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs,
    ) -> "asf_im":
        r"""
        Convert imaginary refractive absorption values ($\beta$) to atomic scattering factors (ASF).

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        absorption : array_like
            Imaginary index of refraction values (i.e. $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_im
            Imaginary atomic scattering factors object.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        return asf_im.from_refractive(
            energies=energies,
            refractive=absorption,
            number_density=number_density,
            density=density,
            formula_mass=formula_mass,
            stoichiometry=stoichiometry,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    @override
    @classmethod
    def from_refractive_index(
        cls: type[Self],
        energies: npt.NDArray[np.floating],
        refractive_index: npt.NDArray[np.floating],
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert imaginary refractive index values (n -> i * $\beta$) to atomic scattering factors (ASF).

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms|units|molecules per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `density` in grams per millilitre (cm^3), and `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive_index : array_like
            Real/imaginary index of refraction values (i.e. $\delta$'s or $\beta$'s).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_im
            Imaginary atomic scattering factors object.
        """
        # Equivalent to the `from_betas` or `from_refractive` method.
        return asf_im.from_refractive(
            energies=energies,
            refractive=refractive_index,
            number_density=number_density,
            density=density,
            formula_mass=formula_mass,
            stoichiometry=stoichiometry,
            scale_to_database=scale_to_database,
            **kwargs,
        )

    def kk_transform(
        self,
        target_energies: npt.NDArray[np.floating] | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
    ) -> "asf_re":
        """
        Kramers-Kronig transform for the imaginary part of the atomic scattering factors.

        Converts `asf_im` to `asp_im` and uses the `KK_PP` method to calculate the transform.

        Parameters
        ----------
        target_energies : npt.NDArray, optional
            The energies at which to calculate the transform.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.

        Returns
        -------
        asf_re
            Atomic scattering factors object with real components.
        """
        asp_im = self.to_atomic_scattering_polynomial()
        return asp_im.kk_transform(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )

    def calculate_complex_polynomial(
        self,
        target_energies: npt.NDArray[np.floating] | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asp_complex":
        """
        Perform a KK transform and create a complex atomic scattering polynomial representation.

        Converts complex factors to a polynomial representation, and then uses the Kramers-Kronig transform
        to create the real factors, also converted to a polynomial representation to form the complex representation.

        Parameters
        ----------
        target_energies : npt.NDArray, optional
            The energies at which to calculate the transform.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `asp_complex` and `atomic_scattering` classes.

        Returns
        -------
        asp_complex
            An atomic scattering polynomial object.
        """
        im_asp = self.to_atomic_scattering_polynomial()
        return im_asp.calculate_complex_polynomial(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
            **kwargs,
        )

    def calculate_complex_factors(
        self,
        target_energies: npt.NDArray[np.floating] | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = DEF_TOL,
        max_iter: int = DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asf_complex":
        """
        Perform a KK transform and create a complex atomic scattering factors representation.

        Converts the real part of the atomic scattering factors to imaginary factors through
        the Kramers-Kronig transform (`kkcalc.kk_transforms.KK_PP`), and then uses both to form a complex representation.

        Parameters
        ----------
        target_energies : npt.NDArray, optional
            The energies at which to calculate the transform.
        improve_accuracy : bool, optional
            Whether to add extra data points to increase resolution.
        stoichiometry : stoichiometry, optional
            The stoichiometry of the material.
        relativistic_correction : float, optional
            The relativistic correction to apply.
        tolerance : float, optional
            The tolerance for the inverse transform. By default `kkcalc.kk_transforms.DEF_TOL`.
        max_iter : int, optional
            The maximum number of iterations for the inverse transform.
            By default `kkcalc.kk_transforms.DEF_ITER`.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `asf_complex` and `atomic_scattering` classes.

        Returns
        -------
        asf_complex
            A complex atomic scattering factor object.
        """
        im: asp_im = self.to_atomic_scattering_polynomial()
        re: asf_re = im.kk_transform(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        if im.energies.shape != re.energies.shape or np.any(im.energies != re.energies):
            # Extend the energies of the imaginary part to match the real part
            im_extended = im.extend_energies(re.energies)
        else:
            # No extension required
            im_extended = im
        # Convert to complex object
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        return asf_complex(
            re=re, im=im_extended.to_atomic_scattering_factors(), **common_kwargs
        )


class asf_complex(asf_abstract, atomic_scattering):
    """
    Container for a pair of atomic scattering factors, reflecting the real and imaginary parts.

    Parameters
    ----------
    re : asf_re | asf
        Real part of the atomic scattering factors.
    im : asf_im | asf
        Imaginary part of the atomic scattering factors.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments for the `atomic_scattering` subclass.
        Does not effect the real or imaginary parts.

    See Also
    --------
    kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial.
    """

    def __init__(
        self, re: asf_re | asf, im: asf_im | asf, **kwargs: Unpack[PROPERTIES_DICT]
    ):  # numpydoc ignore=GL08
        if not np.all(re.energies == im.energies):
            raise ValueError(
                "Real and imaginary parts must have the same energy intervals."
            )
        if not isinstance(re, asf) or not isinstance(im, asf):
            raise TypeError(f"Real and imaginary parts must be subclasses of {asf}.")

        # Use the real then imaginary part properties to update None values
        common_kwargs = re._properties_dict

        # Check properties are the same
        for key in im._properties_dict:
            if key not in common_kwargs or common_kwargs[key] is None:
                common_kwargs[key] = im._properties_dict[key]
            elif common_kwargs[key] != im._properties_dict[key]:
                warnings.warn(
                    f"Property {key} is different between real {re._properties_dict[key]}"
                    + f" and imaginary parts {im._properties_dict[key]} for {self}."
                )
            else:
                # Ignore if the properties are the same
                pass

        # Update properties with kwargs
        common_kwargs.update(kwargs)

        # Convert to appropriate instance objects
        if isinstance(re, asf):
            re = asf_re.from_asf(re)
        if isinstance(im, asf):
            im = asf_im.from_asf(im)

        # Store attributes
        self._re: asf_re = re
        self._im: asf_im = im

        # Initialise atomic scattering object
        atomic_scattering.__init__(self, **common_kwargs)

    # Override the attomic scattering properties to propogate to the real & imaginary parts
    @atomic_scattering.density.setter
    def density(self, density: float | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asf_complex, self.__class__).density.fset(self, density)
        # Propogate to components
        self._re.density = density
        self._im.density = density

    @atomic_scattering.number_density.setter
    def number_density(
        self, number_density: float | None
    ) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asf_complex, self.__class__).number_density.fset(self, number_density)
        # Propogate to components
        self._re.number_density = number_density
        self._im.number_density = number_density

    @atomic_scattering.formula_mass.setter
    def formula_mass(self, formula_mass: float | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asf_complex, self.__class__).formula_mass.fset(self, formula_mass)
        # Propogate to components
        self._re.formula_mass = formula_mass
        self._im.formula_mass = formula_mass

    @atomic_scattering.stoichiometry.setter
    def stoichiometry(
        self, stoich: kk_stoichiometry | str | None
    ) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction from atomic_scattering
        super(asf_complex, self.__class__).stoichiometry.fset(self, stoich)
        # Propogate to components
        self._re.stoichiometry = stoich
        self._im.stoichiometry = stoich

    @atomic_scattering.name.setter
    def name(self, name: str | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asf_complex, self.__class__).name.fset(self, name)
        # Propogate to components
        self._re.name = name
        self._im.name = name

    @asf_abstract.energies.getter
    def energies(self) -> npt.NDArray[np.floating]:  # numpydoc ignore=PR02
        """
        The property for energies corresponding to the atomic scattering factors.

        Parameters
        ----------
        energies : npt.NDArray | npt.ArrayLike
            The photon energies in eV.

        Returns
        -------
        npt.NDArray
            The photon energies in eV.
        """
        return self._re.energies

    @energies.setter
    def energies(
        self, energies: npt.NDArray[np.floating] | npt.ArrayLike
    ) -> None:  # numpydoc ignore=GL08
        energies = np.asarray(energies)
        self._re.energies = energies
        self._im.energies = energies

    @asf_abstract.factors.getter
    def factors(self) -> npt.NDArray[np.complexfloating]:  # numpydoc ignore=PR02
        """
        The atomic scattering factors property.

        Parameters
        ----------
        factors : tuple[npt.NDArray, npt.NDArray] | npt.NDArray[np.complexfloating]
            The atomic scattering factors, either as a tuple of real and imaginary parts,
            or as a complex array. The first index corresponds to energy.
            A second index can correspond to real & imaginary component.

        Returns
        -------
        np.ndarray[np.complexfloating]
            The atomic scattering factors as a complex array.
        """
        return self._re.factors + 1j * self._im.factors

    @factors.setter
    def factors(
        self,
        factors: (
            tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]
            | npt.NDArray[np.complexfloating]
            | npt.ArrayLike
        ),
    ) -> None:  # numpydoc ignore=GL08
        if (
            isinstance(factors, tuple)
            and isinstance(factors[0], np.ndarray)
            and isinstance(factors[1], np.ndarray)
        ):
            self._re.factors = factors[0]
            self._im.factors = factors[1]
        else:
            factors = np.asarray(factors)
            if isinstance(factors, np.ndarray) and factors.dtype == np.complexfloating:
                self._re.factors = factors.real
                self._im.factors = factors.imag
            elif isinstance(factors, np.ndarray) and factors.ndim == 2:
                warnings.warn(
                    "Ambiguous factors provided. Assuming 2D array with real and imaginary parts."
                )
                self._re.factors = factors[:, 0]
                self._im.factors = factors[:, 1]
            else:
                raise ValueError(
                    "Factors must be a tuple of arrays or a complex ndarray."
                )

    @property
    def abs(self) -> npt.NDArray[np.floating]:
        """
        Absolute values of the atomic scattering factors.

        Returns
        -------
        np.ndarray
            Absolute values of the atomic scattering factors.
        """
        return np.abs(self.factors)

    @property
    def phase(self) -> npt.NDArray[np.floating]:
        """
        Phase of the atomic scattering factors.

        Returns
        -------
        np.ndarray
            Phase of the atomic scattering factors.
        """
        return np.angle(self.factors)

    @property
    def re(self) -> "asf_re":
        """
        Real part object of the atomic scattering factors.

        Returns
        -------
        asf_re
            Real part of the atomic scattering factors.
        """
        return self._re

    @property
    def im(self) -> "asf_im":
        """
        Imaginary part object of the atomic scattering factors.

        Returns
        -------
        asf_im
            Imaginary part of the atomic scattering factors.
        """
        return self._im

    @property
    def refractive(self) -> npt.NDArray[np.complexfloating]:
        r"""
        Refractive coefficients (delta, beta) of the atomic scattering factors.

        Uses `conversions.ASF_to_refractive` to calculate the $\delta$ & $\beta$ values after matching energies to segments.
        .. math::
            = \delta + i\beta

        Returns
        -------
        np.ndarray
            The $\delta$ + i*$\beta$ values from the atomic scattering factors.
        """
        return self._re.refractive + 1j * self._im.refractive

    @property
    def refractive_indexes(self) -> npt.NDArray[np.complexfloating]:
        r"""
        Refractive coefficients (delta, beta) of the atomic scattering factors.

        .. math::
            n = 1 - \delta + i\beta

        Returns
        -------
        np.ndarray
            Refractive index of the atomic scattering factors.
        """
        return self.refractive

    @property
    def deltas(self) -> npt.NDArray[np.floating]:
        r"""
        Real part of the refractive index.

        .. math::
            \delta = 1 - n

        Returns
        -------
        npt.NDArray[np.floating]
            The real part of the refractive index.
        """
        return self._re.deltas

    @property
    def betas(self) -> npt.NDArray[np.floating]:
        r"""
        Imaginary part of the refractive index.

        .. math::
            \beta = i\beta

        Returns
        -------
        npt.NDArray[np.floating]
            The imaginary part of the refractive index.
        """
        return self._im.betas

    def to_atomic_scattering_polynomial(
        self, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> "asp_complex":
        """
        Convert factor to polynomial representation.

        Uses the `energies` and `factors` attributes (with length `N`) of the object
        to create an atomic scattering polynomial object of coefficients with length (N-1).

        For an array of polynomial coefficients, use the `atomic_scattering_polynomial` property.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `asp_complex` and `atomic_scattering` object.

        Returns
        -------
        asp
            Atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp_complex : Complex atomic scattering polynomial object.
        """
        # Get the existing properties, update with kwargs
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        # Return asp object
        return asp_complex(
            re=self._re.to_atomic_scattering_polynomial(),
            im=self._im.to_atomic_scattering_polynomial(),
            **common_kwargs,
        )

    # @doc_copy(to_atomic_scattering_polynomial)
    # def to_ASP(self, **kwargs) -> asp_complex:
    #     """
    #     Alias for `to_atomic_scattering_polynomial`.
    #     """
    #     return self.to_atomic_scattering_polynomial(**kwargs)
    to_ASP = to_atomic_scattering_polynomial

    def contrast(
        self, other: "asf_complex"
    ) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]:
        r"""
        The energy-dependent contrast between two sets of atomic scattering factors.

        The energy-dependent contrast magnitude is calculated by the difference
        in imaginary and real parts squared.
        If the objects have different energy domains, only the common energies are considered,
        and a warning is issued.

        .. math::
            contrast ~ \Delta(\delta)^2 + \Delta\beta^2

        Parameters
        ----------
        other : asf_complex
            The other atomic scattering factors object to compare against.

        Returns
        -------
        energies : np.ndarray
            Array of energy values defined for each contrast value.
        contrast : np.ndarray
            The contrast between two atomic scattering polynomials.
        """
        if self.can_calc_refractive and other.can_calc_refractive:
            if (self.energies.shape == other.energies.shape) and np.all(
                self.energies == other.energies
            ):
                re_diff = self.re.refractive - other.re.refractive
                im_diff = self.im.refractive - other.im.refractive
                return self.energies, np.abs(re_diff) ** 2 + np.abs(im_diff) ** 2
            else:
                energy_subset = np.intersect1d(self.energies, other.energies)
                warnings.warn(
                    f"Energy domains do not match. Only common energies {len(energy_subset)} are considered."
                )
                self_ind = np.searchsorted(self.energies, energy_subset)
                other_ind = np.searchsorted(other.energies, energy_subset)
                re_diff = self.re.refractive[self_ind] - other.re.refractive[other_ind]
                im_diff = self.im.refractive[self_ind] - other.im.refractive[other_ind]
                return energy_subset, np.abs(re_diff) ** 2 + np.abs(im_diff) ** 2
        else:
            raise ValueError(
                "Both objects must have refractive values to calculate contrast, but there is not enough (density, stoichiometry) information."
            )

    @property
    def NEXAFS(self) -> np.ndarray:
        r"""
        Convert atomic scattering factors to NEXAFS representation.

        This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
        $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

        .. math::
            \mu_a = 2 r_e \lambda f_2

        Returns
        -------
        np.ndarray
            NEXAFS photoabsorption values corresponding to the `energies` property.

        See Also
        --------
        kkcalc2.models.conversions.ASF_to_NEXAFS : Converts atomic scattering factors to NEXAFS/XANES/Photoabsorption data.
        """
        return self.im.NEXAFS

    def to_NEXAFS(self) -> tuple[np.ndarray, np.ndarray]:
        r"""
        A tuple of energies and NEXAFS photoabsorption values.

        This convesion treats NEXAFS as equivalent to the `atomic photoabsorption cross section`
        $\mu_a$, as defined by Henke (https://henke.lbl.gov/optical_constants/intro.html):

        .. math::
            \mu_a = 2 r_e \lambda f_2

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Tuple of energies (eV) and NEXAFS photoabsorption values.

        See Also
        --------
        kkcalc2.models.conversions.ASF_to_NEXAFS : Converts atomic scattering factors to NEXAFS/XANES/Photoabsorption data.
        """
        return self.energies, self.NEXAFS

    @classmethod
    def from_NEXAFS(
        cls: type[Self],
        energies: npt.NDArray[np.floating],
        NEXAFS: npt.NDArray[np.floating],
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        NEXAFS photoabsorption data to a full complex representation of scattering factors.

        Only requires the imaginary absorption measurement $\beta$, and will use the transform
        to find the equivalent dispersive $\delta$.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        NEXAFS : array_like
            NEXAFS/XANES/photoabsorption data. Infers real and imaginary parts.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_complex
            Atomic scattering factors object.

        Raises
        ------
        ValueError
            If the NEXAFS data is not complex.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        energies = np.asarray(energies)
        if NEXAFS.dtype != np.floating:
            # NEXAFS = NEXAFS.astype(np.complexfloating)
            # warnings.warn(
            #     "NEXAFS data was not complex. Assuming data is real-component only."
            # )
            raise ValueError(f"NEXAFS data must be a set of float. Was {NEXAFS.dtype}.")

        im = asf_im.from_NEXAFS(
            energies=energies,
            NEXAFS=NEXAFS.imag,
            scale_to_database=scale_to_database,
            **kwargs,
        )
        re = im.kk_transform()

        return cls(re=re, im=im, **kwargs)

    @classmethod
    def from_refractive(
        cls: type[Self],
        energies: npt.NDArray,
        refractive: npt.NDArray[np.complexfloating],
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert refractive components ($\delta$ + i*$\beta$) to atomic scattering factors (ASF).

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-f_1 + i f_2\right)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-(f^0 + f^') + i f^{''}\right)

        Assumes `refractive` is a composition of $\delta(E) + i\beta(E)$, not the refractive index.

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive : array_like
            Infers real and imaginary parts from a complex number ($r= \delta + i\beta$).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_complex
            Atomic scattering factors object.

        Raises
        ------
        ValueError
            If the refractive data is not complex.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Convert energy and beta data to numpy arrays.
        energies = np.asarray(energies)

        if refractive.dtype != np.complexfloating:
            # refractive = refractive.astype(np.complexfloating)
            # warnings.warn(
            #     "Beta data was not complex. Assuming data is real-component only."
            # )
            raise ValueError(
                "Refractive data must be complex. Alternatively use `asf_im.from_refractive` or `asf_re_from_refractive` instead."
            )

        # Accumulate keyword arguments
        common_kwargs = {}
        common_kwargs.update(kwargs)
        common_kwargs.update(
            dict(
                number_density=number_density,
                density=density,
                formula_mass=formula_mass,
                stoichiometry=stoichiometry,
                origin_dtype=KK_Datatype.REFRACTIVE_INDEX,
                origin_data=np.c_[energies, refractive],
            )
        )
        # Return asf instances
        re = asf_re.from_refractive(
            energies=energies,
            refractive=refractive.real,
            scale_to_database=scale_to_database,
            **common_kwargs,
        )
        im = asf_im.from_refractive(
            energies=energies,
            refractive=refractive.imag,
            scale_to_database=scale_to_database,
            **common_kwargs,
        )
        # Create complex class.
        return cls(re=re, im=im, **common_kwargs)

    @classmethod
    def from_refractive_index(
        cls: type[Self],
        energies: npt.NDArray,
        refractive_index: npt.NDArray[np.complexfloating],
        *,
        number_density: float | None = None,
        density: float | None = None,
        formula_mass: float | None = None,
        stoichiometry: kk_stoichiometry | str | None = None,
        scale_to_database: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        r"""
        Convert refractive_index (n= 1 - $\delta$ +i $\beta$) to atomic scattering factors (ASF).

        .. math::
            n(E) &= 1 - \delta(E) + i\beta(E)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-f_1 + i f_2\right)
                 &= 1 + \frac{n_a r_e \lambda^2}{2\pi}\left(-(f^0 + f^') + i f^{''}\right)

        Requires some form of material density information to convert to ASF.
        As per positional argument order, the function will use the first available density information.
        This can either be:
        - `number_density` in atoms per millilitre (cm^3),
        - `density` in grams per millilitre (cm^3), and `formula_mass` (molecular mass),
        - `stoichiometry` as a list of elemental symbol, number pairs or string of a formula.

        Parameters
        ----------
        energies : array_like
            Photon energies in eV.
        refractive_index : array_like
            Infers real and imaginary parts of the index of refraction ($n=1-\delta+i\beta$).
        number_density : float, optional
            Material density in atoms per millilitre (cm^3).
        density : float
            Material density in grams per millilitre (cm^3).
        formula_mass : float
            Atomic mass sum of the materials chemical formula (molecular mass).
            Equivalent to providing a `stoichiometry`.
        stoichiometry : stoichiometry | str
            Description of the combination of elements composing the material.
        scale_to_database : bool, optional
            Whether to scale the atomic scattering factors to the database scale.
            Requires a stoichiometry and a designated complexity (i.e. asf_im or asf_re).
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` object.

        Returns
        -------
        asf_complex
            Complex representation of the atomic scattering factors.

        Raises
        ------
        ValueError
            If the refractive data is not complex.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : Common attributes between atomic scattering factor and polynomial models.
        """
        # Convert energy and beta data to numpy arrays.
        energies = np.asarray(energies)

        if refractive_index.dtype != np.complexfloating:
            # refractive = refractive.astype(np.complexfloating)
            # warnings.warn(
            #     "Beta data was not complex. Assuming data is real-component only."
            # )
            raise ValueError(
                "Refractive index data must be complex. Alternatively use `asf_im.from_refractive` or `asf_re_from_refractive` instead."
            )

        # Accumulate keyword arguments
        common_kwargs = {}
        common_kwargs.update(kwargs)
        common_kwargs.update(
            dict(
                number_density=number_density,
                density=density,
                formula_mass=formula_mass,
                stoichiometry=stoichiometry,
                origin_dtype=KK_Datatype.REFRACTIVE_INDEX,
                origin_data=np.c_[energies, refractive_index],
            )
        )
        # Convert to refractive index to refractive component
        delta = 1 - refractive_index.real
        beta = refractive_index.imag

        # Return asf instances
        re = asf_re.from_refractive(
            energies=energies,
            refractive=delta,
            scale_to_database=scale_to_database,
            **common_kwargs,
        )
        im = asf_im.from_refractive(
            energies=energies,
            refractive=beta,
            scale_to_database=scale_to_database,
            **common_kwargs,
        )
        # Create complex class.
        return cls(re=re, im=im, **common_kwargs)

    @classmethod
    def from_asf(
        cls: type[Self],
        energies: npt.NDArray,
        factors: npt.NDArray[np.complexfloating],
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> Self:
        """
        Convert complex atomic scattering factors to a complex object.

        Parameters
        ----------
        energies : npt.NDArray
            Photon energies in eV.
        factors : npt.NDArray[np.complexfloating]
            Complex atomic scattering factors.
        **kwargs : Unpack[PROPERTIES_DICT]
            Keyword arguments for the atomic scattering base class.

        Returns
        -------
        asf_complex
            An instance of the complex atomic scattering factors class.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : The base class for material attributes.
        """
        re = asf_re.from_asf(asf(energies, factors.real, **kwargs))
        im = asf_im.from_asf(asf(energies, factors.imag, **kwargs))
        return cls(re=re, im=im, **kwargs)

    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> "asf_complex":
        """
        Generate a copy of the `asf` object.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT], optional
            Keyword arguments for the atomic scattering subclass.

        Returns
        -------
        type[asf_complex]
            A new `asf_complex` object with the same polynomial coefficients,
            and properties, but unique memory allocation.

        See Also
        --------
        kkcalc2.models.common.atomic_scattering : The base class for material attributes.
        """
        # Copy the object properties
        common_kwargs = self._properties_dict
        for key in common_kwargs:
            if hasattr(common_kwargs[key], "copy"):
                common_kwargs[key] = common_kwargs[key].copy()
        # Update with kwargs
        common_kwargs.update(kwargs)
        # Create a new object
        return self.__class__(re=self.re.copy(), im=self.im.copy(), **common_kwargs)
