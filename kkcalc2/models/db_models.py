"""
A module for database models.

Allows the calculation of atomic scattering factor data, generated for a given stoichiometry.
"""

import numpy as np
import numpy.typing as npt
import scipy.optimize as opt
from typing import Self, override, overload, Unpack

import abc
from kkcalc2.stoich import stoichiometry as kk_stoichiometry

# Import from submodules of models, as models.py will also call these classes.
from kkcalc2.models.polynomials import asp, asp_im, asp_re, asp_complex
from kkcalc2.models.factors import asf, asf_im, asf_re, asf_complex
from kkcalc2 import conversions
from kkcalc2.models.common import (
    PROPERTIES_DICT,
)

# Load the real/imag scattering factors as they vary with energy
from kkcalc2.asf_database import ASF_DATABASE


class asp_db_abstract(asp, metaclass=abc.ABCMeta):
    """
    Abstract class to define the interface for the atomic scattering polynomial object with database data.

    Requires the stoichiometry to be defined, to obtain scattering factor data from the database.
    Also requires optional arguments for energies and coefs, so that the object can be copied.

    Parameters
    ----------
    stoichiometry : stoichiometry | str
        The stoichiometry of the compound, i.e. the elemental composition.
        Forms the recipe for summation of the database scattering factor data if `energies` and `coefs` are not provided.
    energies : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        A 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the
        coefficients for the atomic scattering factors.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `atomic_scattering` base class.
    """

    @abc.abstractmethod
    def __init__(
        self,
        stoichiometry: kk_stoichiometry | str,
        *,
        energies: npt.NDArray[np.float64] | None = None,
        coefs: npt.NDArray[np.float64] | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        pass

    @override
    def __getitem__(self, key: int | slice) -> Self:
        """
        Copy the object, and apply the key to the internal data.

        Parameters
        ----------
        key : int | slice
            The index or slice to apply to the internal data.

        Returns
        -------
        asp_db_abstract
            A copy of the object with the key applied to the internal data.
        """
        # Copy the object
        props = self._properties_dict  # includes the stoichiometry
        new_obj = self.copy(**props)

        # Convert int index to slice
        if isinstance(key, int):
            key = slice(key, key + 1)
        # Get indexes
        start, stop, step = key.indices(len(self))
        # Apply the key to the internal data
        new_obj._energies = new_obj._energies[start : stop + 1 : step]
        new_obj._coefs = new_obj._coefs[start:stop:step]
        return new_obj

    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Create a copy of the current object.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` base class to pass to the copy function.

        Returns
        -------
        asp_db
            A copy of the current object.
        """
        # Copy the object properties
        props = self._properties_dict  # includes the stoichiometry
        for key in props:
            if hasattr(props[key], "copy"):
                props[key] = props[key].copy()
        props.update(kwargs)  # Update with the new kwargs
        if "stoichiometry" not in props or props["stoichiometry"] is None:
            raise ValueError("Stoichiometry must be defined to copy the object")
        # Copy the object
        return self.__class__(
            energies=self.energies.copy(), coefs=self.coefs.copy(), **props
        )

    @classmethod
    def scale_data(
        cls: type["asp_db_abstract"],
        data_e: npt.ArrayLike,
        data_y: npt.ArrayLike,
        stoichiometry: kk_stoichiometry | str,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
    ) -> npt.NDArray:
        r"""
        Scale the user data to the database data.

        .. math::
            \text{scale} = \frac{\text{Range}(db\_merge\_range)}{\text{Range}(data\_merge\_range)}

        Parameters
        ----------
        data_e : npt.ArrayLike
            The energy values of the user data.
        data_y : npt.ArrayLike
            The atomic scattering factor values of the user data.
        stoichiometry : stoichiometry | str
            The stoichiometry of the compound, i.e. the elemental composition.
        merge_domain : tuple[float, float] | None
            The intersection energies (inclusive bounds) of `data_e` and database data.
            If None, the full range of the user data will be used to scale.
        fix_distortions : bool
            Fits a gradient correction to the ASF data, to minimize offsets between the database
            data and the ASF data. Provides the same functionality as used in
            `asp_db_extended.extend_data_with_db`. Correction is calculated using the merge domain,
            but is applied to the full data set.

        Returns
        -------
        factors: numpy.ndarray
            Y-data scaled to the database data (at the merge domain boundaries).

        See Also
        --------
        kkcalc2.models.db_models.asp_db_extended.extend_data_with_db : Method to extend user data with
            database data, where `merge_domain` truncates the data.
        """
        data_e = np.asarray(data_e)
        data_y = np.asarray(data_y)

        # Get stoichiometry
        if isinstance(stoichiometry, str):
            stoichiometry = kk_stoichiometry(stoichiometry)

        if cls is asp_db_abstract:
            raise TypeError(
                "Cannot call `scale_data` on `asp_db_abstract` directly, use `asp_db_re` or `asp_db_im` classes "
                + "depending on whether you want to scale to real or imaginary components respectively."
            )
        db_poly: type[asp_db_abstract] = cls(stoichiometry)

        db_e: npt.NDArray = db_poly.energies
        assert len(db_e.shape) == 1, "Database energies must be 1D"
        db_coefs: npt.NDArray = db_poly.coefs
        assert len(db_coefs.shape) == 2, "Database coefs must be 2D"

        # Check if merge points are defined:
        if merge_domain is None:
            merge_domain = data_e[[0, -1]]  # full range of the data_asf energies
            data_merge_lb_idx = 0
            data_merge_ub_idx = -1
        else:
            if merge_domain[0] >= merge_domain[1]:
                raise ValueError("Merge domain must be in increasing order")
            # Find the indices and values of the data_asf energies that are within the range of the db_asp energies
            data_merge_lb_idx: int = np.argmax(data_e >= merge_domain[0])
            """
            First (lower bound) index of data within the merge domain
            """

            data_merge_ub_idx: int = np.argmax(data_e > merge_domain[1]) - 1
            """Last (upper bound) index of data within of the merge domain"""

            if data_merge_lb_idx == data_merge_ub_idx:
                raise ValueError(
                    f"Data within domain {merge_domain} must contain more than one energy"
                )

        # Use linear interpolation to find corresponding values of the merge domain.
        data_merge_range = np.interp(merge_domain, data_e, data_y)
        """The range of the data_asf energies within the merge domain"""

        # Find the indices of the spans where the db_asp energies are within the range of the data_asf energies.
        first_domain_idx = np.argmax(db_e > merge_domain[0])
        """First index of db_asp energies within the merge domain"""

        db_asp_merge_lb_idx = (
            first_domain_idx - 1 if first_domain_idx > 0 else 0
        )  # Find value before merge/data edge
        """Last index of lower-bound db_asp energies outside the merge domain"""

        db_asp_merge_ub_idx = np.argmax(
            db_e > merge_domain[1]
        )  # Find value after merge/data edge
        """First index of upper-bound db_asp energies outside the merge domain"""

        # Check if the db merge ub is 0 (i.e. merge_domain[1] is always greater than the db_e)
        if db_asp_merge_ub_idx <= db_asp_merge_lb_idx:
            raise ValueError(
                f"Merge domain ({merge_domain[0]},{merge_domain[1]}) must be within the"
                + f"database energy range ({db_e.min()}, {db_e.max()})"
            )

        db_merge_ub_end = data_merge_ub_idx + 1 if data_merge_ub_idx + 1 != 0 else None
        """The upper bound index (exclusive) for slicing the data_y array."""

        # Calculate the corresponding y values using the polynomial coefs
        db_asp_merge_range = tuple(
            asp.eval_asf_on_coefs(
                target_energies=merge_domain, energies=db_e, coefs=db_coefs
            ).tolist()
        )

        ### Calculate the scale difference between the data_asf and db_asp
        # Range(db_asp) / Range(data_asf)
        scale = (db_asp_merge_range[1] - db_asp_merge_range[0]) / (
            data_merge_range[1] - data_merge_range[0]
        )
        scaled_data_y = (data_y - data_merge_range[0]) * scale + db_asp_merge_range[0]

        # Energy values within the merge domain to use for fitting
        energies = data_e[
            data_merge_lb_idx:db_merge_ub_end
        ]  # essential to only use domain data to perform fit.

        if fix_distortions:
            # Perform a fit along the domain
            db_y = asp.eval_asf_on_coefs(
                target_energies=energies,
                energies=db_e,
                coefs=db_coefs,
            )  # Find equivalent values of the db_asp energies to the data energies
            guess_grad = (
                -(data_merge_range[1] - data_merge_range[0])
                / (db_asp_merge_range[1] - db_asp_merge_range[0])
                / data_y[-1]
            )
            fit_func = asp_db_extended.grad_min
            fit_y = scaled_data_y[
                data_merge_lb_idx:db_merge_ub_end
            ]  # essential to only use domain data to perform fit.

            (grad,), _ = opt.leastsq(
                func=fit_func,
                x0=guess_grad,
                args=(energies, fit_y, db_asp_merge_range, db_y),
                # Minimizes difference to the db_y values.
            )
            # Reassign the scaled data, but use full energy domain to generate values.
            scaled_data_y = db_asp_merge_range[0] + asp_db_extended.grad_min(
                grad,
                data_e,
                scaled_data_y,
                db_asp_merge_range,
                db_y=0,  # No minimization here.
                # Adjust merge indexes, so we can scale whole dataset.
                idx0=data_merge_lb_idx,
                idx1=data_merge_ub_idx,
            )

        return scaled_data_y


class asp_db_re(asp_db_abstract, asp_re):
    """
    Uses stochiometry to calculate a real-component piecewise polynomial representation from Henke, Briggs and Lighthill data.

    Generates a summation of scattering factor data given the chemical stoichiometry, then converts to polynomials.

    Parameters
    ----------
    stoichiometry : stoichiometry | str
        The stoichiometry of the compound, i.e. the elemental composition.
        Forms the recipe for summation of the database scattering factor data if `energies` and `coefs` are not provided.
    energies : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the
        shape of the imaginary spectrum in that segment.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `asp_re` and `atomic_scattering` parent classes.

    Attributes
    ----------
    energies : numpy.ndarray
        A 1D numpy array with length `N` listing the photon energies corresponding to `factors`.
    factors : numpy.ndarray
        A 1D numpy array with length `N` listing the real component of the scattering factor at the corresponding energy.

    See Also
    --------
    asf_database : The atomic scattering factor module for KK calc, where data is sourced from Briggs and Lighthill, and Henke et al.
    """

    def __init__(
        self,
        stoichiometry: kk_stoichiometry,
        *,
        energies: npt.NDArray[np.float64] | None = None,
        coefs: npt.NDArray[np.float64] | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        # Run init
        asp_db_abstract.__init__(
            self, stoichiometry, energies=energies, coefs=coefs, **kwargs
        )

        # Get composition
        if isinstance(stoichiometry, str):
            stoichiometry = kk_stoichiometry(stoichiometry)
            kwargs["stoichiometry"] = (
                stoichiometry  # Ensure stoichiometry is set in the kwargs
            )
        comp = stoichiometry.composition

        # Use for creating copies of the object...
        if energies is not None and coefs is not None:
            if len(energies) - 1 == len(coefs):
                # Ensure the stoichiometry is set in the kwargs
                kwargs["stoichiometry"] = stoichiometry
                asp_re.__init__(self, energies=energies, coefs=coefs, **kwargs)
            else:
                raise ValueError(
                    f"Number of energies ({len(energies)}) and coefs ({len(coefs)}) must match."
                )

        else:
            # Get unique energy points for all elements, but limited to values for the real components energies.
            energies = np.unique(
                np.r_[
                    *[
                        ASF_DATABASE[z]["E"][: ASF_DATABASE[z]["Re"].shape[0]]
                        for z, _ in comp
                    ]
                ]
            )

            # Add weighted asf data sets for KK calculation
            re_factors = np.zeros(
                (len(energies))
            )  # Stores summations of real factors at each energy

            # Stores the current energy index for each element, defining factors at intermediate energies.
            counters = np.zeros(len(comp), dtype=int)

            # Iterate over the unique energies
            for i, energy in enumerate(
                energies
            ):  # iterate over the energies of the ASF_DATABASE
                sum_re = 0  # Sum of the real factors at the current energy
                # Sum the real factors at each energy
                for j, (z, n) in enumerate(comp):
                    # Imaginary coefs at current energy
                    db_re = ASF_DATABASE[z]["Re"][
                        counters[j]
                    ]  # the factor at the current energy
                    sum_re += n * db_re  # Multiply by stoichiometry n
                    # Check if the next energy matches the currently used elemental energy, i.e. end of the valid interval.
                    if ASF_DATABASE[z]["E"][counters[j] + 1] == energy:
                        counters[j] += (
                            1  # Increment counter[j] by 1 if the energy matches, to move to the next energy window
                        )
                # Store the sum of the elemental factors at the current energy
                re_factors[i] = sum_re

            # Convert factors to a polynomial
            coefs = conversions.ASF_to_ASP(energies, re_factors)

            # Setup properties
            kwargs["stoichiometry"] = stoichiometry  # Also store the stoichiometry
            kwargs["is_extended"] = True  # We have extended the data
            asp_re.__init__(self, energies=energies, coefs=coefs, **kwargs)


class asp_db_im(asp_db_abstract, asp_im):
    """
    Uses stochiometry to calculate an imaginary-component piecewise polynomial representation from Henke, Briggs and Lighthill data.

    Summation of scattering factor data given the chemical stoichiometry.

    Parameters
    ----------
    stoichiometry : stoichiometry | str
        The stoichiometry of the compound, i.e. the elemental composition.
        Forms the recipe for summation of the database scattering factor data if `energies` and `coefs` are not provided.
    energies : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the
        shape of the imaginary spectrum in that segment.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `asp_im` and `atomic_scattering` parent classes.

    Attributes
    ----------
    energies : numpy.ndarray
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray
        An 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the shape of the imaginary spectrum in that segment.

    See Also
    --------
    asf_database : The atomic scattering factor module for KK calc, where data is sourced from Briggs and Lighthill, and Henke et al.
    kkcalc2.models.polynomials.asp_im : The atomic scattering polynomial object for the imaginary component of the scattering factor.
    kkcalc2.models.common.atomic_scattering : Base class for atomic scattering objects.
    """

    def __init__(
        self,
        stoichiometry: kk_stoichiometry | str,
        *,
        energies: npt.NDArray[np.float64] | None = None,
        coefs: npt.NDArray[np.float64] | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        # Run init
        asp_db_abstract.__init__(
            self, stoichiometry, energies=energies, coefs=coefs, **kwargs
        )

        # Get composition
        if isinstance(stoichiometry, str):
            stoichiometry = kk_stoichiometry(stoichiometry)
        comp = stoichiometry.composition

        if energies is not None and coefs is not None:
            if len(energies) - 1 == len(coefs):
                # Ensure stoichiometry is set in the kwargs
                kwargs["stoichiometry"] = stoichiometry
                asp_im.__init__(self, energies=energies, coefs=coefs, **kwargs)
            else:
                raise ValueError(
                    f"Number of energies ({len(energies)}) and coefs ({len(coefs)}) must match."
                )
        else:
            # Get unique energy points for all elements
            energies = np.unique(np.r_[*[ASF_DATABASE[z]["E"] for z, _ in comp]])

            # Add weighted asf data sets for KK calculation
            im_coefs = np.zeros(
                (len(energies) - 1, 5)
            )  # Stores summations of imaginary coefficients at each energy

            # Stores the current energy index for each element, defining coefficients at intermediate energies.
            counters = np.zeros(len(comp), dtype=int)
            # Iterate over the unique energies
            for i, energy in enumerate(
                energies[1:]
            ):  # iterate over the energies of the ASF_DATABASE, coefs run between N-1 and Nth energy
                sum_im = 0  # Sum of the imaginary coefficients at the current energy
                # Sum the imaginary coefficients at each energy
                for j, (z, n) in enumerate(comp):
                    # Imaginary coefs at current energy
                    db_im_coefs = ASF_DATABASE[z]["Im"][
                        counters[j]
                    ]  # the imaginary piecewise polynomial coefficients
                    sum_im += n * db_im_coefs  # Multiply by stoichiometry n
                    # Check if the next energy matches the currently used elemental energy, i.e. end of the valid interval.
                    if ASF_DATABASE[z]["E"][counters[j] + 1] == energy:
                        counters[j] += (
                            1  # Increment counter[j] by 1 if the energy matches, to move to the next energy window
                        )
                # Store the sum of the imaginary coefficients at the current energy
                im_coefs[i, :] = sum_im

            # Setup properties
            kwargs["stoichiometry"] = stoichiometry  # Also store the stoichiometry
            kwargs["is_extended"] = True  # We have extended the data
            asp_im.__init__(self, energies=energies, coefs=im_coefs, **kwargs)


class asp_db_complex(asp_complex):
    """
    Uses stochiometry to calculate a complex-component piecewise polynomial representation from Henke, Briggs and Lighthill data.

    Summation of scattering factor data given the chemical stoichiometry.

    Parameters
    ----------
    stoichiometry : stoichiometry | str
        The stoichiometry of the compound, i.e. the elemental composition.
        Forms the recipe for summation of the database scattering factor data if `energies` and `coefs` are not provided.
    energies : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray | None, optional
        Keyword argument to enable copy method, by default None.
        An 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the
        shape of the imaginary spectrum in that segment.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments for the `atomic_scattering` base class.

    Attributes
    ----------
    energies : numpy.ndarray
        An N+1 length array listing the starting photon energies of the segments that the spectrum is broken up into.
    coefs : numpy.ndarray
        An 2D numpy array with dimensions (N, 5) in which each row lists the polynomial coefficients describing the shape of the imaginary spectrum in that segment.

    See Also
    --------
    asf_database : The atomic scattering factor module for KK calc, where data is sourced from Briggs and Lighthill, and Henke et al.
    kkcalc2.models.polynomials.asp_complex : The atomic scattering polynomial object for the complex component of the scattering factor.
    kkcalc2.models.common.atomic_scattering : Base class for atomic scattering objects.
    """

    def __init__(
        self,
        stoichiometry: kk_stoichiometry | str,
        *,
        energies: npt.NDArray[np.float64] | None = None,
        coefs: npt.NDArray[np.complex128] | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        # Get composition
        if isinstance(stoichiometry, str):
            stoichiometry = kk_stoichiometry(stoichiometry)

        if energies is not None and coefs is not None:
            if len(energies) - 1 == len(coefs):
                re_db = asp_db_re(
                    stoichiometry, energies=energies, coefs=coefs.real, **kwargs
                )
                im_db = asp_db_im(
                    stoichiometry, energies=energies, coefs=coefs.imag, **kwargs
                )
                asp_complex.__init__(self, re_db, im_db, **kwargs)
            else:
                raise ValueError(
                    f"Number of energies -1 ({len(energies) - 1}) and coefs ({len(coefs)}) must match."
                )
        else:
            # Use asp_re and asp_im to generate the complex component
            kwargs["stoichiometry"] = stoichiometry  # Also store the stoichiometry
            re_db = asp_db_re(**kwargs)
            im_db = asp_db_im(**kwargs)

            # Setup properties
            kwargs["is_extended"] = True  # We have extended the data
            asp_complex.__init__(self, re=re_db, im=im_db, **kwargs)

    @override
    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Create a copy of the database object by copying the energies and coefficients.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` base class to pass to the copy function.
            These can be used to override object properties such as stoichiometry, etc.
            The underlying data will not be modified by these kwargs however.

        Returns
        -------
        Self
            A copy of the current object with the same energies and coefficients.
        """
        # Create a new object with the same energies and coefficients
        energies = self.energies.copy()
        coefs = self.coefs.copy()
        # Ensure the stoichiometry is set in the kwargs
        prop_args = self._properties_dict  # includes the stoichiometry
        for key, value in kwargs.items():
            if value is not None:
                prop_args[key] = value.copy() if hasattr(value, "copy") else value  # type: ignore
        return self.__class__(energies=energies, coefs=coefs, **prop_args)  # type: ignore

    @classmethod
    def scale_data(
        cls: type["asp_db_complex"],
        data_e: npt.ArrayLike,
        data_y: npt.ArrayLike,
        stoichiometry: kk_stoichiometry | str,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
        """
        Scale the user data to the database data.

        Parameters
        ----------
        data_e : npt.ArrayLike
            The energy values of the user data.
        data_y : npt.ArrayLike
            The atomic scattering factor values of the user data.
        stoichiometry : stoichiometry | str
            The stoichiometry of the compound, i.e. the elemental composition.
        merge_domain : tuple[float, float] | None
            The intersection energies of the user data and database data.
            If None, the full range of the user data will be used.
        fix_distortions : bool
            Flag to fix distortions in the user data. Provides the same functionality as used in
            `asp_db_extended.extend_data_with_db`.

        Returns
        -------
        energies: numpy.ndarray
            The truncated energies within the merge domain.
        factors: numpy.ndarray
            Atomic scattering factors scaled to the database data (at the merge domain boundaries).
        """
        # Separate the data into real and imaginary components
        data_e = np.asarray(data_e)
        data_y = np.asarray(data_y, dtype=np.complex128)
        data_re = data_y.real
        data_im = data_y.imag
        # Use the db to scale the data
        energies, data_re = asp_db_re.scale_data(
            data_e, data_re, stoichiometry, merge_domain, fix_distortions
        )
        energies2, data_im = asp_db_im.scale_data(
            data_e, data_im, stoichiometry, merge_domain, fix_distortions
        )
        assert np.all(energies == energies2), (
            "Energies for real and imaginary components do not match after scaling."
        )
        # Combine the data back into a complex array
        data_y = data_re + 1j * data_im
        # Return the scaled data
        return energies, data_y


class asp_db_extended(asp):
    """
    Class for extending an `asp` object with database scattering factor data.

    Merges scattering factor polynomials with the user-provided near-edge data.
    TODO: Implement asp input for dataset, to preserve coefficients.

    Parameters
    ----------
    data_asf : asf | list[asf]
        The atomic scattering factor object.
    database : asp_db | kk_stoichiometry | str
        The atomic scattering potential object, generated for a given material stoichiometry.
        Can also be a `kk_stoichiometry` object or a string representing the stoichiometry,
        which will be converted to an `asp_db` object.
    merge_domain : tuple[float, float] | list[tuple[float, float]] | None
        The range of energies to merge the user data_asf with the db_asp data.
    fix_distortions : bool
        Flag to fix distortions in the user data_asf.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `atomic_scattering` base classes.

    Attributes
    ----------
    dataset_asf : asf
        The original `asf` (atomic scattering factor) object containing the user data,
        used to generate the extended `asp` object.
    database : asp_db
        The original `asp_db` (atomic scattering polynomial) object containing the database data,
        used to extend the data contained in the `asf` object.

    See Also
    --------
    asf_database : The atomic scattering factor module for KK calc, where data is sourced from Briggs and Lighthill, and Henke et al.
    asp_db : The atomic scattering polynomial object for the imaginary component of the scattering factor for a given stoichiometry.
    kkcalc2.models.polynomials.asp_im : The atomic scattering polynomial object for the imaginary component of the scattering factor.
    kkcalc2.models.common.atomic_scattering : Base class for atomic scattering objects.
    """

    @overload
    def __init__(
        self,
        data_asf: asf,
        database: asp_db_abstract,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None: ...  # numpydoc ignore=GL08

    @overload
    def __init__(
        self,
        data_asf: list[asf],
        database: asp_db_abstract,
        merge_domain: list[tuple[float, float] | None] | None = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None: ...  # numpydoc ignore=GL08

    def __init__(
        self,
        data_asf: asf | list[asf],
        database: asp_db_abstract,
        merge_domain: (
            tuple[float, float] | list[tuple[float, float] | None] | None
        ) = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Check if data_asf is a list
        if isinstance(data_asf, list):
            # Check if merge_domain is a list
            if isinstance(merge_domain, list):
                # Check if the lengths match
                if len(data_asf) != len(merge_domain):
                    raise ValueError("Length of data_asf and merge_domain must match")

                # Check each merge domain doesn't overlap
                for i, md1 in enumerate(merge_domain[:-1]):
                    # Collect the first domain
                    if md1 is None:
                        md1_lb, md1_ub = data_asf[i].energies[[0, -1]]
                    else:
                        md1_lb, md1_ub = md1
                    for j, md2 in enumerate(merge_domain[i + 1 :]):
                        # Collect the second domain
                        if md2 is None:
                            md2_lb, md2_ub = data_asf[j + i + 1].energies[[0, -1]]
                        else:
                            md2_lb, md2_ub = md2
                        if (md1_lb > md2_lb and md1_lb < md2_ub) or (
                            md1_ub > md2_lb and md1_ub < md2_ub
                        ):
                            raise ValueError(
                                f"Merge domains must not overlap. #{i} ({md1}) and #{j + i + 1} ({md2}) overlap.)"
                            )

            elif merge_domain is None:
                # Check the data boundaries don't overlap
                for i, d1 in enumerate(data_asf[:-1]):
                    for j, d2 in enumerate(data_asf[i + 1 :]):
                        md1_lb, md1_ub = d1.energies[[0, -1]]
                        md2_lb, md2_ub = d2.energies[[0, -1]]
                        if (md1_lb > md2_lb and md1_lb < md2_ub) or (
                            md1_ub > md2_lb and md1_ub < md2_ub
                        ):
                            raise ValueError(
                                f"ASF data energy domains must not overlap. #{i} ({d1}) and #{j + i + 1} ({d2}) overlap.)"
                            )

            else:
                # Raise an error if merge_domain is not a list
                raise ValueError("data_asf is a list, merge_domain must be a list.")

        else:
            # Check if merge_domain is a list
            if isinstance(merge_domain, list):
                if len(merge_domain) < 2:
                    # Reduce the merge_domain for copying purposes. (see asp_db_extended.copy())
                    if len(merge_domain) == 1:
                        merge_domain = merge_domain[0]
                    else:
                        merge_domain = None
                else:
                    # Raise an error if merge_domain is a list but not the data_asf.
                    raise ValueError(
                        "data_asf is not a list, merge_domain must not be a list."
                    )
                    # Check if the merge_domain is a single tuple / None
            else:
                if merge_domain is not None and not isinstance(merge_domain, tuple):
                    raise ValueError(
                        "`merge_domain` must be a tuple or None if data_asf is not a list."
                    )

            # Put the data_asf / merge_domain into a list to iterate over later
            data_asf = [data_asf]  # Convert to list for iteration
            merge_domain = [merge_domain]  # Convert to list for iteration

        # Check sorted
        for d in data_asf:
            if not np.all(np.diff(d.energies) > 0):
                raise ValueError(f"Data energies for {d} must be sorted")

        # Store construction parameters
        self.merge_domain = merge_domain
        """The merge domain(s) used to choose which `data_asf` values are used to create the extended data."""
        self.fix_distortions = fix_distortions
        """The fix distortions flag used to add extra processing to the provided `data_asf`."""

        # Get the data pointers from the asp object
        asp_e, asp_coefs = database.energies, database.coefs
        if asp_coefs is None:
            raise ValueError("Database coefs must be defined to extend the data")

        # Assign the database data as the default merge data before looping.
        merge_e, merge_coefs = asp_e, asp_coefs
        for i, d_asf in enumerate(data_asf):
            # Merge domain
            if merge_domain is not None and isinstance(merge_domain, list):
                d_merge = merge_domain[i]
            elif merge_domain is None:
                d_merge = merge_domain
            else:
                raise ValueError("Merge domain must be a list of tuples or None")

            ### 1. Alignment of Energy Values:
            # Get the data pointers from the asf object
            data_e, data_y = d_asf.data
            # Extend the data with the merge/database data
            merge_e, merge_coefs = self.extend_data_with_db(
                data_e=data_e,
                data_y=data_y,
                db_e=merge_e,
                db_coefs=merge_coefs,
                merge_domain=d_merge,
                fix_distortions=fix_distortions,
            )

        # Copy the kwargs from the data_asf / asp objects if not None.
        # Data first, for naming priority, then database.
        # Database stoichiometry is more important to reflect extension operation made on user data.
        extra_kwargs = data_asf[0]._properties_dict
        for d_asf in data_asf[1:]:
            # Add the properties from the data_asf object to the kwargs
            for key in d_asf._properties_dict:
                if key not in extra_kwargs:
                    extra_kwargs[key] = d_asf._properties_dict[key]

        # Replace `None` values with the db_asp values
        for key in extra_kwargs:
            if extra_kwargs[key] is None or key == "stoichiometry":
                # Update using the db values.
                extra_kwargs[key] = database._properties_dict[key]

        # Add to the kwargs if not already present, otherwise kwargs takes precedence.
        for key in extra_kwargs.keys():
            if key not in kwargs:
                kwargs[key] = extra_kwargs[key]

        # Update the kwargs to reflect the extended data
        kwargs["is_extended"] = True  # We have extended the data

        # Initialize the asp_im object
        super().__init__(energies=merge_e, coefs=merge_coefs, **kwargs)

        # Store the data_asf object for reference
        self.dataset_asf: asf | list[asf] = (
            data_asf if len(data_asf) > 1 else data_asf[0]
        )
        """
        The original `asf` (atomic scattering factor) object containing the user data,
        used to generate the extended `asp` object.
        """

        self.database_asp: asp_db_abstract = database
        """
        The original `asp` (atomic scattering polynomial) object containing the database data,
        used to extend the `asf` object.
        """
        return

    @staticmethod
    def extend_data_with_db(
        data_e: npt.NDArray,
        data_y: npt.NDArray,
        db_e: npt.NDArray,
        db_coefs: npt.NDArray,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
    ) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Merge the user data (factors) with the database data (polynomial coefs).

        Parameters
        ----------
        data_e : npt.NDArray
            The energy values of the user data.
        data_y : npt.NDArray
            The atomic scattering factor values of the user data.
        db_e : npt.NDArray
            The energy values of the database data.
        db_coefs : npt.NDArray
            The atomic scattering factor polynomial coefficients of the database data.
        merge_domain : tuple[float, float] | None, optional
            The intersection energies (inclusive bounds) of `data_e` and database data.
            If None, the full range of the user data will be used to scale.
        fix_distortions : bool, optional
            Fits a gradient correction to the ASF data, to minimize offsets between the database
            data and the ASF data. Correction is calculated using the merge domain,
            but is applied to the full data set.

        Returns
        -------
        npt.NDArray
            The merged energy values.
        npt.NDArray
            The merged atomic scattering factor polynomial coefficients.

        Raises
        ------
        ValueError
            Merge domain must be in increasing order.
        ValueError
            Data within the provided energy domain must contain more than one datapoint.
        """

        # Check if merge points are defined:
        if merge_domain is None:
            merge_domain = data_e[[0, -1]]  # full range of the data_asf energies
            data_merge_lb_idx = 0
            data_merge_ub_idx = -1
        else:
            if merge_domain[0] >= merge_domain[1]:
                raise ValueError("Merge domain must be in increasing order")
            # Find the indices and values of the data_asf energies that are within the range of the db_asp energies
            data_merge_lb_idx: int = np.argmax(data_e >= merge_domain[0])
            """
            First (lower bound) index of data within (inclusive) the merge domain
            """

            data_merge_ub_idx: int = np.argmax(data_e > merge_domain[1]) - 1
            """Last (upper bound) index of data within (inclusive) of the merge domain"""

            if data_merge_lb_idx == data_merge_ub_idx:
                raise ValueError(
                    f"Data within domain {merge_domain} must contain more than one energy"
                )

        # Use linear interpolation to find corresponding values of the merge domain.
        data_merge_range = np.interp(merge_domain, data_e, data_y)
        """The range of the data_asf energies within the merge domain"""

        # Find the indices of the spans where the db_asp energies are within the range of the data_asf energies.
        first_domain_idx = np.argmax(db_e > merge_domain[0])
        """First index of db_asp energies within the merge domain"""

        db_asp_merge_lb_idx = (
            first_domain_idx - 1 if first_domain_idx > 0 else 0
        )  # Find value before merge/data edge
        """Last index of lower-bound db_asp energies outside the merge domain"""

        db_asp_merge_ub_idx = np.argmax(
            db_e > merge_domain[1]
        )  # Find value after merge/data edge
        """First index of upper-bound db_asp energies outside the merge domain"""

        # Check if the db merge ub is 0 (i.e. merge_domain[1] is always greater than the db_e)
        if db_asp_merge_ub_idx <= db_asp_merge_lb_idx:
            raise ValueError(
                f"Merge domain ({merge_domain[0]},{merge_domain[1]}) must be within the"
                + f"database energy range ({db_e.min()}, {db_e.max()})"
            )

        db_merge_ub_end = data_merge_ub_idx + 1 if data_merge_ub_idx + 1 != 0 else None
        """The upper bound index (exclusive) for slicing the data_y array."""

        # Calculate the corresponding y values using the polynomial coefs
        db_asp_merge_range = tuple(
            asp.eval_asf_on_coefs(
                target_energies=merge_domain, energies=db_e, coefs=db_coefs
            ).tolist()
        )

        ### Calculate the scale difference between the data_asf and db_asp
        # Range(db_asp) / Range(data_asf)
        scale = (db_asp_merge_range[1] - db_asp_merge_range[0]) / (
            data_merge_range[1] - data_merge_range[0]
        )
        scaled_data_y = (data_y - data_merge_range[0]) * scale + db_asp_merge_range[0]

        # The energy values within the merge domain used for fitting
        energies = data_e[
            data_merge_lb_idx:db_merge_ub_end
        ]  # essential to only use domain data to perform fit.

        if fix_distortions:
            # Perform a fit along the domain
            db_y = asp.eval_asf_on_coefs(
                target_energies=energies,
                energies=db_e,
                coefs=db_coefs,
            )  # Find equivalent values of the db_asp energies to the data energies
            guess_grad = (
                -(data_merge_range[1] - data_merge_range[0])
                / (db_asp_merge_range[1] - db_asp_merge_range[0])
                / data_y[-1]
            )
            fit_func = asp_db_extended.grad_min
            fit_y = scaled_data_y[
                data_merge_lb_idx:db_merge_ub_end
            ]  # essential to only use domain data to perform fit.
            (grad,), _ = opt.leastsq(
                func=fit_func,
                x0=guess_grad,
                args=(energies, fit_y, db_asp_merge_range, db_y),
            )
            # Reassign the scaled data
            merge_data_e = energies
            merge_data_y = db_asp_merge_range[0] + asp_db_extended.grad_min(
                grad,
                energies,
                fit_y,
                db_asp_merge_range,
                db_y=0,
                idx0=0,
                idx1=-1,
            )
        else:
            # Construct the merge data to use
            merge_data_e = energies
            merge_data_y = scaled_data_y[data_merge_lb_idx:db_merge_ub_end]

        # Add merge domain to the merge data if not already present
        if merge_domain[0] != merge_data_e[0]:
            merge_data_e = np.r_[merge_domain[0], merge_data_e]
            merge_data_y = np.r_[db_asp_merge_range[0], merge_data_y]
        if merge_domain[1] != merge_data_e[-1]:
            merge_data_e = np.r_[merge_data_e, merge_domain[1]]
            merge_data_y = np.r_[merge_data_y, db_asp_merge_range[1]]

        # Convert factors to coefficients
        merge_data_coefs = conversions.ASF_to_ASP(
            energies=merge_data_e, factors=merge_data_y
        )

        # Add the db sections to the merge data
        merge_e = merge_data_e
        merge_coefs = merge_data_coefs

        # Boundary already added, so finish at idx-1.
        if db_asp_merge_lb_idx > 0:
            merge_e = np.r_[db_e[0 : db_asp_merge_lb_idx + 1], merge_e]
            merge_coefs = np.r_[db_coefs[0 : db_asp_merge_lb_idx + 1], merge_coefs]
        # Boundary already added, so start at idx+1 for energies, and idx for coefs.
        if db_asp_merge_ub_idx < len(db_e):
            merge_e = np.r_[merge_e, db_e[db_asp_merge_ub_idx:]]
            merge_coefs = np.r_[merge_coefs, db_coefs[db_asp_merge_ub_idx - 1 :]]

        return merge_e, merge_coefs

    @staticmethod
    def grad_min(
        grad,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        db_merge_range: tuple[float, float],
        db_y: npt.NDArray[np.float64] | float | int,
        idx0: int = 0,
        idx1: int = -1,
    ) -> npt.NDArray[np.float64]:
        r"""
        Minimizing function to  a general gradient of the data to the database.

        .. math::
            \min\left(\frac{(y - y_{i_0}) - m (x - x_{i_0})}{y_{i_1} - y_{i_0} - m (x_{i_1} - x_{i_0})} \cdot (f_{db, max} - f_{db, min}) - f_{db}\right)

        Parameters
        ----------
        grad : float
            The gradient to fit.
        x : npt.NDArray[np.float64]
            The x data values to fit.
        y : npt.NDArray[np.float64]
            The y data values to fit.
        db_merge_range : tuple[float, float]
            The range of the database atomic scattering factors at the domain endpoints.
            This maps onto the amplitude of the gradient line.
        db_y : npt.NDArray[np.float64] | float | int
            The database atomic scattering factor values (not coefs) to fit.
        idx0 : int, optional
            The starting index where the gradient difference is calculated as zero, by default 0.
        idx1 : int, optional
            The ending index where the gradient difference is normalised as one, by default -1.

        Returns
        -------
        npt.NDArray[np.float64]
            The difference between the fitted data and the database data.
        """
        # Initial gradient from 0 point to a final value.
        data_grad_diff = (y - y[idx0]) - grad * (x - x[idx0])  # 0 to some number
        # Gradient over the whole domain.
        data_grad_diff_total = (y[idx1] - y[idx0]) - grad * (
            x[idx1] - x[idx0]
        )  # some number
        # Normalize the gradient as a function of the final total gradient.
        norm_grad_diff = (
            data_grad_diff / data_grad_diff_total
        )  # Evolves from 0 to 1, from idx0 to idx1.
        # Range of the database values over the domain.
        db_range = db_merge_range[1] - db_merge_range[0]  # Range of the database values
        # Difference between the gradient data scaled to the database range, and the database values.
        return norm_grad_diff * db_range - db_y

    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Create a copy of the current object.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments to pass to `atomic_scattering` base classes.

        Returns
        -------
        asp_db_extended
            A copy of the current object.
        """
        # Copy the object properties
        kwargs = self._properties_dict
        for key in kwargs:
            if hasattr(kwargs[key], "copy"):
                kwargs[key] = kwargs[key].copy()
        # Create the copy
        obj = self.__class__(
            data_asf=self.dataset_asf.copy(),
            database=self.database_asp.copy(),
            merge_domain=self.merge_domain,
            fix_distortions=self.fix_distortions,
            **kwargs,
        )
        # Need to copy the energies and coefs to ensure exact match
        # even though constructor should generate same data.
        # These may have updated with extend/truncate operations.
        obj._energies = self.energies.copy()
        obj._coefs = self.coefs.copy()
        return obj


class asp_db_im_extended(asp_db_extended, asp_im):
    """
    The extended imaginary-component atomic scattering polynomial object.

    Forms an imaginary part extension of atomic scattering factor data, using the database data.

    Parameters
    ----------
    data_asf : asf | asf_im | list[asf | asf_im]
        The atomic scattering factors.
    database : asp_db_im | kk_stoichiometry | str
        The database atomic scattering polynomial, generated for a given material stoichiometry.
        Can also be a `kk_stoichiometry` object or a string representing the stoichiometry,
        which will be converted to an `asp_db` object.
    merge_domain : tuple[float, float] | None
        The intersection energies to merge the user data_asf with the db_asp data (in eV).
        By default None, using full data domain.
    fix_distortions : bool
        Flag to fix distortions in the user data_asf. By default False.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `atomic_scattering` base classes.
    """

    def __init__(
        self,
        data_asf: asf | asf_im | list[asf | asf_im],
        database: asp_db_im | kk_stoichiometry | str,
        merge_domain: tuple[float, float] | list[tuple[float, float]] | None = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Convert the database to an asp_db_im object
        if isinstance(database, str):
            stoichiometry = kk_stoichiometry(database)
            im_db = asp_db_im(stoichiometry)
        elif isinstance(database, kk_stoichiometry):
            im_db = asp_db_im(database)
        elif isinstance(database, asp_db_im):
            im_db = database
        else:
            raise ValueError(
                "Database must be a stoichiometry, string, or asp_db_im object"
            )

        # Construct the extended object
        super().__init__(
            data_asf=data_asf,
            database=im_db,
            merge_domain=merge_domain,
            fix_distortions=fix_distortions,
            **kwargs,
        )


class asp_db_re_extended(asp_db_extended, asp_re):
    """
    The extended real-component atomic scattering polynomial object.

    Forms a real part extension of atomic scattering factor data, using the database data.

    Parameters
    ----------
    data_asf : asf | asf_re | list[asf | asf_re]
        The atomic scattering factors.
    database : asp_db_re | kk_stoichiometry | str
        The database atomic scattering polynomial, generated for a given material stoichiometry.
        Can also be a `kk_stoichiometry` object or a string representing the stoichiometry,
        which will be converted to an `asp_db` object.
    merge_domain : tuple[float, float] | None
        The intersection energies to merge the user data_asf with the db_asp data (in eV).
        By default None, using full data domain.
    fix_distortions : bool
        Flag to fix distortions in the user data_asf. By default False.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `atomic_scattering` base classes.
    """

    def __init__(
        self,
        data_asf: asf | asf_re | list[asf | asf_re],
        database: asp_db_re | kk_stoichiometry | str,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Convert the database to an asp_db_re object
        if isinstance(database, str):
            stoichiometry = kk_stoichiometry(database)
            re_db = asp_db_re(stoichiometry)
        elif isinstance(database, kk_stoichiometry):
            re_db = asp_db_re(database)
        elif isinstance(database, asp_db_re):
            re_db = database
        else:
            raise ValueError(
                "Database must be a stoichiometry, string, or asp_db_im object"
            )

        super().__init__(data_asf, re_db, merge_domain, fix_distortions, **kwargs)


class asp_db_complex_extended(asp_db_extended, asp_complex):
    """
    The extended complex-component atomic scattering polynomial object.

    Forms a complex part extension of atomic scattering factor data, using the database data.

    Parameters
    ----------
    data_asf : asf | asf_complex | list[asf | asf_complex]
        The atomic scattering factors.
    database : asp_db_complex | kk_stoichiometry | str
        The database atomic scattering polynomial, generated for a given material stoichiometry.
        Can also be a `kk_stoichiometry` object or a string representing the stoichiometry,
        which will be converted to an `asp_db` object.
    merge_domain : tuple[float, float] | None
        The intersection energies to merge the user data_asf with the db_asp data (in eV).
        By default None, using full data domain.
    fix_distortions : bool
        Flag to fix distortions in the user data_asf. By default False.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments to pass to `atomic_scattering` base classes.
    """

    def __init__(
        self,
        data_asf: asf_complex | list[asf_complex],
        database: asp_db_complex | kk_stoichiometry | str,
        merge_domain: tuple[float, float] | None = None,
        fix_distortions: bool = False,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> None:  # numpydoc ignore=GL08
        # Convert the database to an asp_db_complex object
        if isinstance(database, str):
            stoichiometry = kk_stoichiometry(database)
            complex_db = asp_db_complex(stoichiometry)
        elif isinstance(database, kk_stoichiometry):
            complex_db = asp_db_complex(database)
        elif isinstance(database, asp_db_complex):
            complex_db = database
        else:
            raise ValueError(
                "Database must be a stoichiometry, string, or asp_db_im object"
            )

        super().__init__(data_asf, complex_db, merge_domain, fix_distortions, **kwargs)


if __name__ == "__main__":
    ## Test various formulas
    # Setup graph
    import matplotlib.pyplot as plt

    plots = plt.subplots(2, 2)
    fig: plt.Figure = plots[0]
    ax: plt.Axes = plots[1][0][0]
    ax2: plt.Axes = plots[1][0][1]
    ax3: plt.Axes = plots[1][1][0]
    ax4: plt.Axes = plots[1][1][1]

    P3MEET = "C9H12O6S2"  # C9H11O3S
    CARBON = "C"
    ANTIMONY = "Sb"
    BISMUTH = "Bi"
    TELLURIUM = "Te"
    SELINIUM = "Se"
    SULFUR = "S"

    # compounds = [P3MEET, CARBON, SULFUR, ANTIMONY, BISMUTH, TELLURIUM, SELINIUM]
    # compounds = [P3MEET, CARBON, SULFUR]
    compounds = [ANTIMONY, BISMUTH, TELLURIUM, SELINIUM]

    for compound in compounds:
        stoich = kk_stoichiometry(compound)
        stoich_asp = stoich.asp_im()

        # Convert all energies to asf:
        energies = stoich_asp.energies
        stoich_asf = stoich_asp.to_asf()

        # Graph the asf
        scat = ax.scatter(
            energies, stoich_asf.factors, s=1, alpha=0.5, label=f"{compound} ASF"
        )
        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel("ASF Data")
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Plot the polynomials
        for i, e1 in enumerate(energies[:-1]):
            e2 = energies[i + 1]
            x = np.linspace(e1, e2, 100)
            x_asf = stoich_asp.eval_asf(target_energies=x)
            ax.plot(
                x,
                x_asf,
                linewidth=0.5,
                c=scat.get_edgecolor(),
                label=f"'{compound}' Polynomial" if i == 0 else None,
            )
    ax.set_title("Atomic Scattering Factors of Elements and Compounds")
    ax.legend()

    # Create a merge of physical data and database data
    POLYSTYRENE = "CH"
    PS_NAME = "Polystyrene"
    ps_stoich = kk_stoichiometry(POLYSTYRENE)
    asp_db_PS = asp_db_im(ps_stoich)

    # Import Data
    import os

    data_dir = os.path.join(os.path.dirname(__file__), "../../examples/data")
    data_file = os.path.normpath(os.path.join(data_dir, "PS_004_-dc.txt"))
    data_PS = np.genfromtxt(data_file, skip_header=4)

    # Convert to KK Calc objects
    from kkcalc2.models import asf_im, asf_re

    assert data_PS.shape[1] == 2, "Data file must have two columns"
    asf_PS = asf_im(
        energies=data_PS[:, 0], factors=data_PS[:, 1], stoichiometry=ps_stoich
    )

    # Combine the data with the database
    asp_db_PS_extended = asp_db_im_extended(
        data_asf=asf_PS,
        database=asp_db_PS,
        merge_domain=(280, 320),
        # fix_distortions=False
    )

    asp_db_PS_extended_fixed = asp_db_im_extended(
        data_asf=asf_PS,
        database=asp_db_PS,
        merge_domain=(280, 320),
        fix_distortions=True,
    )

    extended_asf = asp_db_PS_extended.to_atomic_scattering_factors()
    extended_asf_fixed = asp_db_PS_extended_fixed.to_atomic_scattering_factors()
    ax2.plot(
        extended_asf.energies, extended_asf.factors, label=f"{PS_NAME} Extended ASF"
    )
    ax2.plot(
        extended_asf_fixed.energies,
        extended_asf_fixed.factors,
        label=f"{PS_NAME} Extended ASF Fixed",
    )
    db_asf = asp_db_PS.to_asf()
    ax2.plot(asp_db_PS.energies, db_asf.factors, label=f"{PS_NAME} DB ASF")
    ax2.set_xlim(270, 330)
    # ax2.set_xscale("log")
    # ax2.set_ylim(450, 900)

    ax2.set_title(PS_NAME + " Imaginary Extension")
    ax2.legend()

    ### Axis 3:
    # Perform Real Extension
    asp_db_PS_re = asp_db_re(ps_stoich)
    ax3.plot(asp_db_PS_re.energies, asp_db_PS_re.asf, label=f"{PS_NAME} DB Real ASF")
    # Transform the imag database
    asp_db_transf = asp_db_PS.kk_transform()
    ax3.plot(
        asp_db_transf.energies,
        asp_db_transf.factors,
        label=f"{PS_NAME} kkTransform (Imag Database)",
    )
    # Transform the extended data
    asp_transf = asp_db_PS_extended.kk_transform()
    ax3.plot(
        asp_transf.energies,
        asp_transf.factors,
        label=f"{PS_NAME} kkTransform (Imag Extended Data)",
    )
    # Extend the transformed data
    data_trans = asf_PS.kk_transform(max_iter=5)
    extended_re_data = asp_db_re_extended(data_asf=data_trans, database=asp_db_PS_re)
    ax3.plot(
        extended_re_data.energies,
        extended_re_data.asf,
        label=f"{PS_NAME} Real-Extended (transformed imag data)",
    )

    # Axes:
    ax3.set_title(PS_NAME + " Real Extensions")
    ax3.set_xlabel("Energy [eV]")
    ax3.set_ylabel("ASF Data")
    # ax3.set_xscale("log")
    ax3.set_xlim(270, 330)
    ax3.legend()

    fig.tight_layout()
    plt.show()
