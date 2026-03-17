"""
'Piecewise polynomial' representation models of scattering factors.
"""

import abc
import numpy as np
import numpy.typing as npt
import warnings
from collections.abc import Iterator
from typing import TYPE_CHECKING, Self, overload, Unpack, override
import scipy.constants as sc

from kkcalc2.models.common import (
    atomic_scattering,
    atomic_scattering_abstract,
    PROPERTIES_DICT,
    PROPERTIES_DICT_NO_STOICH,
)
from kkcalc2 import conversions, transforms
from kkcalc2.stoich import stoichiometry as kk_stoichiometry

if TYPE_CHECKING:
    from kkcalc2.models.factors import (
        asf as asf_type,
        asf_im,
        asf_re,
        asf_complex,
        asf_abstract,
    )

has_pandas: bool
"""Flag to check if pandas is available."""
try:
    import pandas as pd

    has_pandas = True
except ImportError:
    has_pandas = False


class asp_abstract(atomic_scattering_abstract, metaclass=abc.ABCMeta):
    """
    Abstract class for a piecewise polynomial representation of atomic scattering factors.

    See Also
    --------
    kkcalc2.models.common.atomic_scattering_abstract : Base interface for atomic scattering.
    """

    @property
    @abc.abstractmethod
    def coefs(self) -> npt.NDArray:
        """
        Abstract property for the polynomial coefficients defining scattering factors between energy intervals.

        Returns
        -------
        npt.NDArray | None
            The polynomial coefficients for the scattering factors, with shape `(N, M)`,
            where N is the number of segments and `M` is the number of polynomial coefficients.
        """
        pass

    @property
    @abc.abstractmethod
    def energies(self) -> npt.NDArray:
        """
        Abstract property for the energy intervals defining the polynomial coefficients.

        Returns
        -------
        npt.NDArray
            The energy values defining the intervals for the polynomial coefficients.
            Has length `N+1`, where `N` is the number of segments.
        """
        pass

    @property
    @abc.abstractmethod
    def orders(self) -> npt.NDArray | None:
        """
        Abstract property for the polynomial orders of the scattering factors.

        Returns
        -------
        npt.NDArray | None
            The polynomial orders for the scattering factors, with length `M`. If None,
            then kkcalc2 internally assumes the polynomial orders are by default [1, 0, -1, -2, -3].
        """
        pass

    @staticmethod
    # @doc_copy(conversions.ASP_to_ASF)
    def coefs_to_atomic_scattering_factors(
        energies: npt.NDArray, coefs: npt.NDArray, orders: npt.NDArray | None = None
    ) -> npt.NDArray:
        r"""
        Alias for static method `conversions.ASP_to_ASF`.

        Calculates the atomic scattering factors from polynomial `coefs` defined between `energies`.

        Parameters
        ----------
        energies : npt.NDArray
            The energy values of length `N+1` defining the `N` intervals for the polynomial coefficients.
        coefs : npt.NDArray
            The polynomial coefficients of shape `(N, M)` for the scattering factors,
            defined on the intervals of `energies` where `M` is the number of coefficients.
        orders : npt.NDArray | None, optional
            The polynomial orders for the scattering factors, with length `M`. If None,
            coefficients are assumed to be in the order [1, 0, -1, -2, -3]. By default None.

        Returns
        -------
        npt.NDArray
            The atomic scattering factors calculated from the polynomial coefficients.
        """
        return conversions.ASP_to_ASF(energies=energies, coefs=coefs, orders=orders)

    @property
    def atomic_scattering_factors(self) -> npt.NDArray:
        """
        Calculate `N+1` atomic scattering factors from the `N` piecewise polynomial coefficients.

        Returns
        -------
        npt.NDArray
            The atomic scattering factors calculated from the polynomial coefficients.
        """
        return self.coefs_to_atomic_scattering_factors(
            energies=self.energies, coefs=self.coefs, orders=self.orders
        )

    # @property
    # @doc_copy(atomic_scattering_factors)
    # def asf(self) -> npt.NDArray:
    #     """
    #     Alias for `atomic_scattering_factors`.
    #     """
    #     return self.atomic_scattering_factors

    asf = atomic_scattering_factors  # Alias for atomic_scattering_factors

    @abc.abstractmethod
    def to_atomic_scattering_factors(
        self,
        target_energies: npt.ArrayLike | npt.NDArray | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asf_abstract":
        """
        Convert the piecewise polynomial representation to an atomic scattering factor object.

        Parameters
        ----------
        target_energies : npt.NDArray | None, optional
            Energy values at which to calculate the atomic scattering factors.
            By default None, which uses the object's energies.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `asf` or `atomic_scattering` classes.

        Returns
        -------
        asf
            An atomic scattering factor object with the same polynomial coefficients as the piecewise polynomial.
        """
        pass

    # @doc_copy(to_atomic_scattering_factors)
    # @abc.abstractmethod
    # def to_asf(self) -> "asf_abstract":
    #     """
    #     Alias for `to_atomic_scattering_factors`.
    #     """
    #     pass
    to_asf = to_atomic_scattering_factors

    @staticmethod
    def eval_asf_on_coefs(
        target_energies: npt.ArrayLike,
        energies: npt.ArrayLike,
        coefs: npt.ArrayLike,
        orders: npt.NDArray | None = None,
    ) -> npt.NDArray:
        """
        Calculate the atomic scattering factors at the `target_energies`.

        Uses the provided polynomial `coefs` defined over `energies` intervals.
        Can also provide the polynomial orders for the coefficients.

        Parameters
        ----------
        target_energies : npt.NDArray
            The energies of length `L` at which to calculate the atomic scattering factors.
        energies : npt.NDArray
            The energy intervals of length `N` defining the polynomial coefficients.
        coefs : npt.NDArray
            The coefficients of shape `(N, M)` the polynomial defining the scattering factors,
            where `M` is the number of polynomial coefficients.
        orders : npt.NDArray | None
            The polynomial orders for the scattering factors, with length `M`. If None,
            coefficients are assumed to be in the order [1, 0, -1, -2, -3]. By default None.

        Returns
        -------
        npt.NDArray
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.

        Raises
        ------
        ValueError
            If the target energies are outside the defined energy domain.
        """
        target_energies = np.asarray(target_energies)
        energies = np.asarray(energies)
        coefs = np.asarray(coefs)
        if orders is not None:
            orders = np.asarray(orders)

        assert np.all(energies[:-1] <= energies[1:]), (
            "Energies must be in increasing order."
        )
        # Find where the energies are located in the object's energies.
        indices = np.asarray(
            np.searchsorted(energies, target_energies) - 1
        )  # subtract to transfer from spans (N+1) to coefficients (N).
        if -1 in indices or len(energies) - 1 in indices:
            # Check if all searchsorted invalid indexes are defined.
            invalid = np.where((indices < 0) | (indices == len(energies) - 1))
            inval_energies = target_energies[invalid]
            for inv_e in inval_energies:
                if inv_e == target_energies[0]:
                    indices[invalid] = 0  # First defined polynomial.
                elif inv_e == target_energies[-1]:
                    indices[invalid] = len(energies) - 2  # Last defined polynomial.
                else:
                    raise ValueError(
                        f"Some energies {target_energies[invalid]} "
                        + f"are outside the defined energy range ({energies.min()}, {energies.max()})."
                    )

        # Collate coefficients corresponding to the energies.
        target_coefs = coefs[indices]
        # Calculate the ASF values at the given energies.
        factors = asp_abstract.coefs_to_atomic_scattering_factors(
            energies=target_energies, coefs=target_coefs, orders=orders
        )
        return factors

    @overload
    def __call__(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray:  # numpydoc ignore=GL08
        pass

    @overload
    def __call__(self, target_energies: float | int) -> float:  # numpydoc ignore=GL08
        pass

    def __call__(
        self, target_energies: npt.NDArray | npt.ArrayLike | float | int | None = None
    ) -> npt.NDArray | float:
        r"""
        Calculate scattering factors from object polynomial coefficients at desired `energies`.

        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | npt.NDArray | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        # Type check
        if target_energies is None:
            # If no energies or coefficients are provided, use the object's values to return the intrinsic ASF values.
            return self.eval_asf(self.energies)

        if not isinstance(target_energies, (int, float)):
            target_energies = np.asarray(target_energies)
            factors = self.eval_asf_on_coefs(
                target_energies=target_energies,
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )
        else:
            target_energies = np.array([target_energies])
            # Remove the singleton dimension from the output.
            factors = self.eval_asf_on_coefs(
                target_energies=target_energies,
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )[0]
        return factors

    eval_asf = __call__  # Alias for __call__

    # # @doc_copy(eval_asf)
    # def __call__(
    #     self, target_energies: npt.NDArray | float | None = None
    # ) -> npt.NDArray | float:
    #     """
    #     Callable alias for `evaluate_energies`.
    #     """
    #     return self.eval_asf(target_energies)

    def __iter__(self) -> Iterator[tuple[tuple[float, float], np.ndarray]]:
        """
        Provide each segment and piecewise polynomial coefficients of the energy-dependent scattering amplitude.

        Yields
        ------
        segment : tuple[float, float]
            The energy interval for which the polynomial coefficients are valid.
        poly_coefs : np.ndarray
            The polynomial coefficients for the scattering factor in the given energy interval.
        """
        for i in range(len(self.energies) - 1):
            yield (self.energies[i], self.energies[i + 1]), self.coefs[i]

    def __getitem__(self, key: int | slice) -> Self:
        """
        A new `asp` object truncated to the specified index.

        Contains the same energy spans and corresponding polynomial coefficients as the original,
        but sliced to the specified key.

        Parameters
        ----------
        key : int | slice
            The index or slice to select the segment of the polynomial coefficients, and corresponding energy interval.

        Returns
        -------
        type[asp_abstract]
            A new `asp` object with the same energy spans and corresponding polynomial coefficients,
            but sliced to the specified index.
        """
        # Collect the kwargs
        kwargs = self._properties_dict

        # Convert int index to slice
        if isinstance(key, int):
            key = slice(key, key + 1)

        # Slice the energies and coefficients
        start, stop, step = key.indices(len(self))
        energies = self.energies[start : stop + 1 : step]
        coefs = self.coefs[start:stop:step]
        return self.__class__(
            energies=energies, coefs=coefs, orders=self.orders, **kwargs
        )

    def dataframe(self) -> "pd.DataFrame":
        """
        Generate a Pandas representation of the coefficients list, useful for display.

        Returns
        -------
        pd.DataFrame
            A dataframe of the energies spanning the coefficients at given orders.

        Raises
        ------
        ImportError
            If pandas is not available but the method is called.
        """
        if not has_pandas:
            raise ImportError("Pandas is required for this method.")
        orders = self.orders
        if orders is not None:
            return pd.DataFrame(
                np.c_[self.energies[:-1], self.energies[1:], *self.coefs.T],
                columns=["Energy LB", "Energy UB", *[f"A{order}" for order in orders]],
            )
        else:
            if self.energies.shape[0] == self.coefs.shape[0] + 1:
                return pd.DataFrame(
                    np.c_[self.energies[:-1], self.energies[1:], *self.coefs.T],
                    columns=["Energy LB", "Energy UB", "A1", "A0", "A-1", "A-2", "A-3"],
                )
            else:
                idx = min(self.energies.shape[0] - 1, self.coefs.shape[0])
                return pd.DataFrame(
                    np.c_[
                        self.energies[:idx],
                        self.energies[1 : idx + 1],
                        *self.coefs.T[:, :idx],
                    ],
                    columns=["Energy LB", "Energy UB", "A1", "A0", "A-1", "A-2", "A-3"],
                )

    to_pandas = dataframe  # Alias for dataframe method.

    def __str__(self) -> str:
        """
        A string representation of polynomial object.

        Returns
        -------
        str
            A string representation of the object properties.
        """
        header1: str = "ASP" + ("" if self.name is None else f" '{self.name}'")
        header2: str = (
            f"{self.coefs.shape[0]} en segments, {self.coefs.shape[1]} coefficients."
        )
        return f"{header1}({header2})"

    def __repr__(self, **kwargs) -> str:
        """
        A string representation of the coefficient list.

        Uses Pandas if available, to create a string representation of the coefficient list.
        Rows displayed are the first and last 5 if more than 10 rows.

        If Pandas is not avaialble, manually

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `pd.dataFrame.to_string` method.
            I.e. `max_rows` is defaulted to 6.

        Returns
        -------
        str
            A string representation of the coefficients.
        """
        if has_pandas:
            header1: str = "ASP" + ("" if self.name is None else f" '{self.name}'")
            # Create a default max_rows if not provided.
            if "max_rows" not in kwargs:
                kwargs["max_rows"] = 6
            return header1 + "\n" + self.dataframe().to_string(**kwargs)
        else:
            # Manually show the first and last 5
            header1: str = "ASP" + ("" if self.name is None else f" '{self.name}'")
            header2: str = "Energy0\tEnergy1\t" + "\t".join(
                [
                    "C_" + str(e)
                    for e in (
                        self.orders.tolist()
                        if self.orders is not None
                        else [1, 0, -1, -2, -3]
                    )
                ]
            )
            header3: str = "".join(["-"] * len(header2))

            data_head: list[str] = []
            data_tail: list[str] = []
            M = self.coefs.shape[1]
            for i in range(0, M):
                tail_line = []
                head_line = []
                # Add energy values.
                head_line.extend(self.energies[i : i + 2].tolist())
                tail_line.extend(
                    self.energies[-6 + i : -4 + i if i - 4 != 0 else None].tolist()
                )
                # Add the coefficients
                head_line.extend(self.coefs[i].tolist())
                tail_line.extend(self.coefs[-5 + i].tolist())
                # Add to the data
                data_head.append("\t".join([f"{val:3f}" for val in head_line]))
                data_tail.append("\t".join([f"{val:3f}" for val in tail_line]))

            return "\n".join([header1, header2, header3, *data_head, "...", *data_tail])

    def __len__(self) -> int:
        """
        The number of segments (N) in the polynomial.

        Not to be confused with the number of points (N+1).

        Returns
        -------
        int
            N - The number of segments in the polynomial representation.
        """
        return self.coefs.shape[0]

    @abc.abstractmethod
    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Generate a copy of the `asp` object.

        Parameters
        ----------
        **kwargs
            Any keyword arguments for the constructors to update the copy properties.

        Returns
        -------
        type[asp_abstract]
            A new `asp` object with the same atomic scattering polynomial coefficients
            and properties, but unique memory allocation.
        """
        pass


class asp(asp_abstract, atomic_scattering):
    """
    Atomic scattering polynomial.

    A generic piecewise polynomial representation of scattering factors.
    Allows the evaluation of the scattering factors at specified energies, by calling
    the object or using the `evaluate_energies` method.

    Parameters
    ----------
    energies : npt.ArrayLike
        The energy values of length `N+1` defining the `N` intervals for the polynomial coefficients.
    coefs : npt.ArrayLike
        The polynomial coefficients of shape `(N, M)` for the scattering factors,
        defined on the intervals of `energies` where `M` is the number of coefficients.
    orders : npt.ArrayLike | None, optional
        The polynomial orders for the scattering factors. If None, then kkcalc2 internally
        assumes the polynomial orders are by default [1, 0, -1, -2, -3]. By default None.
        Must have length `M` if provided.
    **kwargs : Unpack[PROPERTIES_DICT], optional
        Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` such as:
        - `number_density` : float
        - `density` : float
        - `stoich` : stoichiometry
        - `formula_mass` : float
        - `name` : str

    Raises
    ------
    ValueError
        If the energies are not in increasing order, the dimensions of `energies` and `coefs` do not match,
        or the dimensions of `orders` do not match the number of coefficients.

    See Also
    --------
    kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
    """

    def __init__(
        self,
        energies: npt.ArrayLike,
        coefs: npt.ArrayLike,
        orders: npt.ArrayLike | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        # Initialise atomic scattering object
        atomic_scattering.__init__(self, **kwargs)

        # Convert inputs to numpy arrays is not already
        energies = np.asarray(energies)
        coefs = np.asarray(coefs)
        if energies.ndim != 1:
            raise ValueError("Energies must be a 1D array.")
        if coefs.ndim != 2:
            raise ValueError("Coefficients must be a 2D array.")

        # Check energies are monotonic
        diff_sign = np.diff(energies) > 0  # True = Positive, False = Negative.
        if not np.all(diff_sign):
            raise ValueError(
                "Energies must be in increasing order. Indexes of non-monotonic values: ",
                np.where(~diff_sign)[0],
            )

        # Check input dimensions match
        if len(energies) != len(coefs) + 1:
            raise ValueError(
                "Pairs of energies define the intervals for each set of polynomial coefficients. "
                + f"Number of coefficients ({len(coefs)}) does not match the number of energies ({len(energies)} - 1)."
            )

        # Check orders if provided
        if orders is not None:
            orders = np.asarray(orders)
            if orders.ndim != 1:
                raise ValueError("Orders must be a 1D array.")
            if len(orders) != coefs.shape[1]:
                raise ValueError(
                    "Number of orders must match the number of coefficients."
                )

        # Store attributes
        self._energies = energies
        self._coefs = coefs
        self._orders = orders

    @asp_abstract.energies.getter
    def energies(self) -> npt.NDArray:  # numpydoc ignore=PR02
        """
        Attribute for the interval energy values, between which the `coefs` are defined.

        Setting the `energies` will discard the existing `coefs` if they do not match the new length.

        Parameters
        ----------
        energies : npt.NDArray
            The energy values of length N+1 defining the N intervals for the polynomial coefficients.
            If `coefs` exist without shape `(N, M)`, then the coefs are discarded.

        Returns
        -------
        npt.NDArray
            An array of energy values with length N+1, where N is the number of segments.
        """
        return self._energies

    @energies.setter
    def energies(self, energies: npt.NDArray) -> None:  # numpydoc ignore=GL08
        self._energies = energies
        # Wipe coefficients if the energies are changed to a different length.
        if self.coefs is not None and len(energies) != len(self.coefs) + 1:
            warnings.warn(
                f"({self}) Energies have changed length. Coefficients set to `None`."
            )
            self._coefs = None

    def extend_energies(
        self, new_energies: npt.NDArray, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> Self:
        """
        Extend the `energies` of the `asp` object to include `new_energies` values.

        Uses the existing interval to generate a new `asp` object with the same
        polynomial coefficients, but defined on the existing and `new_energies`.
        All `new_energies` must be within the `energies` domain.

        Parameters
        ----------
        new_energies : npt.NDArray
            The new energy values to extend the intervals.
            Existing `energies` must be a subset of `new_energies`.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asp
            A new `asp` object with the same polynomial coefficients, but defined on the new `energies`.

        Raises
        ------
        ValueError
            If `new_energies` are not a subset of existing `energies`.
        """
        # Check energies are a subset
        en_min, en_max = self.energies.min(), self.energies.max()
        if not np.all([en >= en_min and en <= en_max for en in new_energies]):
            raise ValueError("Existing energies must be a subset of the new energies.")
        # Get the class and creation kwargs
        cls = self.__class__
        # props = self._properties_dict
        # update the properties with the kwargs
        # props.update(kwargs)

        # Find new energies not in the existing energies
        nocoef_energies = np.setdiff1d(new_energies, self.energies)
        if len(nocoef_energies) == 0:
            # If no new energies, just return a copy of the object
            # obj = self.copy(**props)
            obj = self.copy(**kwargs)
            return obj

        # Find the indices of the existing energies in the new energies
        indices = np.searchsorted(self.energies, nocoef_energies)
        # Collect coefs: -1 because the coefficients are defined on the previous index interval.
        new_coefs = np.array([self.coefs[i - 1] for i in indices])
        if len(new_coefs.shape) == 1:
            new_coefs = new_coefs.reshape(-1, 1)  # Ensure new_coefs is 2D
        # Create new energy and coefficients
        unsorted_en = np.r_[self.energies[:-1], nocoef_energies]
        sort_indices = np.argsort(
            unsorted_en
        )  # Sort the combination of old (except the last bound) and new energies
        energies = unsorted_en[sort_indices]
        coefs = np.r_[self.coefs, new_coefs][
            sort_indices[:-1]
        ]  # Exclude the last bound

        # Check if class is asp_db, in which case the constructor cannot take energies.
        from kkcalc2.models import asp_db_abstract, asp_db_extended

        if issubclass(cls, (asp_db_abstract, asp_db_extended)):
            obj = self.copy(**kwargs)  # Use kwargs not common_kwargs due to copying.
            # Set internal attributes to avoid length warnings
            assert len(energies) - 1 == len(coefs)
            obj._energies = energies
            obj._coefs = coefs
            return obj
        else:
            common_kwargs = self._properties_dict
            common_kwargs.update(kwargs)
            return cls(
                energies=energies, coefs=coefs, orders=self.orders, **common_kwargs
            )

    def truncate_energies(
        self, domain: tuple[float, float], **kwargs: Unpack[PROPERTIES_DICT]
    ) -> Self:
        """
        Truncate the `energies` of the `asp` object to the `domain` values.

        Uses the existing interval to generate a new `asp` object with the same
        polynomial coefficients, but defined between `domain` values.
        `domain` values must be within the existing `energies` values.

        Parameters
        ----------
        domain : tuple[float, float]
            The domain values to truncate the intervals.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asp
            A new `asp` object with the same polynomial coefficients, but truncated to the `domain`.
        """
        # Re-order the domain values if needed
        if domain[0] >= domain[1]:
            domain = (domain[1], domain[0])
        # Check domain is within the existing energies
        if not np.all(
            [en >= self.energies.min() and en <= self.energies.max() for en in domain]
        ):
            raise ValueError("Domain values must be within the existing energies.")

        # Get the class
        cls = self.__class__

        print(self.__class__.__name__, self.energies.shape)

        if domain[0] in self.energies and domain[1] in self.energies:
            # If the domain is already in the energies, just return a copy of the object
            obj = self.copy(**kwargs)
            # Update the energies and coefficients to the domain
            index = (self.energies >= domain[0]) & (self.energies <= domain[1])
            new_energies = obj.energies[index]
            new_coefs = obj.coefs[index[:-1] & index[1:]]
            # Set internal attributes to avoid length warnings
            assert len(new_energies) - 1 == len(new_coefs)
            obj._energies = new_energies
            obj._coefs = new_coefs
            return obj
        else:
            # Add domain values in the existing energies if not already present:
            nocoef_energies = np.setdiff1d(domain, self.energies)
            # Find the indices of the existing energies in the new energies
            indices = np.searchsorted(self.energies, nocoef_energies)
            # Collect coefs: -1 because the coefficients are defined on the previous index interval.
            new_coefs = np.array([self.coefs[i - 1] for i in indices])
            if len(new_coefs.shape) == 1:
                new_coefs = new_coefs.reshape(-1, 1)  # Ensure new_coefs is 2D
            # Create new energy and coefficients
            sort_indices = np.argsort(
                np.r_[self.energies[:-1], nocoef_energies]
            )  # Sort the combination of old (except the last bound) and new energies
            energies = np.r_[self.energies, nocoef_energies][sort_indices]
            coefs = np.r_[self.coefs, new_coefs][sort_indices]
            # Truncate the energies and coefficients to the domain
            index = (energies >= domain[0]) & (energies <= domain[1])
            energies = energies[index]
            coefs = coefs[index[:-1] & index[1:]]

            # Check if class is asp_db, in which case the constructor cannot take energies.
            from kkcalc2.models import asp_db_abstract, asp_db_extended

            if issubclass(cls, (asp_db_abstract, asp_db_extended)):
                obj = self.copy(**kwargs)
                # Set internal attributes to avoid length warnings
                assert len(energies) - 1 == len(coefs)
                obj._energies = energies
                obj._coefs = coefs
                return obj
            else:
                common_kwargs = self._properties_dict
                common_kwargs.update(kwargs)
                return cls(
                    energies=energies, coefs=coefs, orders=self.orders, **common_kwargs
                )

    @asp_abstract.coefs.getter
    def coefs(self) -> npt.NDArray | None:  # numpydoc ignore=PR02
        """
        Return the polynomial coefficients for the scattering factor, defined on the intervals of `energies`.

        Parameters
        ----------
        coefs : npt.NDArray
            The polynomial coefficients of shape `(N, M)` for the scattering factors.

        Returns
        -------
        npt.NDArray
            A 2D array, where rows correspond to the segments defined by `energies`, and columns are the polynomial coefficients.
        """
        return self._coefs

    @coefs.setter
    def coefs(self, coefs: npt.NDArray) -> None:  # numpydoc ignore=GL08
        if len(self.energies) - 1 != len(coefs):
            raise ValueError(
                f"Number of coefficients ({len(coefs)}) must match the number of energy intervals ({len(self.energies) - 1})."
            )
        self._coefs = coefs

    @asp_abstract.orders.getter
    def orders(self) -> npt.NDArray | None:
        """
        Return the polynomial orders for the scattering factors.

        Returns
        -------
        npt.NDArray | None
            A 1D array of polynomial orders, with length M, where M is the number of coefficients.
            If None, then kkcalc2 internally assumes the polynomial orders are by default [1, 0, -1, -2, -3].
        """
        return self._orders

    @override
    def to_atomic_scattering_factors(
        self,
        target_energies: npt.ArrayLike | npt.NDArray | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asf_type":
        """
        Convert the piecewise polynomial representation to an atomic scattering factor object.

        Parameters
        ----------
        target_energies : npt.NDArray | None, optional
            Energy values at which to calculate the atomic scattering factors.
        **kwargs
            Additional keyword arguments for the `asf` or `atomic_scattering` classes.

        Returns
        -------
        asf
            An atomic scattering factor object with the same polynomial coefficients as the piecewise polynomial.

        See Also
        --------
        kkcalc2.models.factors.asf : Atomic scattering factor object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        from kkcalc2.models.factors import asf as asf_type

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        if target_energies is None:
            return asf_type(
                energies=self.energies,
                factors=self.atomic_scattering_factors,
                **common_kwargs,
            )
        else:
            target_energies = np.asarray(target_energies)
            return asf_type(
                energies=target_energies,
                factors=self.eval_asf(target_energies),
                **common_kwargs,
            )

    # @doc_copy(to_atomic_scattering_factors)
    # def to_asf(self, **kwargs) -> "asf_type":
    #     """
    #     Alias for `to_atomic_scattering_factors`.
    #     """
    #     return self.to_atomic_scattering_factors(**kwargs)

    to_asf = to_atomic_scattering_factors  # Alias for to_atomic_scattering_factors

    @overload
    def eval_refractive(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_refractive(
        self, target_energies: float | int
    ) -> float | complex: ...  # numpydoc ignore=GL08

    def eval_refractive(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray | float | complex:
        r"""
        Determine energy-dependent refractive index component.

        Determined from the object atomic scattering polynomial coefficients at desired `energies`.
        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        if not self.can_calc_refractive:
            raise AttributeError(f"{self} cannot calculate delta/beta values.")
        # Type check
        if target_energies is None:
            # If no energies or coefficients are provided, use the object's values to return the intrinsic ASF values.
            return conversions.ASF_to_refractive(
                self.energies,
                self.eval_asf(self.energies),
                number_density=self.number_density,
                density=self.density,
                formula_mass=self.formula_mass,
                stoichiometry=self.stoichiometry,
            )

        if not isinstance(target_energies, (int, float)):
            target_energies = np.asarray(target_energies)
            factors = self.eval_asf_on_coefs(
                target_energies=target_energies,
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )
        else:
            # Remove the singleton dimension from the output.
            factors = self.eval_asf_on_coefs(
                target_energies=np.array([target_energies]),
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )[0]
        return conversions.ASF_to_refractive(
            target_energies,
            factors,
            number_density=self.number_density,
            density=self.density,
            formula_mass=self.formula_mass,
            stoichiometry=self.stoichiometry,
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
        # Check keys are valid properties
        for key in kwargs:
            if key not in PROPERTIES_DICT.__annotations__.keys():
                raise ValueError(f"Invalid property: {key}.")

        # Copy the object properties
        common_kwargs = self._properties_dict
        for key in common_kwargs:
            if hasattr(common_kwargs[key], "copy"):
                common_kwargs[key] = common_kwargs[key].copy()
        # Update the common kwargs with provided values.
        common_kwargs.update(kwargs)
        # Create a new object
        return self.__class__(
            energies=self.energies.copy(),
            coefs=self.coefs.copy(),
            orders=self.orders.copy() if self.orders is not None else None,
            **common_kwargs,
        )


class asp_im(asp):
    """
    Identical to `asp`, but reserved for the imaginary component of the atomic scattering.

    Enables kk algorithms to convert to real and complex representations of the atomic scattering factors.

    Parameters
    ----------
    energies : npt.ArrayLike
        The energy values of length `N+1` defining the `N` intervals for the polynomial coefficients.
    coefs : npt.ArrayLike
        The polynomial coefficients of shape `(N, M)` for the scattering factors,
        defined on the intervals of `energies` where `M` is the number of coefficients.
    orders : npt.ArrayLike | None, optional
        The polynomial orders for the scattering factors. If None, then kkcalc2 internally
        assumes the polynomial orders are by default [1, 0, -1, -2, -3]. By default None.
        Must have length `M` if provided.
    **kwargs : Unpack[PROPERTIES_DICT], optional
        Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` base class.
    """

    def __init__(
        self,
        energies: npt.ArrayLike,
        coefs: npt.ArrayLike,
        orders: npt.ArrayLike | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        super().__init__(energies, coefs, orders, **kwargs)

    @classmethod
    def from_asp(
        cls: type["asp_im"], asp: asp, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> "asp_im":
        """
        Convert an undesignated `asp` object to a type of `asp_im` object.

        Parameters
        ----------
        asp : asp
            Atomic scattering polynomial object.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asp_im
            An imaginary-part designated atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        common_kwargs = asp._properties_dict
        common_kwargs.update(kwargs)
        return cls(
            energies=asp.energies, coefs=asp.coefs, orders=asp.orders, **common_kwargs
        )

    def to_atomic_scattering_factors(
        self,
        target_energies: (
            npt.NDArray[np.floating | np.integer] | npt.ArrayLike | None
        ) = None,
        **kwargs,
    ) -> "asf_im":
        """
        Convert the piecewise polynomial representation to an atomic scattering factor object.

        Parameters
        ----------
        target_energies : npt.NDArray[np.floating | np.int_] | None, optional
            Energy values at which to calculate the atomic scattering factors.
            By default None, which uses the object's energies.
        **kwargs
            Additional keyword arguments for the `asf_im` or `atomic_scattering` classes.

        Returns
        -------
        asf
            An atomic scattering factor object with the same polynomial coefficients as the piecewise polynomial.

        See Also
        --------
        kkcalc2.models.factors.asf_im : Atomic scattering factor object for the imaginary part.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        # Use kwargs
        from kkcalc2.models.factors import asf_im

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        # Create the asf object
        if target_energies is None:
            return asf_im(
                energies=self.energies,
                factors=self.atomic_scattering_factors,
                **common_kwargs,
            )
        else:
            target_energies = np.asarray(target_energies)
            return asf_im(
                energies=target_energies,
                factors=self.eval_asf(target_energies),
                **common_kwargs,
            )

    # @doc_copy(to_atomic_scattering_factors)
    # def to_asf(self, energies: npt.NDArray | None = None, **kwargs) -> "asf_im":
    #     """
    #     Alias for `to_atomic_scattering_factors`.
    #     """
    #     return self.to_atomic_scattering_factors(energies, **kwargs)

    to_asf = to_atomic_scattering_factors  # Alias for to_atomic_scattering_factors

    def kk_transform(
        self,
        target_energies: npt.NDArray[np.floating | np.integer] | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
    ) -> "asf_re":
        """
        Generate the real part of the atomic scattering factors.

        Uses `kk_algorithms.KK_PP` on the imaginary polynomial coefficients to calculate the real part of
        the atomic scattering factors.
        Can only provide the `stoichiometry` parameter or the `relativistic_correction` parameter, not both.

        Parameters
        ----------
        target_energies : npt.NDArray[np.floating | np.int] | None, optional
            The energies at which to calculate the real part of the atomic scattering factors.
            If None, then the object's energies are used, by default None.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calculate the relativistic correction.
        relativistic_correction : float | None, optional
            The relativistic correction factor to apply to the calculation, by default None.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.

        Returns
        -------
        asp_real
            An `asf_re` object that represents the real part of the atomic scattering factors.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        # Check parameters for/to-define relativistic correction
        energies = self.energies
        coefs = self.coefs
        assert energies is not None and coefs is not None

        if stoichiometry is not None and relativistic_correction is not None:
            raise ValueError(
                "Cannot provide both stoichiometry and relativistic correction."
            )
        elif (
            stoichiometry is None
            and relativistic_correction is None
            and self.stoichiometry is None
        ):
            raise ValueError(
                "Must provide either stoichiometry or relativistic correction."
            )
        # Check argument stoichiometry before using the object's stoichiometry.
        elif stoichiometry is not None:
            relativistic_correction = stoichiometry.relativistic_correction
        elif self.stoichiometry is not None:
            stoichiometry = self.stoichiometry
            relativistic_correction = self.stoichiometry.relativistic_correction
        # At this point, relativistic_correction must be defined.
        assert relativistic_correction is not None

        real_factors = transforms.KK_PP(
            target_energies=(
                target_energies if target_energies is not None else energies
            ),
            energies=energies,
            imag_coefs=coefs,
            relativistic_correction=relativistic_correction,
        )

        # Collate data
        imp_energies = target_energies if target_energies is not None else energies
        imp_real_factors = real_factors

        # Perform accuracy improvement if requested
        if improve_accuracy and max_iter > 0:
            imp_energies, imp_real_factors = transforms.improve_accuracy(
                energies=imp_energies,
                real_asf=imp_real_factors,
                imag_coefs=self.coefs,
                relativistic_correction=relativistic_correction,
                tolerance=tolerance,
                max_iter=max_iter,
            )

        # Import asf_re and create object
        from kkcalc2.models.factors import asf_re

        kwargs = self._properties_dict
        return asf_re(energies=imp_energies, factors=imp_real_factors, **kwargs)

    def calculate_complex_polynomial(
        self,
        target_energies: npt.NDArray | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asp_complex":
        """
        Generate a complex representation of the atomic scattering factors.

        Transforms (`kk_transform`) the imaginary atomic scattering polynomial
        to real factors, and then uses both to form a complex representation.

        Parameters
        ----------
        target_energies : npt.NDArray | None, optional
            The energies at which to calculate the real part of the atomic scattering factors.
            If None, then the object's energies are used, by default None.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calcualte the relativistic correction.
        relativistic_correction : float, optional
            The relativistic correction factor to apply to the calculation, by default False.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.
        **kwargs
            Additional keyword arguments for the `asp_complex` or `atomic_scattering` classes.

        Returns
        -------
        asp_complex
            An atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp_complex : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        kkcalc2.models.polynomials.asp_im.kk_transform : KK transform method.
        """
        from kkcalc2.models.polynomials import asp_complex

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        common_kwargs.update(stoichiometry=stoichiometry)
        # Set default target energies to the object's energies if None
        target_energies = np.asarray(
            target_energies if target_energies is not None else self.energies
        )
        re = self.kk_transform(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        ).to_ASP()
        im = self.extend_energies(re.energies)

        # if the energies don't match, truncate.
        if im.energies.shape != re.energies.shape or not np.all(
            im.energies == re.energies
        ):
            common_e_max = min(im.energies.max(), re.energies.max())
            common_e_min = max(im.energies.min(), re.energies.min())
            im = im.truncate_energies(domain=(common_e_min, common_e_max))
            re = re.truncate_energies(domain=(common_e_min, common_e_max))

        # Create complex object
        return asp_complex(re=re, im=im, **common_kwargs)

    def calculate_complex_factors(
        self,
        target_energies: npt.NDArray | npt.ArrayLike | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT_NO_STOICH],
    ) -> "asf_complex":
        """
        Generate a complex atomic scattering factor object.

        Transforms (`kk_transform`) the imaginary atomic scattering factors to real factors,
        and then uses both to form a complex representation.

        Parameters
        ----------
        target_energies : npt.NDArray | None, optional
            The energies at which to calculate the real part of the atomic scattering factors.
            If None, then the object's energies are used, by default None.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calcualte the relativistic correction.
        relativistic_correction : float, optional
            The relativistic correction factor to apply to the calculation, by default None.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` class.

        Returns
        -------
        asf_complex
            A complex atomic scattering factor object.

        See Also
        --------
        kkcalc2.models.factors.asf_complex : Complex atomic scattering factor object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        kkcalc2.models.polynomials.asp_im.kk_transform : KK transform method.
        """
        from kkcalc2.models.factors import asf_complex

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)  # type: ignore #PROPERTIES_DICT_NO_STOICH is a subset of PROPERTIES_DICT
        # Set default target energies to the object's energies if None
        target_energies = np.asarray(
            target_energies if target_energies is not None else self.energies
        )
        # Calculate the KK transform
        re = self.kk_transform(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        # Evaluate the complex factors at the same energies
        im_asf = self.to_asf(target_energies=re.energies)
        return asf_complex(re=re, im=im_asf, **common_kwargs)

    @overload
    def eval_NEXAFS(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_NEXAFS(
        self, target_energies: float | int
    ) -> float: ...  # numpydoc ignore=GL08

    def eval_NEXAFS(
        self, target_energies: npt.NDArray | npt.ArrayLike | float | int | None = None
    ) -> npt.NDArray | float:
        r"""
        Calculate scattering factors from object polynomial coefficients at desired `energies`.

        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        # Type check
        if target_energies is None:
            # If no energies or coefficients are provided, use the object's values to return the intrinsic ASF values.
            return self.eval_NEXAFS(self.energies)

        if not isinstance(target_energies, (int, float)):
            target_energies = np.asarray(target_energies)
            factors = self.eval_asf_on_coefs(
                target_energies=target_energies,
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )
        else:
            target_energies = np.array([target_energies])
            # Remove the singleton dimension from the output.
            factors = self.eval_asf_on_coefs(
                target_energies=target_energies,
                energies=self.energies,
                coefs=self.coefs,
                orders=self.orders,
            )[0]
        return conversions.ASF_to_NEXAFS(target_energies, factors)

    @overload
    def eval_betas(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_betas(
        self, target_energies: float | int
    ) -> float | complex: ...  # numpydoc ignore=GL08

    def eval_betas(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray | float | complex:
        r"""
        Determine the energy-dependent, imaginary, absorption refractive index component ($\beta$).

        Calculated from the object atomic scattering polynomial coefficients at desired `target_energies`.

        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        return self.eval_refractive(target_energies=target_energies)

    @overload
    def attenuation_length(
        self, energies: npt.NDArray
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def attenuation_length(
        self, energies: float | int
    ) -> float: ...  # numpydoc ignore=GL08

    def attenuation_length(
        self, energies: npt.ArrayLike | npt.NDArray | int | float
    ) -> npt.NDArray | float:
        r"""
        Calculate the attenuation_length for the material at (a) specified energies.

        Length is in angstroms. The attenuation_length is the distance at which light
        is attenuated to an intensity of 1/e.

        .. math::
            I(x) = I_0 e^{-\mu x}
            \mu = 2 \rho \r_e \lambda f_2

        where
        - :math:`f2` is the imaginary atomic scattering factor,
        - :math:`\rho` is the density of the material (g/cm³),
        - :math:`\r_e` is the classical electron radius (2.81794e-15 m),
        - :math:`\lambda` is the wavelength of the radiation (in angstroms),

        Parameters
        ----------
        energies : array_like | npt.NDArray | int | float
            1D array (or singular float) of `M` energies in eV.

        Returns
        -------
        npt.NDArray | float
            The critical angle at energy (or energies) `energies` in radians.
        """
        en = np.asarray(energies)
        if not self.density:
            raise AttributeError(
                "Density information is required to calculate the attenuation length."
            )

        # Calculate the critical angle
        # r_e = sc.value("classical electron radius")  # in meters
        wavelength = sc.h * sc.c / (en * sc.e)  # in meters
        # att_len = self.density * 2 * r_e * self.eval_asf(en) / wavelength
        alpha = 4 * np.pi * self.eval_betas(en) / wavelength  # absorption coefficient
        att_len = 1 / alpha  # penetration depth
        return att_len * 1e10  # Convert m to angstroms.


class asp_re(asp):
    """
    Identical to the `asp` class, but reserved for the real component.

    Enables kk algorithms to convert to real and complex representations of the atomic scattering factors.

    Parameters
    ----------
    energies : npt.ArrayLike
        The energy values of length `N+1` defining the `N` intervals for the polynomial coefficients.
    coefs : npt.ArrayLike
        The polynomial coefficients of shape `(N, M)` for the scattering factors,
        defined on the intervals of `energies` where `M` is the number of coefficients.
    orders : npt.ArrayLike | None, optional
        The polynomial orders for the scattering factors. If None, then kkcalc2 internally
        assumes the polynomial orders are by default [1, 0, -1, -2, -3]. By default None.
        Must have length `M` if provided.
    **kwargs : Unpack[PROPERTIES_DICT], optional
        Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` base class.
    """

    def __init__(
        self,
        energies: npt.ArrayLike,
        coefs: npt.ArrayLike,
        orders: npt.ArrayLike | None = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        super().__init__(energies, coefs, orders, **kwargs)

    @staticmethod
    def from_asp(asp: asp, **kwargs: Unpack[PROPERTIES_DICT]) -> "asp_re":
        """
        Convert an undesignated `asp` object to an `asp_re` object.

        Parameters
        ----------
        asp : asp
            The real part of the atomic scattering factor.
        **kwargs
            Additional keyword arguments for the `asp_re` or `atomic_scattering` classes.

        Returns
        -------
        asp_im
            The imaginary part of the atomic scattering factor.

        See Also
        --------
        kkcalc2.models.polynomials.asp : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        common_kwargs = asp._properties_dict
        common_kwargs.update(kwargs)
        return asp_re(energies=asp.energies, coefs=asp.coefs, **common_kwargs)

    def to_atomic_scattering_factors(
        self,
        target_energies: (
            npt.NDArray[np.floating | np.integer] | npt.ArrayLike | None
        ) = None,
        **kwargs: Unpack[PROPERTIES_DICT],
    ) -> "asf_re":
        """
        Convert the piecewise polynomial representation to an atomic scattering factor object.

        Parameters
        ----------
        target_energies : npt.NDArray[np.floating | np.int_] | np.ArrayLike | None, optional
            Energy values at which to calculate the atomic scattering factors.
            By default None, then the object's energies are used.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` classes.

        Returns
        -------
        asf
            An atomic scattering factor object with the same polynomial coefficients as the piecewise polynomial.

        See Also
        --------
        kkcalc2.models.factors.asf_re : Atomic scattering factor object for the real part.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        from kkcalc2.models.factors import asf_re

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        if target_energies is None:
            return asf_re(
                energies=self.energies,
                factors=self.atomic_scattering_factors,
                **common_kwargs,
            )
        else:
            target_energies = np.asarray(target_energies)
            return asf_re(
                energies=target_energies,
                factors=self.eval_asf(target_energies),
                **common_kwargs,
            )

    # @doc_copy(to_atomic_scattering_factors)
    # def to_asf(self, energies: npt.NDArray | None = None, **kwargs) -> "asf_re":
    #     """
    #     Alias for `to_atomic_scattering_factors`.
    #     """
    #     return self.to_atomic_scattering_factors(energies, **kwargs)

    to_asf = to_atomic_scattering_factors  # Alias for to_atomic_scattering_factors

    def kk_transform_inv(
        self,
        target_energies: npt.ArrayLike | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
    ) -> "asf_im":
        """
        Transform the real polynomial coefficients to imaginary factors.

        Applies the inverse Kramers Kronig transform (`kk_transforms.KK_PP_inv`) to calculate the
        imaginary part of the atomic scattering factors.

        Only one of `stoichiometry` or `relativistic_correction` parameter should be provided, not both.

        Parameters
        ----------
        target_energies : npt.ArrayLike | None, optional
            The energies at which to calculate the imaginary part of the atomic scattering factors.
            If None, then the object's energies are used, by default None.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calcualte the relativistic correction.
        relativistic_correction : float, optional
            The relativistic correction factor to apply to the calculation, by default False.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.

        Returns
        -------
        asf_im
            An `asf_im` object that represents the imaginary part of the atomic scattering factors.
        """
        target_energies = np.asarray(target_energies)

        # Check parameters for/to-define relativistic correction
        if stoichiometry is not None and relativistic_correction is not None:
            raise ValueError(
                "Cannot provide both stoichiometry and relativistic correction."
            )
        elif (
            stoichiometry is None
            and relativistic_correction is None
            and self.stoichiometry is None
        ):
            raise ValueError(
                f"Must provide either stoichiometry or relativistic correction, unless defined on {self}."
            )
        # Check argument stoichiometry before using the object's stoichiometry.
        elif stoichiometry is not None:
            relativistic_correction = stoichiometry.relativistic_correction
        elif self.stoichiometry is not None:
            stoichiometry = self.stoichiometry
            relativistic_correction = self.stoichiometry.relativistic_correction
        else:
            # The relativistic correction is not none.
            assert relativistic_correction is not None

        # Calculate the imaginary part of the atomic scattering factors
        imag_factors = transforms.KK_PP_inv(
            target_energies=(
                target_energies if target_energies is not None else self.energies
            ),
            energies=self.energies,
            real_coefs=self.coefs,
            relativistic_correction=relativistic_correction,
        )

        # Collate "improved" data
        imp_energies = target_energies if target_energies is not None else self.energies
        imp_imag_factors = imag_factors
        imp_real_coefs = self.coefs

        # Perform accuracy improvement if requested
        if improve_accuracy:
            imp_energies, imp_imag_factors = transforms.improve_accuracy_inv(
                energies=imp_energies,
                real_coefs=imp_real_coefs,
                imag_asf=imp_imag_factors,
                relativistic_correction=relativistic_correction,
                tolerance=tolerance,
                max_iter=max_iter,
            )

        # Import asf_im and create object
        from kkcalc2.models.factors import asf_im

        common_kwargs = self._properties_dict
        return asf_im(energies=imp_energies, factors=imp_imag_factors, **common_kwargs)

    def calculate_complex_polynomial(
        self,
        target_energies: npt.NDArray | npt.ArrayLike | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT_NO_STOICH],
    ) -> "asp_complex":
        """
        Calculate a complex polynomial representation of the scattering factors.

        Transforms the real part of the atomic scattering factors into imaginary factors, and then uses both
        to form a complex polynomial representation.

        Uses `kk_transform_inv` default options to perform the transform.

        Parameters
        ----------
        target_energies : npt.ArrayLike | None, optional
            The energies at which to calculate the imaginary part of the atomic scattering factors.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calcualte the relativistic correction.
        relativistic_correction : float, optional
            The relativistic correction factor to apply to the calculation, by default False.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.
        **kwargs
            Additional keyword arguments for the `asp_complex` or `atomic_scattering` classes.

        Returns
        -------
        asp_complex
            An atomic scattering polynomial object.

        See Also
        --------
        kkcalc2.models.polynomials.asp_complex : Atomic scattering polynomial object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        """
        from kkcalc2.models.polynomials import asp_complex

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)  # type: ignore #PROPERTIES_DICT_NO_STOICH is a subset of PROPERTIES_DICT
        common_kwargs.update(stoichiometry=stoichiometry)
        # Set the default target energies to the object's energies if None
        target_energies = np.asarray(
            target_energies if target_energies is not None else self.energies
        )
        # Calculate the imaginary part
        im = self.kk_transform_inv(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        # Extend the energies
        re = self.extend_energies(im.energies)
        return asp_complex(re=re, im=im.to_ASP(), **common_kwargs)

    def calculate_complex_factors(
        self,
        target_energies: npt.ArrayLike | npt.NDArray | None = None,
        improve_accuracy: bool = True,
        stoichiometry: kk_stoichiometry | None = None,
        relativistic_correction: float | None = None,
        tolerance: float = transforms.DEF_TOL,
        max_iter: int = transforms.DEF_ITER,
        **kwargs: Unpack[PROPERTIES_DICT_NO_STOICH],
    ) -> "asf_complex":
        """
        Generate a complex representation by applying the Kramers-Kronig inverse transform.

        Applies the transform to the  real part of the atomic scattering factors,
        to generate imaginary factors.

        Parameters
        ----------
        target_energies : npt.ArrayLike | None, optional
            The energies at which to calculate the imaginary part of the atomic scattering factors.
        improve_accuracy : bool, optional
            If True, then the algorithm will attempt to improve the accuracy of the calculation, by default True.
            This uses the `kk_algorithms.improve_accuracy` algorithm.
        stoichiometry : stoichiometry | None, optional
            The stoichiometry object for the material, by default None
            Used to calcualte the relativistic correction.
        relativistic_correction : float, optional
            The relativistic correction factor to apply to the calculation, by default False.
            Can also be calculated by providing the `stoich` parameter.
        tolerance : float, optional
            Used if `improve_accuracy` is enabled. The tolerance for the accuracy improvement algorithm, by default 1e-2.
        max_iter : int, optional
            Used if `improve_accuracy` is enabled. The maximum number of iterations for the accuracy improvement algorithm, by default 50.
        **kwargs : Unpack[PROPERTIES_DICT_NO_STOICH]
            Additional keyword arguments for `atomic_scattering` classes.

        Returns
        -------
        asf_complex
            A complex atomic scattering factor object.

        See Also
        --------
        kkcalc2.models.factors.asf_complex : Complex atomic scattering factor object.
        kkcalc2.models.common.atomic_scattering : Base class for atomic scattering factors.
        kkcalc2.models.polynomials.asp_re.kk_transform_inv : Inverse KK transform method.
        """
        from kkcalc2.models.factors import asf_complex

        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)  # type: ignore #PROPERTIES_DICT_NO_STOICH is a subset of PROPERTIES_DICT
        common_kwargs.update(stoichiometry=stoichiometry)
        # Set default target energies to the object's energies if None
        target_energies = np.asarray(
            target_energies if target_energies is not None else self.energies
        )
        # Calculate the KK transform
        im = self.kk_transform_inv(
            target_energies=target_energies,
            improve_accuracy=improve_accuracy,
            stoichiometry=stoichiometry,
            relativistic_correction=relativistic_correction,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        # Evaluate the complex factors at the same energies
        return asf_complex(
            re=self.to_atomic_scattering_factors(target_energies=im.energies),
            im=im,
            **common_kwargs,
        )

    @overload
    def critical_angle(
        self, energies: npt.NDArray[np.floating | np.integer]
    ) -> npt.NDArray[np.floating]: ...  # numpydoc ignore=GL08

    @overload
    def critical_angle(
        self, energies: float | int
    ) -> float: ...  # numpydoc ignore=GL08

    def critical_angle(
        self,
        energies: npt.ArrayLike | npt.NDArray[np.floating | np.integer] | int | float,
    ) -> npt.NDArray[np.floating] | float:
        r"""
        Calculate the critical angle for the material at (a) specified energies.

        Angle is in radians. The critical angle is the angle of incidence (from
        the horizon) at which light enters the denser medium at 90 degrees.

        .. math::
            \theta_c = \sqrt{2 \delta}

        Parameters
        ----------
        energies : array_like | npt.NDArray[np.floating | np.int_] | int | float
            1D array (or singular float) of `M` energies in eV.

        Returns
        -------
        npt.NDArray[np.floating] | float
            The critical angle at energy (or energies) `energies` in radians.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        en = np.asarray(energies)
        if self.density is None:
            raise ValueError(
                "Density must be provided via the `asp_re.density` attribute."
            )

        # Calculate the critical angle
        c_angle = np.sqrt(2 * self.eval_refractive(en))
        if c_angle.ndim == 0:
            return c_angle.item()  # If a single value is provided, return as a scalar.
        else:
            return c_angle

    @overload
    def eval_deltas(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_deltas(
        self, target_energies: float | int
    ) -> float | complex: ...  # numpydoc ignore=GL08

    def eval_deltas(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray | float | complex:
        r"""
        Determine the energy-dependent, real, dispersive refractive index component ($\delta$).

        Determined from the object atomic scattering polynomial coefficients at desired `target_energies`.
        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        return self.eval_refractive(target_energies=target_energies)


class asp_complex(asp_abstract, atomic_scattering):
    """
    Container for a combined pair (real & imaginary) of atomic scattering polynomials.

    The default behaviour is to truncate the energies to the common domain, as the
    real and imaginary databases often have different energy domains.

    Parameters
    ----------
    re : asp_re | asp
        The real part of the atomic scattering factor.
    im : asp_im | asp
        The imaginary part of the atomic scattering factor.
    truncate : bool, optional
        If True, then the energies will be truncated to the common domain
        between the real and imaginary parts. Default is True.
    **kwargs : Unpack[PROPERTIES_DICT]
        Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` class.
        Default values are copied from the real part object unless `None` (then the imaginary part object).
        Provided values will override the defaults.
    """

    def __init__(
        self,
        re: asp_re | asp,
        im: asp_im | asp,
        truncate: bool = True,
        **kwargs: Unpack[PROPERTIES_DICT],
    ):  # numpydoc ignore=GL08
        if np.any(re.energies.shape != im.energies.shape) or np.any(
            re.energies != im.energies
        ):
            if not truncate:
                raise ValueError(
                    "Real and imaginary parts must have the same energy domains."
                )
            # While this condition isn't essential, better to have the same energy intervals.
            re_min, re_max, im_min, im_max = (
                re.energies.min(),
                re.energies.max(),
                im.energies.min(),
                im.energies.max(),
            )

            # Truncate to the common interval
            min_energy: float = max(re_min, im_min)
            max_energy: float = min(re_max, im_max)
            all_energies = set(re.energies.tolist() + im.energies.tolist())
            all_energies = [
                en for en in all_energies if en >= min_energy and en <= max_energy
            ]

            # Check if re is a subset of im
            if all(
                [
                    en in im.energies and en >= im_min and en <= im_max
                    for en in re.energies
                ]
            ):
                im = im.extend_energies(re.energies)  # Fill in any additional energies
                im = im.truncate_energies(
                    (min_energy, max_energy)
                )  # Truncate to the common interval
                re = re.extend_energies(im.energies)  # Fill in any additional energies

            # Check if im is a subset of re
            elif all(
                [
                    en in re.energies and en >= re_min and en <= re_max
                    for en in im.energies
                ]
            ):
                re = re.extend_energies(im.energies)  # Fill in any additional energies
                re = re.truncate_energies(
                    (min_energy, max_energy)
                )  # Truncate to the common interval
                im = im.extend_energies(re.energies)  # Fill in any additional energies

            else:  # if they are not subsets, then truncate to the common interval
                re = re.truncate_energies(
                    (min_energy, max_energy)
                )  # Truncate to the common interval
                im = im.truncate_energies(
                    (min_energy, max_energy)
                )  # Truncate to the common interval
                re = re.extend_energies(im.energies)  # Fill in any additional energies
                im = im.extend_energies(re.energies)  # Fill in any additional energies

        if not isinstance(re, asp) or not isinstance(im, asp):
            raise ValueError(f"Real and imaginary parts must be of type {asp}.")

        if (
            re.orders is not None
            and im.orders is not None
            and not np.all(re.orders == im.orders)
        ):
            warnings.warn("Real and imaginary parts have different polynomial orders.")

        # Use the real then imaginary part properties to update None values
        common_kwargs = re._properties_dict

        # Check properties are the same
        for key in im._properties_dict:
            if key not in common_kwargs or common_kwargs[key] is None:
                common_kwargs[key] = im._properties_dict[key]
            elif common_kwargs[key] != im._properties_dict[key]:
                warnings.warn(
                    f"Property {key} is different between real `{re._properties_dict[key]}`"
                    + f" and imaginary parts `{im._properties_dict[key]}` when instantiating `asp_complex`."
                )
            else:
                # Ignore if the properties are the same
                pass

        # Update properties with kwargs
        common_kwargs.update(kwargs)

        # Convert to appropriate instance objects
        if isinstance(re, asp):
            re = asp_re.from_asp(re)
        if isinstance(im, asp):
            im = asp_im.from_asp(im)

        # Store attributes
        self._re: asp_re = re
        self._im: asp_im = im

        # Initialise atomic scattering object
        atomic_scattering.__init__(self, **common_kwargs)

    # Override the attomic scattering properties to propogate to the real & imaginary parts
    @atomic_scattering.density.setter
    def density(self, density: float | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asp_complex, self.__class__).density.fset(self, density)
        # Propogate to components
        self._re.density = density
        self._im.density = density

    @atomic_scattering.number_density.setter
    def number_density(
        self, number_density: float | None
    ) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asp_complex, self.__class__).number_density.fset(self, number_density)
        # Propogate to components
        self._re.number_density = number_density
        self._im.number_density = number_density

    @atomic_scattering.formula_mass.setter
    def formula_mass(self, formula_mass: float | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asp_complex, self.__class__).formula_mass.fset(self, formula_mass)
        # Propogate to components
        self._re.formula_mass = formula_mass
        self._im.formula_mass = formula_mass

    @atomic_scattering.stoichiometry.setter
    def stoichiometry(
        self, stoich: kk_stoichiometry | str | None
    ) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction from atomic_scattering
        super(asp_complex, self.__class__).stoichiometry.fset(self, stoich)
        # Propogate to components
        self._re.stoichiometry = stoich
        self._im.stoichiometry = stoich

    @atomic_scattering.name.setter
    def name(self, name: str | None) -> None:  # numpydoc ignore=GL08
        # Repeat the same instruction
        super(asp_complex, self.__class__).name.fset(self, name)
        # Propogate to components
        self._re.name = name
        self._im.name = name

    @asp_abstract.energies.getter
    def energies(self) -> npt.NDArray:
        """
        The energy intervals for the polynomial coefficients.

        Returns
        -------
        npt.NDArray
            The energy values defining the intervals for the polynomial coefficients.
            Has length `N+1`, where `N` is the number of segments.
        """
        if self._re.energies.shape == self._im.energies.shape and np.all(
            self._re.energies == self._im.energies
        ):
            return self._re.energies
        else:
            return np.union1d(self._re.energies, self._im.energies)

    @property
    def coefs(self) -> npt.NDArray:
        """
        The complex polynomial coefficients for the scattering factors.

        Returns
        -------
        npt.NDArray
            A complex 2D array of shape `(N, M)`, where `N` is the number of segments and `M` is the number of polynomial coefficients.
        """
        return self._re.coefs + 1j * self._im.coefs

    @property
    def orders(self) -> npt.NDArray | None:
        """
        The polynomial orders for the scattering factors, if provided, otherwise `None`.

        Returns
        -------
        npt.NDArray | None
            A 1D array of polynomial orders, with length `M`, where `M` is the number of coefficients.
        """
        re_orders = self._re.orders
        im_orders = self._im.orders
        if re_orders is None and im_orders is None:
            return None
        elif re_orders is None:
            return im_orders
        elif im_orders is None:
            return re_orders
        elif re_orders == im_orders:
            return re_orders
        else:
            raise ValueError(
                "Real and imaginary parts have different polynomial orders, orders cannot be determined."
            )

    @property
    def re(self) -> asp_re:
        """
        The real part object of the atomic scattering polynomial.

        Returns
        -------
        asp_re
            The real part component of the atomic scattering polynomial.
        """
        return self._re

    @property
    def im(self) -> asp_im:
        """
        The imaginary part object of the atomic scattering polynomial.

        Returns
        -------
        asp_im
            The imaginary part component of the atomic scattering polynomial.
        """
        return self._im

    def to_atomic_scattering_factors(
        self, target_energies: npt.ArrayLike | npt.NDArray | None = None, **kwargs
    ) -> "asf_complex":
        """
        Generate an atomic scattering factor object from the piecewise polynomial representation.

        Parameters
        ----------
        target_energies : npt.NDArray | None
            Energy values at which to calculate the atomic scattering factors.
            If None, then the object's energies are used.
        **kwargs
            Additional keyword arguments for the `asf_complex` or `atomic_scattering` classes.

        Returns
        -------
        asf
            An atomic scattering factor object with the same polynomial coefficients as the piecewise polynomial.
        """
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        from kkcalc2.models.factors import asf_complex

        energies = target_energies if target_energies is not None else self.energies
        return asf_complex(
            re=self.re.to_atomic_scattering_factors(energies, **common_kwargs),
            im=self.im.to_atomic_scattering_factors(energies, **common_kwargs),
            **common_kwargs,
        )

    to_asf = to_atomic_scattering_factors  # Alias for to_atomic_scattering_factors

    @overload
    def eval_refractive(  # numpydoc ignore=GL08
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray[np.complex128]: ...

    @overload
    def eval_refractive(
        self,
        target_energies: float | int,  # numpydoc ignore=GL08
    ) -> complex: ...

    def eval_refractive(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray[np.complex128] | complex:
        r"""
        Calculate refractive coefficients from object polynomial coefficients at desired `energies`.

        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.
        .. math::
            = \delta + i\beta

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        if not self.can_calc_refractive:
            raise AttributeError(
                f"{self} cannot calculate delta/beta values; requires density information."
            )

        # Run eval_betas on the real and imaginary parts
        deltas_re = self.re.eval_refractive(target_energies)
        betas_im = self.im.eval_refractive(target_energies)
        return deltas_re + 1j * betas_im

    @overload
    def eval_refractive_index(  # numpydoc ignore=GL08
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray[np.complexfloating]: ...

    @overload
    def eval_refractive_index(  # numpydoc ignore=GL08
        self, target_energies: float | int
    ) -> complex: ...

    def eval_refractive_index(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray[np.complexfloating] | complex:
        r"""
        Calculate the refractive index from the atomic scattering factors at desired `energies`.

        For x-ray energies, the refractive index is calculated as:

        .. math::
            n = 1 - \delta + i\beta

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray[np.complex128] | complex
            The refractive index at energy (or energies) `energies`.
        """
        if not self.can_calc_refractive:
            raise AttributeError(
                f"{self} cannot calculate delta/beta values; requires density information."
            )
        result = self.eval_refractive(target_energies)
        return 1 - result.real + 1j * result.imag

    @overload
    def eval_betas(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_betas(
        self, target_energies: float | int
    ) -> float | complex: ...  # numpydoc ignore=GL08

    def eval_betas(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray | float | complex:
        r"""
        Determine the energy-dependent, imaginary, absorption refractive index component ($\beta$).

        Calculated from the object atomic scattering polynomial coefficients at desired `target_energies`.

        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        if not self.can_calc_refractive:
            raise AttributeError(
                f"{self} cannot calculate delta/beta values; requires density information."
            )
        return self._im.eval_refractive(target_energies=target_energies)

    @overload
    def eval_deltas(
        self, target_energies: npt.NDArray | None
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def eval_deltas(
        self, target_energies: float | int
    ) -> float | complex: ...  # numpydoc ignore=GL08

    def eval_deltas(
        self, target_energies: npt.NDArray | float | int | None = None
    ) -> npt.NDArray | float | complex:
        r"""
        Determine the energy-dependent, real, dispersive refractive index component ($\delta$).

        Determined from the object atomic scattering polynomial coefficients at desired `target_energies`.
        Uses `coefs_to_atomic_scattering_factors` to calculate the ASF values after matching energies to segments.

        Parameters
        ----------
        target_energies : array_like | float, optional
            1D array (or singular float) of `M` energies in eV.
            If None then the energies defined in the object are used.

        Returns
        -------
        npt.NDArray | float
            The magnitude of the atomic scattering factors at energy (or energies) `energies`.
            Dimensions are `M` if `energies` is an array, otherwise a float if `energies` is a float value.
        """
        if not self.can_calc_refractive:
            raise AttributeError(
                f"{self} cannot calculate delta/beta values; requires density information."
            )
        return self._re.eval_deltas(target_energies=target_energies)

    def contrast(self, other: "asp_complex") -> tuple[npt.NDArray, npt.NDArray]:
        r"""
        The energy-dependent contrast magnitude between two complex atomic scattering polynomials.

        The contrast is the difference in imaginary and real parts squared.
        If the objects have different energy domains, only the common domain is considered.

        .. math::
            contrast ~ \Delta(\delta)^2 + \Delta\beta^2

        Parameters
        ----------
        other : asp_complex
            The other atomic scattering polynomial object to compare the contrast.

        Returns
        -------
        energies : np.ndarray
            Array of energy values defined for each contrast value.
        contrast : np.ndarray
            The contrast between two atomic scattering polynomials.
        """
        if self.can_calc_refractive and other.can_calc_refractive:
            energies_self, energies_other = self.energies, other.energies

            if energies_self.shape == energies_other.shape and np.all(
                energies_self == energies_other
            ):
                # All energies are the same, no need to check. Perform a direct calculation.
                energies = energies_self
                betas_self = conversions.ASF_to_refractive(
                    energies=energies,
                    factors=conversions.ASP_to_ASF(energies, self.coefs, self.orders),
                    number_density=self.number_density,
                )
                betas_other = conversions.ASF_to_refractive(
                    energies=energies,
                    factors=conversions.ASP_to_ASF(energies, other.coefs, self.orders),
                    number_density=other.number_density,
                )
                contrast_real = (betas_self.real - betas_other.real) ** 2
                contrast_imag = (betas_self.imag - betas_other.imag) ** 2
                contrast = contrast_real + contrast_imag
                common_energies = energies
            else:
                # Full energies
                energies = np.sort(np.unique(np.r_[energies_self, energies_other]))
                # Get the common lowebound and upperbound
                lower_bound = (
                    energies_other[0]
                    if energies_self[0] < energies_other[0]
                    else energies_self[0]
                )
                upper_bound = (
                    energies_other[-1]
                    if energies_self[-1] > energies_other[-1]
                    else energies_self[-1]
                )
                # Find bounds of the common energies
                lb_idx = np.argmax(energies == lower_bound)
                ub_idx = np.argmax(energies == upper_bound)
                ub_idx = (
                    ub_idx + 1 if ub_idx < len(energies) - 1 else ub_idx
                )  # Include the upper bound

                # Get the common energies
                common_energies = energies[lb_idx:ub_idx]

                # Find the indicies of the starting energies
                self_idx = np.argmax(
                    energies_self > lower_bound
                )  # larger than starting energy
                other_idx = np.argmax(energies_other > lower_bound)
                self_idx = (
                    self_idx - 1 if self_idx > 0 else self_idx
                )  # at or below starting energy
                other_idx = other_idx - 1 if other_idx > 0 else other_idx
                # Iterate over the common energies, and calculate the contrast
                contrast = np.zeros(ub_idx - lb_idx)
                for i in range(ub_idx - lb_idx):
                    # Get values
                    energy = common_energies[i]
                    # Collect the beta values
                    betas_self = conversions.ASF_to_refractive(
                        energies=energy,
                        factors=conversions.ASP_to_ASF(
                            energy, self.coefs[self_idx], orders=self.orders
                        ),
                        number_density=self.number_density,
                        density=self.density,
                        formula_mass=self.formula_mass,
                        stoichiometry=self.stoichiometry,
                    )
                    betas_other = conversions.ASF_to_refractive(
                        energies=energy,
                        factors=conversions.ASP_to_ASF(
                            energy, other.coefs[other_idx], orders=other.orders
                        ),
                        number_density=other.number_density,
                        density=other.density,
                        formula_mass=other.formula_mass,
                        stoichiometry=other.stoichiometry,
                    )
                    # Calculate the contrast
                    contrast_real = (betas_self.real - betas_other.real) ** 2
                    contrast_imag = (betas_self.imag - betas_other.imag) ** 2
                    contrast[i] = contrast_real + contrast_imag

                    # Update the indices, if next energy is reached
                    if energy >= energies_self[self_idx + 1]:
                        self_idx += 1
                    if energy >= energies_other[other_idx + 1]:
                        other_idx += 1
            return common_energies, contrast
        else:
            raise ValueError(
                "Both objects must have beta values to calculate contrast."
            )

    def copy(self, **kwargs: Unpack[PROPERTIES_DICT]) -> Self:
        """
        Generate a copy of the `asp_complex` object.

        Parameters
        ----------
        **kwargs : Unpack[PROPERTIES_DICT]
            Keyword arguments for `asp` and `atomic_scattering` constructors.

        Returns
        -------
        type[asp_complex]
            A new `asp_complex` object with the same polynomial coefficients,
            and properties, but unique memory allocation.
        """
        # Copy the object properties
        common_kwargs = self._properties_dict
        for key in common_kwargs:
            if hasattr(common_kwargs[key], "copy"):
                common_kwargs[key] = common_kwargs[key].copy()
        # Update kwargs
        common_kwargs.update(kwargs)
        return self.__class__(re=self.re.copy(), im=self.im.copy(), **common_kwargs)

    def extend_energies(
        self, energies: npt.NDArray, **kwargs: Unpack[PROPERTIES_DICT]
    ) -> Self:
        """
        Extend the atomic scattering polynomial to include new energy values.

        Parameters
        ----------
        energies : npt.NDArray
            The new energy values to extend the atomic scattering polynomial.
        **kwargs : Unpack[PROPERTIES_DICT]
            Additional keyword arguments for the `atomic_scattering` classes.

        Returns
        -------
        asp_complex
            A new `asp_complex` object with the same polynomial coefficients,
            and properties, but unique memory allocation and additional energy values.
        """
        im_extend = self.im.extend_energies(energies)
        re_extend = self.re.extend_energies(energies)
        common_kwargs = self._properties_dict
        common_kwargs.update(kwargs)
        return self.__class__(re=re_extend, im=im_extend, **common_kwargs)

    @overload
    def critical_angle(
        self, energies: npt.NDArray[np.floating | np.integer]
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def critical_angle(
        self, energies: float | int
    ) -> float: ...  # numpydoc ignore=GL08

    def critical_angle(
        self, energies: npt.NDArray[np.floating | np.integer] | int | float
    ) -> npt.NDArray[np.floating] | float:
        r"""
        Calculate the critical angle for the material at (a) specified energies.

        Angle is in radians. The critical angle is the angle of incidence (from
        the horizon) at which light enters the denser medium at 90 degrees.

        .. math::
            \theta_c = \sqrt{2 \delta}

        Parameters
        ----------
        energies : npt.NDArray[np.floating | np.int_] | int | float
            1D array (or singular float) of `M` energies in eV.

        Returns
        -------
        npt.NDArray | float
            The critical angle at energy (or energies) `energies` in radians.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        return self.re.critical_angle(energies=energies)

    @overload
    def attenuation_length(
        self, energies: npt.NDArray
    ) -> npt.NDArray: ...  # numpydoc ignore=GL08

    @overload
    def attenuation_length(
        self, energies: float | int
    ) -> float: ...  # numpydoc ignore=GL08

    def attenuation_length(
        self, energies: npt.NDArray[np.floating | np.integer] | int | float
    ) -> npt.NDArray[np.floating] | float:
        r"""
        Calculate the attenuation_length for the material at (a) specified energies.

        Length is in angstroms. The attenuation_length is the distance at which light
        is attenuated to an intensity of 1/e.

        .. math::
            I(x) = I_0 e^{-\mu x}
            \mu = 2 \rho \r_e \lambda f_2

        where
        - :math:`f2` is the imaginary atomic scattering factor,
        - :math:`\rho` is the density of the material (g/cm³),
        - :math:`\r_e` is the classical electron radius (2.81794e-15 m),
        - :math:`\lambda` is the wavelength of the radiation (in angstroms),

        Parameters
        ----------
        energies : npt.NDArray[np.floating | np.int_] | int | float
            1D array (or singular float) of `M` energies in eV.

        Returns
        -------
        npt.NDArray | float
            The critical angle at energy (or energies) `energies` in radians.
        """
        # TODO: Use np.integer instead of np.int_ for type hinting, when numpydoc supports it.
        return self.im.attenuation_length(energies=energies)
