"""
Module to define classes and properties relating to the chemical composition of materials.

Includes
"""

# Stdlib
from __future__ import (
    annotations,
)  # Required to allow union of string and class in type hints
import re
from typing import Self, TYPE_CHECKING, TypeAlias, Iterable, Unpack

# External
has_periodictable: bool
"""Flag to indicate if the periodictable module is available."""
try:
    import periodictable as pt
    from periodictable.formulas import Formula as Formula
    from periodictable.core import Element as Element

    has_periodictable = True
except ImportError:
    has_periodictable = False

# Internal
from kkcalc2.asf_database import ASF_DATABASE, ASFElement  # noqa: E402

if TYPE_CHECKING:
    # Do not compile at runtime due to circular import.
    from kkcalc2.models.db_models import asp_db_im, asp_db_re, asp_db_complex
    from periodictable.formulas import Formula
    from kkcalc2.models.common import (
        PROPERTIES_DICT,
    )

# Generate a list of atomic elements. Should already be sorted from the periodictable module.
ELEMENTS: list[tuple[str, int, float]]
"""A list of tuples containing the atomic symbol, atomic number and atomic mass of each element."""
if has_periodictable:
    # Also contains N=0, i.e. neutral, as the first element. So ELEMENTS[1] = H.
    ELEMENTS = (
        # Neutron first element, not included in the __iter__ method currently... changed in 2024...
        [(pt.elements[0].symbol, pt.elements[0].number, pt.elements[0].mass)]
        + [(element.symbol, element.number, element.mass) for element in pt.elements]
    )
    assert ELEMENTS[1][0] == "H", f"Second element should be H, was {ELEMENTS[1][0]}"
else:
    # Use the asf database
    db = ASF_DATABASE
    atomic_nums = sorted(db.keys())
    """Atomic Numbers"""
    atomic_syms = []
    """Atomic Symbols"""
    atomic_masses = []
    """Atomic Masses"""
    for Z in atomic_nums:
        a: ASFElement = db[Z]
        atomic_syms.append(a["symbol"])
        atomic_masses.append(a["mass"])

    ELEMENTS = [
        (
            "n",
            0,
            1.008,
        ),  # Neutron first element, so H is ELEMENTS[1], consistent with periodictable.
        *zip(atomic_syms, atomic_nums, atomic_masses),
    ]
    assert ELEMENTS[1][0] == "H", f"Second element should be H, was {ELEMENTS[1][0]}"


def relativistic_correction_eq(composition: list[tuple[int, float]]) -> float:
    r"""
    Calculate a relativistic correction to the Kramers-Kronig transform owing to the elemental composition.

    Automatically calculable for a `stoichiometry` using `stoichiometry.relativistic_correction`. Each element contributes
    (z - (z/82.5)**2.37) * n to the correction, where z is the atomic number and n is the relative stoichiometry.

    .. math::
        \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i

    Parameters
    ----------
    composition : list[tuple[int, float]]
        A list of tuples, where each tuple contains the atomic number and the counts of an element.
        Counts may be fractional.

    Returns
    -------
    float
        The relativistic corection to the Kramers-Kronig transform.
    """
    return sum([(z - (z / 82.5) ** 2.37) * n for z, n in composition])


CompositionAlias: TypeAlias = (
    Iterable[tuple[int, float]]
    | "Formula"
    | str
    | "stoichiometry"
    | dict[str, int | float]
)


class stoichiometry:
    """
    Defines the stoichiometry of a chemical compound.

    Internally uses a periodictable.formulas.Formula object,
    or a list of tuples to represent the composition of a compound.

    Parameters
    ----------
    composition : list[tuple[int, float]] | Formula | str | stoichiometry | dict[str, float]
        The stoichiometry of the compound, i.e. the elemental composition.
        Can be a list of tuples, a Formula object, a string,
        a dictionary of element symbols to counts, or another stoichiometry object.

        Examples:
        - [(6, 9), (1, 12), (8, 6), (16, 2)] for C9H12O6S2
        - "C9H12O6S2" for C9H12O6S2
        - pt.formula("C9H12O6S2") for C9H12O6S2
        - "(A)1.2(B)0.8" for a combined composition.

    Attributes
    ----------
    TYPING : TypeAlias
        Type alias for the composition parameter and property.

    See Also
    --------
    stoichiometry.from_chemical_formula : Convert a chemical formula string to a stoichiometry object.
    periodictable.formulas.Formula : Formula object from the periodictable package.
    """

    COMPOSITION_TYPING = CompositionAlias
    """Type alias for the composition parameter."""

    def __init__(self, composition: CompositionAlias) -> None:  # numpydoc ignore=GL08
        self._composition: list[tuple[int, float]] | "Formula"
        """A list of tuples, where each tuple contains the atomic number and the counts of an element."""
        if isinstance(composition, type(self)):
            # Copy the formula / list.
            c = composition._composition
            if isinstance(c, Formula):
                self._composition = Formula(c)
            else:
                self._composition = c.copy()
        elif has_periodictable and isinstance(composition, Formula):
            self._composition = composition
        elif isinstance(composition, str):
            # Convert string to composition.
            c = stoichiometry._parse_chemical_formula(composition)
            c = stoichiometry._consolidate_elements(c)
            self._composition = c
        elif isinstance(composition, dict):
            # Convert dict to composition.
            c = []
            for elem, n in composition.items():
                if isinstance(elem, str):
                    # Convert element symbol to atomic number.
                    elem = stoichiometry._element_to_atomic_number(elem)
                if elem < 1 or elem > 92:
                    raise ValueError("Atomic number out of range.")
                if n < 0:
                    raise ValueError("Negative stoichiometry.")
                c.append((elem, n))
            c = stoichiometry._consolidate_elements(c)
            self._composition = c
        else:
            # Test if iterable
            iterable = False
            try:
                iter(composition)
                iterable = True
            except TypeError:
                raise ValueError(f"Invalid stoichiometry {composition}.")
            if iterable:
                # Check validity of composition, collect duplicate elements
                # Also operates for dicts, as they have the .items() method.
                final_comp: list[tuple[int, int | float]] = []
                for elem, n in composition:
                    if elem < 1 or elem > 92:
                        raise ValueError("Atomic number out of range.")
                    if n < 0:
                        raise ValueError("Negative stoichiometry.")
                    # Check if element is already accounted for
                    exists: bool = False
                    for i, (elem2, n2) in enumerate(final_comp):
                        if elem == elem2:
                            final_comp[i] = (elem, n + n2)
                            exists = True
                            continue
                    if not exists:
                        final_comp.append((elem, n))
                self._composition = final_comp

    def __eq__(self, other: Self | str) -> bool:
        """
        Compare the stoichiometry of two compounds.

        Comparison is made by calling the `composition` property, rather than the `_composition` attribute.

        Parameters
        ----------
        other : stoichiometry | str
            The stoichiometry to compare with the current stoichiometry.
            Can also be a string representation of a stoichiometry.

        Returns
        -------
        bool
            True if the stoichiometry of the compounds are equal, False otherwise.
        """
        if isinstance(other, self.__class__):
            return self.composition == other.composition
        elif isinstance(other, str):
            # Try to convert the string to a stoichiometry object.
            try:
                stoich = stoichiometry(other)
                return self.composition == stoich.composition
            except ValueError:
                # Try to convert
                return str(self) == other
        return False

    def __req__(self, other: Self | str) -> bool:
        """
        Compare the stoichiometry of two compounds.

        Comparison is made by calling the `composition` property, rather than the `_composition` attribute.

        Parameters
        ----------
        other : stoichiometry | str
            The stoichiometry to compare with the current stoichiometry.
            Can also be a string representation of a stoichiometry.

        Returns
        -------
        bool
            True if the stoichiometry of the compounds are equal, False otherwise.
        """
        return self.__eq__(other)

    def __str__(self) -> str:
        """
        Generate a string representation of the stoichiometry.

        Returns
        -------
        str
            A string representation of the stoichiometry.

        Examples
        --------
        >>> str(stoichiometry("C9H12") + stoichiometry("O6S2"))
        'C9H12O6S2'
        """
        return "".join(
            [
                # Get the element symbol
                ELEMENTS[element[0]][0]
                # Get the number of atoms
                + (
                    # Float
                    str(
                        element[1]
                    )  # Could possibly use fractions.Fraction here to display 1/3 etc.
                    # Check if the number is an integer and != 1.
                    if (element[1] * 10) % 10 != 0
                    # Integer
                    else (str(int(element[1])) if element[1] != 1 else "")
                )
                for element in self._composition
            ]
        )

    def __repr__(self) -> str:
        return f"stoichiometry({self._composition})"

    @property
    def elements(self) -> list[str]:
        """
        Return a string list of the elements present in the stoichiometry.

        Returns
        -------
        list[str]
            A list of the elemental symbols present in the stoichiometry.
        """
        return [ELEMENTS[element][0] for element, _ in self.composition]

    def __len__(self) -> float | int:
        """
        Return the summed number of each element in the stoichiometry.

        Returns
        -------
        float | int
            The summed number of each element in the stoichiometry. Can be fractional.
        """
        return sum([int(count) for _, count in self.composition])

    def __add__(self, other: Self | str) -> Self:
        """
        Combine two stoichiometries, or a stoichiometry and a string.

        Parameters
        ----------
        other : stoichiometry | str
            The stoichiometry to combine with the current stoichiometry.
            Can also be a string representation of a stoichiometry.

        Returns
        -------
        stoichiometry
            A new stoichiometry object with the composition of the two combined.
        """
        if isinstance(other, str):
            return self.__class__(self.composition + stoichiometry(other).composition)
        return self.__class__(self.composition + other.composition)

    def __radd__(self, other: Self | str) -> Self:
        """
        Combine two stoichiometries, or a stoichiometry and a string.

        Parameters
        ----------
        other : stoichiometry | str
            The stoichiometry to combine with the current stoichiometry.
            Can also be a string representation of a stoichiometry.

        Returns
        -------
        stoichiometry
            A new stoichiometry object with the composition of the two combined.
        """
        return self.__add__(other)

    def __iadd__(self, other: Self | str) -> None:
        """
        Add another stoichiometry to the current stoichiometry.

        Parameters
        ----------
        other : stoichiometry
            The stoichiometry to add to the current stoichiometry.

        Returns
        -------
        stoichiometry
            The current stoichiometry object with the composition of the two combined.
        """
        initial_comp = self.composition
        if isinstance(other, str):
            other_comp = stoichiometry(other).composition
        # Check validity of composition, collect duplicate elements
        final_comp: list[tuple[int, int | float]] = []
        for elem, n in initial_comp + other_comp:
            # Check if element is already accounted for
            exists: bool = False
            for i, (elem2, n2) in enumerate(final_comp):
                if elem == elem2:
                    final_comp[i] = (elem, n + n2)
                    exists = True
                    continue
            if not exists:
                final_comp.append((elem, n))
        self._composition = final_comp

    def __mul__(self, other: float) -> Self:
        """
        Multiply the stoichiometry by a scalar.

        Parameters
        ----------
        other : float
            The scalar to multiply the stoichiometry by.

        Returns
        -------
        Self
            A new stoichiometry object with the composition multiplied by the scalar.
        """
        return self.__class__(
            [(elem, count * other) for elem, count in self.composition]
        )

    def __rmul__(self, other: float) -> Self:
        """
        Reflection operator for multiplication.

        Parameters
        ----------
        other : float
            The scalar to multiply the stoichiometry by.

        Returns
        -------
        Self
            A new stoichiometry object with the composition multiplied by the scalar.
        """
        return self.__mul__(other)

    def __truediv__(self, other: float) -> Self:
        """
        Calculate the true division of a stoichiometry object by a scalar.

        Parameters
        ----------
        other : float
            The scalar to divide the stoichiometry by.

        Returns
        -------
        stoichiometry
            A new stoichiometry object with the composition divided by the scalar.

        Examples
        --------
        >>> stoichiometry("C9H12O6S2") / 2
        C4.5H6O3S1
        """
        return self.__mul__(1 / other)

    def __floordiv__(self, other: float) -> Self:
        """
        Calculate the floor division of a stoichiometry object by a scalar.

        Parameters
        ----------
        other : float
            The scalar by which to divide the stoichiometry.

        Returns
        -------
        stoichiometry
            A new stoichiometry object with the composition divided by the scalar.

        Examples
        --------
        >>> stoichiometry("C9H12O6S2") // 2
        C4H6O3S1
        """
        return self.__class__(
            [(elem, int(count // other)) for elem, count in self.composition]
        )

    def copy(self) -> Self:
        """
        Create a copy of stoichiometry.

        Returns
        -------
        stoichiometry
            A copy of the stoichiometry object, with a unique composition reference.
        """
        return self.__class__(self.composition.copy())

    @property
    def composition(self) -> list[tuple[int, float]]:
        """
        The stoichiometry of the compound, i.e. the elemental composition.

        Returns
        -------
        list[tuple[int, float]]
            A list of tuples, where each tuple contains the atomic number and the counts of an element.
            Counts may be fractional.
        """
        if has_periodictable and isinstance(self._composition, Formula):
            c = []
            element: Element
            count: float
            for element, count in self._composition.atoms.items():
                c.append((element.number, count))
            return c
            # return [(element.number, count) for element, count in self._composition.atoms.items()]
        elif isinstance(self._composition, list):
            return self._composition.copy()
        else:
            raise ValueError("Composition is not a valid type.")

    @property
    def relativistic_correction(self) -> float:
        r"""
        Calculate the relativistic correction to the Kramers-Kronig transform owing to the elemental composition.

        Uses `stoich.relativistic_correction_eq`. Each element contributes (z - (z/82.5)**2.37) * n to the correction,
        where z is the atomic number and n is the relative stoichiometry.

        .. math::
            \mathcal{Z}^\star = \sum_i (Z_i - (Z_i/82.5)^{2.37}) \cdot n_i

        Returns
        -------
        float
            The relativistic corection to the Kramers-Kronig transform.
        """
        return relativistic_correction_eq(composition=self.composition)

    @property
    def formula_mass(self) -> float:
        """
        The sum of atomic masses.

        Returns
        -------
        float
            The sum of atomic masses for the given stoichiometry.
        """
        if has_periodictable:
            return sum(
                [
                    number * pt.elements[element].mass
                    for element, number in self.composition
                ]
            )
        else:
            return sum(
                [number * ELEMENTS[element][2] for element, number in self.composition]
            )

    def atomic_scattering_polynomial_im(
        self, **kwargs: Unpack["PROPERTIES_DICT"]
    ) -> "asp_db_im":
        """
        Generate a piecewise polynomial of the imaginary atomic scattering factors for the given stoichiometry.

        Uses the energy-dependent atomic scattering factor data from the Henke, Briggs and Lighthill database.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` such as:
            - `number_density` : float
            - `density` : float
            - `stoich` : stoichiometry
            - `formula_mass` : float
            - `name` : str

        Returns
        -------
        asp_db_im
            An object representing the piecewise polynomial calculated from the summation of scattering factor data.
        """
        from kkcalc2.models.db_models import asp_db_im

        return asp_db_im(self, **kwargs)

    asp_im = (
        atomic_scattering_polynomial_im  # Alias for atomic_scattering_polynomial_im
    )

    def atomic_scattering_polynomial_re(
        self, **kwargs: Unpack["PROPERTIES_DICT"]
    ) -> "asp_db_re":
        """
        Generate a piecewise polynomial of the real atomic scattering factors for the given stoichiometry.

        Uses the energy-dependent atomic scattering factor data from the Henke, Briggs and Lighthill database.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` such as:
            - `number_density` : float
            - `density` : float
            - `stoich` : stoichiometry
            - `formula_mass` : float
            - `name` : str

        Returns
        -------
        asp_db_re
            An object representing the dispersive piecewise polynomial calculated from the summation of scattering factor data.
        """
        from kkcalc2.models.db_models import asp_db_re

        return asp_db_re(self, **kwargs)

    asp_re = (
        atomic_scattering_polynomial_re  # Alias for atomic_scattering_polynomial_re
    )

    def asp_complex(self, **kwargs: Unpack["PROPERTIES_DICT"]) -> "asp_db_complex":
        """
        Generate a piecewise polynomial of the complex atomic scattering factors for the given stoichiometry.

        Uses the energy-dependent atomic scattering factor data from the Henke, Briggs and Lighthill database.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments for the `kkcalc2.models.common.atomic_scattering` such as:
            - `number_density` : float
            - `density` : float
            - `stoich` : stoichiometry
            - `formula_mass` : float
            - `name` : str

        Returns
        -------
        asp_db_complex
            An object representing the complex piecewise polynomial calculated from the summation of scattering factor data.
        """
        from kkcalc2.models.db_models import asp_db_complex

        return asp_db_complex(self, **kwargs)

    atomic_scattering_polynomial_complex = asp_complex

    @staticmethod
    def _consolidate_elements(
        composition: list[tuple[int, float]],
    ) -> list[tuple[int, float]]:
        """
        Consolidate a list of elements and quantities into a unique list of elements and quantities.

        Parameters
        ----------
        composition : list[tuple[int, float]]
            A list of tuples, where each tuple contains the atomic number and the counts of an element.

        Returns
        -------
        list[tuple[int, float]]
            A list of tuples, where each tuple contains the atomic number and the counts of an element.
        """
        # Setup a dictionary to store the composition
        consolidated = {}
        for element, count in composition:
            if element in consolidated:
                consolidated[element] += count
            else:
                consolidated[element] = count
        return [(element, count) for element, count in consolidated.items()]

    @staticmethod
    def _parse_chemical_formula(
        formula: str, recursion: bool = True
    ) -> list[tuple[int, float]]:
        """
        Convert a chemical compound string into a list of elements and quantities.

        Parameters
        ----------
        formula : str
            A string consisting of element symbols, numbers and parentheses.
        recursion : bool, optional
            Flag to enable recursion in the parsing of the formula string, by default True.

        Returns
        -------
        list[tuple[int, int]]
            A list of tuples, where each tuple contains the atomic number and the counts of an element.
        """
        # Setup a list to store the composition
        composition = []
        ## Regex explaination:
        # ?P<groupname> is a named group to capture.
        # Here we 1st capture either an element symbol or a parenthesized group.
        # Then we capture a number (if present) and the remainder of the formula.
        # <Paren> or <Remainder> groups are then also processed by a recursive call.
        # +? is a non-greedy match, to capture the smallest possible group.
        search = re.compile(
            "".join(
                [
                    r"((?P<Element>[A-Z][a-z]?)|\((?P<Paren>.+?)\))",
                    r"(?P<Number>\d*(\.\d+)?)(?P<Remainder>.*)",
                ]
            )
        )
        # Perform the search on the formula
        m = re.search(search, formula)
        if m is None:
            raise ValueError(f"No formula match: {formula}")
        # Process the search.
        if len(m.group("Number")) != 0:
            Number = float(m.group("Number"))
        else:
            Number = 1.0
        if m.group("Element") is not None:
            Z = stoichiometry._element_to_atomic_number(m.group("Element"))
            if Z != 0:
                composition.append((Z, Number))
        elif len(m.group("Paren")) > 0:
            composition += [
                (x[0], x[1] * Number)
                for x in stoichiometry._parse_chemical_formula(
                    m.group("Paren"), recursion=recursion
                )
            ]
        if len(m.group("Remainder")) != 0:
            composition += stoichiometry._parse_chemical_formula(
                m.group("Remainder"), recursion=recursion
            )
        return composition

    @staticmethod
    def from_chemical_formula(
        formula: str, recursion: bool = True, use_peroidictable: bool = True
    ) -> "stoichiometry":
        """
        Parse a chemical formula string to obtain a stoichiometry.

        Parameters
        ----------
        formula : str
            A string consisting of element symbols, numbers and parentheses.
        recursion : bool, optional
            Whether to use recursion to parse the formula, by default True.
        use_peroidictable : bool, optional
            Whether to use the periodictable module to parse the formula, by default True.

        Returns
        -------
        stoichiometry
            A stoichiometry object representing the composition of the formula.
        """
        if use_peroidictable:
            return stoichiometry(pt.formula(formula))
        else:
            # Parse the formula string
            composition = stoichiometry._parse_chemical_formula(
                formula=formula, recursion=recursion
            )
            # Consolidate the elements
            composition = stoichiometry._consolidate_elements(composition)
            # Create the stoichiometry object
            return stoichiometry(composition)

    @staticmethod
    def _element_to_atomic_number(SymbolString: str) -> int:
        """
        Replace list of elemental symbols with the corresponding atomic numbers.

        Parameters
        ----------
        SymbolString : str
            An elemental symbol (i.e. "H", "C", "O", etc.).

        Returns
        -------
        int
            The function returns an integer atomic number corresponding to the input symbol.
            Zero is returned when the string is not recognised.
        """
        for i in range(len(ELEMENTS)):
            if ELEMENTS[i][0] == SymbolString:
                return ELEMENTS[i][1]
        raise ValueError(f"`{SymbolString}` is not a known element!")

    @staticmethod
    def _atomic_number_to_element(Z: int) -> str:
        """
        Replace list of atomic numbers with the corresponding elemental symbols.

        Parameters
        ----------
        Z : int
            Integer representing an atomic number.

        Returns
        -------
        str
            The function returns a string elemental symbol corresponding to the input atomic number.
        """
        # Z'th list item should match the element.
        if ELEMENTS[Z][1] == Z:
            return ELEMENTS[Z][0]
        # If not, search for the element index. This should not be necessary.
        for i in range(len(ELEMENTS)):
            if ELEMENTS[i][1] == Z:
                return ELEMENTS[i][0]
        raise ValueError(f"Element #{Z} is not a known atomic number to kkcalc!")


if __name__ == "__main__":
    # Test the stoichiometry class
    P3MEET1 = "C9H12O6S2"  # C9H11O3S
    P3MEET2 = "(C9H12O6S2)0.1(C9H11O3S)0.9"
    if has_periodictable:
        P3MEET3 = pt.formula("C9H12O6S2")

    compounds = [P3MEET1, P3MEET2] + ([P3MEET3] if has_periodictable else [])
    data_titles = [
        "Stoichiometry",
        "Composition",
        "Relativistic Correction",
        "Formula Mass",
    ]
    data = []

    for c, compound in enumerate(compounds):
        stoich = stoichiometry(compound)
        comp = stoich.composition
        for i, (atom, count) in enumerate(comp):
            if type(count) is float and int(count) != count:
                comp[i] = (atom, float(f"{count:.2f}"))  # Round to 3 decimal places
        data.append(
            [compound, comp, stoich.relativistic_correction, stoich.formula_mass]
        )
        if c == 1:
            print(
                f"Testing bracketed formula: {stoich.composition[0]} == {9 * 0.1 + 9 * 0.9}? {stoich.composition[0][1] == 9 * 0.1 + 9 * 0.9}"
            )

    import pandas as pd

    df = pd.DataFrame(data, columns=data_titles)
    print(df)
