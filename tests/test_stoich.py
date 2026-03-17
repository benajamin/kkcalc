"""Tests for the `stoich` module."""

import pytest
import periodictable.formulas as ptf
import periodictable as pt
from kkcalc2 import stoichiometry as kk_stoich


class basic_stoichs:
    """Basic functions for creating and testing `stoichiometry` objects."""

    # Define chemical formulae for testing
    POLYMER_PS: str = "C8H8"
    """Polystyrene (PS) chemical formula."""
    POLYMER_P3HT: str = "C10H14S"
    """Poly(3-hexylthiophene-2,5-diyl) (P3HT) chemical formula."""
    POLYMER_PMMA: str = "C5H8O2"
    """Poly(methyl methacrylate) (PMMA) chemical formula."""


class TestStoichiometryInit:
    """Tests for the creation of integer `stoichiometry` objects."""

    @pytest.mark.parametrize(
        "compound,expected",
        [
            (basic_stoichs.POLYMER_PS, [(6, 8), (1, 8)]),
            (basic_stoichs.POLYMER_P3HT, [(6, 10), (1, 14), (16, 1)]),
        ],
    )
    def test_composition(self, compound, expected) -> None:
        """
        Tests the reduction of a str representation.
        """
        result = kk_stoich._parse_chemical_formula(compound)
        result = kk_stoich._consolidate_elements(result)
        assert result == expected

    @pytest.mark.parametrize(
        "input_data",
        [
            # Use some complex polymer such as N2200 (PNDI2OD-T2) for testing
            # Iterable of tuples[int, float]
            [(6, 62), (1, 88), (7, 2), (8, 4), (16, 2)],
            # Periodictable formula contructor
            ptf.formula("C62H88N2O4S2"),
            # Periodictable Formula class
            ptf.Formula(
                (
                    (62, pt.elements[6]),
                    (88, pt.elements[1]),
                    (2, pt.elements[7]),
                    (4, pt.elements[8]),
                    (2, pt.elements[16]),
                )
            ),
            # String
            "C62H88N2O4S2",
            # Another stoichiometry object
            kk_stoich("C62H88N2O4S2"),
            # Dictionary of element symbols to counts
            {"C": 62, "H": 88, "N": 2, "O": 4, "S": 2},
        ],
    )
    def test_format_compilations(self, input_data) -> None:
        """
        Tests the constructor compilations of different types of input.

        Uses N2200 (PNDI2OD-T2) as a test case, which has the chemical formula C62H88N2O4S2.
        """
        stoich_obj = kk_stoich(input_data)
        str_repr = "C62H88N2O4S2"
        str_obj = kk_stoich(str_repr)
        # C62H88N2O4S2
        # Test against the string representation.
        assert stoich_obj == str_obj

    @pytest.mark.parametrize(
        "compound",
        [
            basic_stoichs.POLYMER_PS,
            basic_stoichs.POLYMER_P3HT,
            basic_stoichs.POLYMER_PMMA,
        ],
    )
    def test_stiochiometry(self, compound) -> None:
        """
        Tests the creation of a simple integer `stiochiometry`.
        """
        obj = kk_stoich(compound)
        assert str(obj) == compound


class fractional_stoichs(basic_stoichs):
    """Basic functions for creating and testing fractional `stoichiometry` objects."""

    BLOCK_COPOLYMER_PS_PMMA: str = (
        f"({basic_stoichs.POLYMER_PS})1.0({basic_stoichs.POLYMER_PMMA})1.2"
    )
    """Block copolymer of PS and PMMA 1:1.2 chemical formula."""


class TestStoichiometryFractional:
    """Tests for the creation of fractional `stoichiometry` objects."""

    @pytest.mark.parametrize(
        "frac_str,compound",
        [
            (fractional_stoichs.BLOCK_COPOLYMER_PS_PMMA, "C14H17.6O2.4"),
        ],
    )
    def test_copolymer(self, frac_str, compound) -> None:
        """
        Tests the creation of a fractional `stoichiometry`, where elements get grouped.
        """
        obj = kk_stoich(frac_str)
        assert str(obj) == compound


class TestStoichometryOperations:
    """Tests for the numerical and functional operators of `stoichiometry` objects."""

    @pytest.mark.parametrize(
        "compound",
        [
            fractional_stoichs.POLYMER_PS,
            fractional_stoichs.POLYMER_P3HT,
            fractional_stoichs.POLYMER_PMMA,
            fractional_stoichs.BLOCK_COPOLYMER_PS_PMMA,
        ],
    )
    def test_stiochiometry_numeric(self, compound) -> None:
        """
        Tests the addition, multiplication and division of stoichiometry objects
        """
        # Create the integer components
        polymer = kk_stoich(compound)

        # Test multiplication operations
        assert 0.5 * polymer == kk_stoich(f"({polymer})0.5")

        # Test division
        assert polymer / 2 == kk_stoich(f"({polymer})0.5")
        composition = polymer.composition
        for i in range(len(composition)):
            composition[i] = (composition[i][0], composition[i][1] / 2)
        assert polymer / 2 == kk_stoich(composition)

        # Test floor division
        composition = polymer.composition
        for i in range(len(composition)):
            composition[i] = (composition[i][0], composition[i][1] // 2)
        assert (polymer // 2) == kk_stoich(composition)

    def test_block_numeric(self) -> None:
        """
        Tests the addition of copolymers to form a block
        """
        # Create the integer components
        ps = kk_stoich(fractional_stoichs.POLYMER_PS)
        pmma = kk_stoich(fractional_stoichs.POLYMER_PMMA)
        # Create the block
        block = kk_stoich(fractional_stoichs.BLOCK_COPOLYMER_PS_PMMA)
        # Test numeric operations
        assert (1.0 * ps + 1.2 * pmma) == block

    @pytest.mark.parametrize(
        "compound",
        [
            fractional_stoichs.POLYMER_PS,
            fractional_stoichs.POLYMER_P3HT,
            fractional_stoichs.POLYMER_PMMA,
            fractional_stoichs.BLOCK_COPOLYMER_PS_PMMA,
        ],
    )
    def test_stoichiometry_equivalence(self, compound) -> None:
        """
        Tests the equivalence of various forms of `stiochiometry` objects.
        """
        #
        polymer = kk_stoich(compound)

        # Test equivalence operations
        assert polymer == compound
        assert compound == polymer  # reverse
        assert polymer == polymer.copy()
        assert polymer == 1.0 * polymer
        # Check quantity of elements affects equivalence
        assert polymer != kk_stoich(
            [(elem[0], elem[1] * 1.1) for elem in polymer.composition]
        )

    @pytest.mark.parametrize(
        "compound,length",
        [
            (fractional_stoichs.POLYMER_PS, 16),
            (fractional_stoichs.POLYMER_P3HT, 25),
            (fractional_stoichs.POLYMER_PMMA, 15),
            (fractional_stoichs.BLOCK_COPOLYMER_PS_PMMA, 33),
        ],
    )
    def test_stoichiometry_len(self, compound, length) -> None:
        """
        Tests the atomic length of a `stoichiometry` object.
        """
        # Create the component
        polymer = kk_stoich(compound)
        # Test length
        assert len(polymer) == length
