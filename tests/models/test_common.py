"""
Model tests for polynomial and factor representations.
"""

import pytest
import warnings
from kkcalc2 import models, stoichiometry as kk_stoich

from ..test_stoich import fractional_stoichs as fs


class TestCommon:
    @pytest.mark.parametrize(
        "kwargs, msgs",
        [
            (dict(), []),
            (dict(name="Sample"), []),
            (dict(name="Sample", number_density=1.1), []),
            (dict(name="Sample", number_density=1.1, stoichiometry="CH"), []),
            (
                dict(
                    name="Sample",
                    number_density=1.1,
                    stoichiometry="CH",
                    density=1.8,
                ),
                ["Competing information"],
            ),
        ],
    )
    def test_instantiation_atomic_scattering(
        self, kwargs: dict, msgs: list[str]
    ) -> None:
        """Tests the creation of an `atomic_scattering` object, with expected errors."""
        # Create the object
        with warnings.catch_warnings(record=True) as w:
            _ = models.atomic_scattering(**kwargs)

        # Check each msg is included by at least one warning
        for msg in msgs:
            assert len(w) > 0
            assert any(msg in str(warn.message) for warn in w)

    # Define some mass values for testing
    MASS_VALUES = {
        "1": 1.00784,  # Hydrogen
        "2": 4.0026,  # Helium
        "3": 6.94,  # Lithium
        "4": 9.0122,  # Beryllium
        "5": 10.81,  # Boron
        "6": 12.011,  # Carbon
        "7": 14.007,  # Nitrogen
        "8": 15.999,  # Oxygen
        "16": 32.06,  # Sulfur
        "17": 35.45,  # Chlorine
    }

    @pytest.mark.parametrize("composition", [fs.POLYMER_PS])
    def test_propogation(self, composition: kk_stoich.COMPOSITION_TYPING) -> None:
        """Tests the propogation of the `atomic_scattering` object."""
        # Create the object
        atomic_scattering = models.atomic_scattering(
            name="Sample",
            number_density=1.1,
            stoichiometry=composition,
        )
        # Parse the stoichiometry:
        stoich = kk_stoich(composition)
        # Get the propogated formula mass
        fm = atomic_scattering.formula_mass
        # Calculate the fomula mass from known values
        fm_calc = sum(
            [
                self.MASS_VALUES[str(atom)] * counts
                for atom, counts in stoich.composition
            ]
        )

        assert f"{fm:0.2f}" == f"{fm_calc:0.2f}", (
            f"Formula mass {fm} != {fm_calc} at 2 decimal places."
        )

        old_density = atomic_scattering.density
        # Modify the formula mass
        atomic_scattering.formula_mass = fm * 2
        # Remove the stoichiometry
        atomic_scattering.stoichiometry = None
        # Check the density has updated
        assert old_density * 2 == atomic_scattering.density, (
            "Density did not update correctly after formula mass change."
        )


class TestPolynomial:
    """Tests the functionality for the polynomial models."""
