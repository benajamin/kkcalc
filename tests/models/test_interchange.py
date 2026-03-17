"""
Model tests for conversions between polynomial and factor representations.

I.e. Between `asf` and `asp` objects.
"""

import numpy as np
import pytest

# import warnings
from kkcalc2.models import (
    asp_db_im,
    asp_db_re,
    # asp_db_im_extended,
    # asp_db_re_extended,
    asp_db_complex,
    # asp_db_complex_extended,
)

from ..test_stoich import basic_stoichs as bs


class TestAsfToAsp:
    """Tests factors to polynomial representation conversion."""


class TestAspToAsf:
    """Tests polynomial representation to factors conversion."""


class TestDbAspToAsf:
    """Tests database polynomial representation to factors conversion."""

    @pytest.mark.parametrize(
        "model",
        [
            asp_db_im,
            asp_db_re,
            asp_db_complex,
        ],
    )
    def test_to_asf(self, model):
        stoich = bs.POLYMER_P3HT
        density = 1.33
        # Create the polynomial object
        poly = model(stoich)
        poly.density = density
        # Convert to asf
        asf = poly.to_asf()
        # Check the asf object
        assert np.isclose(asf.density, density, rtol=1e-5), (
            "The density should match the input density."
        )
        assert asf.stoichiometry == stoich, (
            "The stoichiometry should match the input stoichiometry."
        )
        if isinstance(poly, (asp_db_im, asp_db_complex)):
            assert asf.betas is not None, "The asf object should have betas defined."
        if isinstance(poly, (asp_db_re, asp_db_complex)):
            assert asf.deltas is not None, "The asf object should have deltas defined."
