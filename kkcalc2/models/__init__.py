"""
This module is the main entry point for the models of the Kramers-Kronig Calculator.

It consists of the following submodules:
        * common: Contains the common material models used by the other models
        * factors: Contains the models for the atomic scattering factors (`asf`)
        * polynomials: Contains the models for polynomial respresentations of the `asf`.
        * conversions: A library of functions for the conversions between datatypes.
        * db_models: Contains the models for accessing and using the atomic scattering
                     factor database.

Additionally, the database models are imported from the `kkcalc2.asf_database.db_models` module,
for model completeness.
"""

# Import the conversions module
from kkcalc2 import conversions
from kkcalc2 import transforms as transforms

# Import the common base models
from kkcalc2.models.common import (
    atomic_scattering_abstract,
    atomic_scattering,
    PROPERTIES_DICT,
    PROPERTIES_DICT_NO_STOICH,
)

# Import the usage models
import kkcalc2.models.polynomials as polynomials
import kkcalc2.models.factors as factors
import kkcalc2.models.common as common
import kkcalc2.models.db_models as db_models
from kkcalc2.models.polynomials import asp, asp_im, asp_re, asp_complex, asp_abstract
from kkcalc2.models.factors import (
    asf_abstract,
    asf,
    asf_im,
    asf_re,
    asf_complex,
    KK_Datatype,
)

# Import the database models
from kkcalc2.models.db_models import (
    asp_db_abstract,
    asp_db_im,
    asp_db_re,
    asp_db_complex,
    asp_db_extended,
    asp_db_im_extended,
    asp_db_re_extended,
    asp_db_complex_extended,
)

__all__ = [
    # Conversions module
    "conversions",
    # Kramers-Kronig Transforms module
    "transforms",
    # Class modules
    "polynomials",
    "factors",
    "common",
    "db_models",
    # Common models and types
    "atomic_scattering_abstract",
    "atomic_scattering",
    "PROPERTIES_DICT",
    "PROPERTIES_DICT_NO_STOICH",
    # Atomic Scattering Factor models and types
    "asf_abstract",
    "asf",
    "asf_im",
    "asf_re",
    "asf_complex",
    "KK_Datatype",
    # Atomic Scattering Polynomial models and types
    "asp_abstract",
    "asp",
    "asp_im",
    "asp_re",
    "asp_complex",
    # Database models
    "asp_db_abstract",
    "asp_db_im",
    "asp_db_re",
    "asp_db_complex",
    "asp_db_extended",
    "asp_db_im_extended",
    "asp_db_re_extended",
    "asp_db_complex_extended",
]
