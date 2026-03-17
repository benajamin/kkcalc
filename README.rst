=====================================
 KKCalc2
=====================================

``kkcalc2`` is a comprehensive toolkit for calculating Kramers Kronig transforms of X-ray absorption/dispersion data, and is built to the feature-rich standards of `xraysoftmat <https://github.com/xraysoftmat>`_.

|tool-semver| |tool-black| |tool-ruff| |tool-numpydoc|

|PyPI Version2| |readthedocs| |Coveralls| |Pre-commit|

|PyTest| |Linting| |Documentation|

.. |PyPI Version| image:: https://img.shields.io/pypi/v/kkcalc?label=KKCalc&logo=pypi
   :target: https://pypi.org/project/kkcalc/
   :alt: pypi
.. |PyPI Version2| image:: https://img.shields.io/pypi/v/kkcalc2?label=KKCalc2&logo=pypi
    :target: https://pypi.org/project/kkcalc2/
    :alt: pypi
.. |PyTest| image:: https://github.com/xraysoftmat/kkcalc/actions/workflows/tests.yml/badge.svg
    :alt: PyTest
    :target: https://github.com/xraysoftmat/kkcalc/actions/workflows/test.yml
.. |Linting| image:: https://github.com/xraysoftmat/kkcalc/actions/workflows/linting.yml/badge.svg
    :alt: Linting
    :target: https://github.com/xraysoftmat/kkcalc/actions/workflows/linting.yml
.. |Documentation| image:: https://github.com/xraysoftmat/kkcalc/actions/workflows/docs.yml/badge.svg
    :alt: Documentation
    :target: https://github.com/xraysoftmat/kkcalc/actions/workflows/docs.yml
.. |Coveralls| image:: https://img.shields.io/coverallsCoverage/github/xraysoftmat/kkcalc?branch=v2&label=Coveralls
    :alt: Coverage Status
    :target: https://coveralls.io/github/xraysoftmat/kkcalc?branch=v2
.. |Pre-commit| image:: https://results.pre-commit.ci/badge/github/xraysoftmat/kkcalc/v2.svg
    :alt: pre-commit.ci status
    :target: https://results.pre-commit.ci/latest/github/xraysoftmat/KKCalc/v2
.. |readthedocs| image:: https://img.shields.io/readthedocs/kkcalc?version=latest&style=flat&label=ReadtheDocs
    :alt: Documentation
    :target: https://kkcalc.readthedocs.io/

.. |tool-semver| image:: https://img.shields.io/badge/versioning-Python%20SemVer-blue.svg
    :alt: Python SemVer
    :target: https://python-semantic-release.readthedocs.io/en/stable/
.. |tool-black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :alt: Code style: black
    :target: https://github.com/psf/black
.. |tool-ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
    :alt: Ruff
    :target: https://github.com/astral-sh/ruff
.. |tool-numpydoc| image:: https://img.shields.io/badge/doc_style-numpydoc-blue.svg
    :alt: Code doc: numpydoc
    :target: https://github.com/numpy/numpydoc

Introduction
############

``kkcalc2`` is an open-source `python` package to calculate the Kramers-Kronig (inverse) transform of X-ray absorption (dispersion) data:

.. note: mathjax github requires double backslashes to properly compute.

$$f_2(E) = \\frac{2}{\\pi} P \\int_{0}^{\\infty}\\frac{x f_1(x)}{x^2 - E^2} dx + \\mathcal{Z}^\\star$$

where $f_1$ and $f_2$ are the real and imaginary parts of the complex index of refraction, respectively, $\\mathcal{Z}^\\star$ is the relativistic correction, and $P$ denotes the Cauchy principal value at ($x=E$).

``kkcalc2`` uses a polynomial representation algorithm developed by Watts [1]_.

This package provides an object oriented API, to evaluate optical constants (index of refraction, absorption and dispersion, etc.), extend measurement spectra with databases, or can be accessed through a PyQT6 GUI interface. Documentation can be found at  `readthedocs <https://kkcalc.rtfd.org/>`_, and releases (including documentation and executable builds) can be found at `github <https://github.com/xraysoftmat/kkcalc>`_.

References
==========

.. [1] Benjamin Watts, "Calculation of the Kramers-Kronig transform of X-ray spectra by a piecewise Laurent polynomial method", *Opt. Express* **22**, (2014) 23628-23639. `DOI:10.1364/OE.22.023628 <https://doi.org/10.1364/OE.22.023628>`_

.. We use the optical constants databases from Henke et al. [2]_, Biggs and Lighthill [3]_.

..
    .. [2] B.L. Henke, E.M. Gullikson, and J.C. Davis, "X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92", *Atomic Data and Nuclear Data Tables* **54** (2) (1993) 181-342 `DOI:10.1006/adnd.1993.1013 <https://doi.org/10.1006/adnd.1993.1013>`_.

    .. [3] F. Biggs, and R. Lighthill, "Analytical approximations for X-ray cross-sections III", *Sandia Report* SAND87-0070 UC-34 (1988). `DOI:10.2172/7124946 <https://doi.org/10.2172/7124946>`_
