#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""Unit tests for module `kk.py`."""

from . import kk
import unittest


class KK_Tests(unittest.TestCase):

	"""Unit tests for the `kk` module."""

	def setUp(self):
		pass

	def tearDown(self):
		pass

	def test_calc_relativistic_correction_known_values(self):
		"""Test `calc_relativistic_correction` for known values."""
		known_values = (
			((1, 1, 1, 1, 1), (1, 1, 1, 1, 1), 4.9998565),
			((1, 2, 3, 4, 5), (1, 1, 1, 1, 1), 14.9973659),
			((1, 2, 3, 4, 5), (5, 4, 3, 2, 1), 34.9952628)
			)
		for Z, stoichiometry, expected in known_values:
			result = kk.calc_relativistic_correction(Z, stoichiometry)
			self.assertTrue(abs(expected - result) < 1E-6,
							"Expected %.6f, got %.6f" % (expected, result))

	def test_coeffs_to_ASF_known_values(self):
		"""Test `coeffs_to_ASF` for known values."""
		known_values = (
			(1, (1, 1, 1, 1, 1), 5),
			(2, (1, 1, 1, 1, 1), 3),
			(3, (1, 2, 3, 4, 5), 6),
			(4, (5, 4, 3, 2, 1), 24)
			)
		for E, coeffs, expected in known_values:
			result = kk.coeffs_to_ASF(E, coeffs)
			self.assertTrue(abs(expected - result) < 1E-6,
							"Expected %.6f, got %.6f" % (expected, result))


if __name__ == "__main__":
	unittest.main()
