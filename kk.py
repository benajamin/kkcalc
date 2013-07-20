#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""This module implements the Kramers-Kronig transformation."""

import logging
import math
import numpy


def calc_relativistic_correction(Z, stoichiometry):
	"""Calculate the relativistic correction to the
	Kramers-Kronig transform.

	Parameters:
	-----------
	Z : array of integers
		The list of elements identified by their atomic number
		(not mass).
	stoichiometry : array of integers
		The list of element counts (i.e. the relative proportions
		of the elements).

	Returns
	-------
	This function returns a ``float`` holding the relativistic
	corection to the Kramers-Kronig transform.

	"""
	correction = 0
	for z, s in zip(Z, stoichiometry):
		correction += (z - (z/82.5)**2.37) * s
	return correction


def coeffs_to_ASF(E, coeffs):
	"""Calculate Henke scattering factors from polynomial coefficients.

	Parameters
	----------
	E : float
		Energy in eV
	coeffs : array of floats
		PECS in cm^2/atom

	Returns
	-------
	The function returns the magnitude of the imaginary
	part of the atomic scattering factors at energy `E`.

	"""
	return coeffs[0]*E + coeffs[1] + coeffs[2]/E + coeffs[3]/(E**2) + coeffs[4]/(E**3)



def KK_PP(Energy, imaginary_spectrum, relativistic_correction):
	"""Calculate Kramers-Kronig transform with "Piecewise Polynomial"
	algorithm plus the Biggs and Lighthill extended data.

	Parameters
	----------
	Energy : numpy vector of `float`
		Set of photon energies describing intervals for which each row of `imaginary_spectrum` is valid
	imaginary_spectrum : two-dimensional `numpy.array` of `float`
		The array consists of five columns of polynomial coefficients: A_1, A_0, A_-1, A_-2, A_-3
	relativistic_correction : float
		The relativistic correction to the Kramers-Kronig transform.
		You can calculate the value using the `calc_relativistic_correction` function.

	Returns
	-------
	This function returns the real part of the scattering factors.

	"""
	logger = logging.getLogger(__name__)
	logger.info("Calculate Kramers-Kronig transform using piecewise-polynomial algorithm")
	X1 = Energy[0:-1]
	X2 = Energy[1:]
	E = numpy.tile(Energy, (len(Energy)-1, 1)).T
	Full_coeffs = imaginary_spectrum.T
	Ident = numpy.identity(len(E))  # Use this to annul illegal operations
	temp = (1-(Ident[:, 1:]+Ident[:, 0:-1]))
	Symb_1 = (1-(Ident[:, 1:]+Ident[:, 0:-1]))*((Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)+(Full_coeffs[0, :]*E**2+Full_coeffs[1, :]*E+Full_coeffs[2, :]+Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2-E+Ident[:, 1:])/(X1-E+Ident[:, 0:-1])))-(Full_coeffs[3, :]/E+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))+Full_coeffs[4, :]/E*(X2**-1-X1**-1))
	Symb_2 =                                   (-Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)+(Full_coeffs[0, :]*E**2-Full_coeffs[1, :]*E+Full_coeffs[2, :]-Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2+E)/(X1+E)))                            +(Full_coeffs[3, :]/E-Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))-Full_coeffs[4, :]/E*(X2**-1-X1**-1)
	Symb_B = numpy.sum(Symb_2-Symb_1, axis=1)  # Sum areas for approximate integral
	# Patch singularities
	X1 = E[0:-2, 0]
	XE = E[1:-1, 0]
	X2 = E[2:, 0]
	C1 = Full_coeffs[:, 0:-1]
	C2 = Full_coeffs[:, 1:]
	Symb_singularities = numpy.zeros(len(Energy))
	Symb_singularities[1:-1] = (C2[0, :]*XE**2+C2[1, :]*XE+C2[2, :]+C2[3, :]*XE**-1+C2[4, :]*XE**-2)*numpy.log(numpy.absolute((X2-XE)/(X1-XE)))+(C2[0, :]*XE+C2[1, :])*(X2-XE)+0.5*C2[0, :]*(X2**2-XE**2)-(C2[3, :]*XE**-1+C2[4, :]*XE**-2)*numpy.log(numpy.absolute(X2/XE))+C2[4, :]*XE**-1*(X2**-2-XE**-2)
	Symb_singularities[1:-1] = Symb_singularities[1:-1]+(C1[0, :]*XE+C1[1, :])*(XE-X1)+0.5*C1[0, :]*(XE**2-X1**2)-(C1[3, :]*XE**-1+C1[4, :]*XE**-2)*numpy.log(numpy.absolute(XE/X1))+C1[4, :]*XE**-1*(XE**-2-X1**-2)
	# Finish things off
	KK_Re = (Symb_B-Symb_singularities) / (math.pi*Energy) + relativistic_correction
	logger.debug("Done!")
	return KK_Re

