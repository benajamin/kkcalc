#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module implements the Kramers-Kronig transformation.
"""

import logging
import math
import numpy
import scipy.fftpack


# Constants
CLASSICAL_ELECTRON_RADIUS = 2.81794029957951365441605230194258e-15  # meters
PLANCKS_CONSTANT = 4.1356673310e-15  # eV*seconds
SPEED_OF_LIGHT = 2.99792458e8  # meters per second
AVOGADRO_CONSTANT = 6.02214129e23  # no unit


def calc_relativistic_correction(Z, stoichiometry):
	"""Calculate the relativistic correction to the
	Kramers-Kronig transform.

	Parameters:
	-----------
	Z : array_like
		TODO: explain parameter
	stoichiometry : array_like
		TODO: explain parameter

	Returns
	-------
	This function returns a ``float`` holding the relativistic
	corection to the Kramers-Kronig transform.

	"""
	correction = 0
	for i in xrange(len(Z)):
		correction += (Z[i] - (Z[i]/82.5)**2.37) * stoichiometry[i]
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
	Henke scattering factor.

	"""
	return coeffs[0]*E + coeffs[1] + coeffs[2]/E + coeffs[3]/(E**2) + coeffs[4]/(E**3)


def KK_FFT(merged_Im, relativistic_correction):
	"""Calculate Kramers-Kronig transform with FFT algorithm.

	Parameters
	----------
	merged_Im : two dimensional array_like
		TODO: explain with more detail
	relativistic_correction : float
		The relativistic correction to the Kramers-Kronig transform.
		You can calculate the value using the
		`calc_relativistic_correction` function.

	Returns
	-------
	TODO

	"""
	logger = logging.getLogger(__name__)
	logger.info("Calculate Kramers-Kronig transform (FFT)")
	E_step = merged_Im[1:-1, 0] - merged_Im[0:-2, 0]
	FFT_step = min(E_step)/2
	even_E = numpy.arange(0, merged_Im[-1, 0], FFT_step)
	N_E = len(even_E)
	N_pow2 = 2**math.ceil(math.log(N_E, 2))
	logger.debug("Use step size of %d eV (%d + %d points)" %
				 (FFT_step, N_E, N_pow2-N_E))
	logger.debug("Interpolate (Part 1/4)")
	even_Im = numpy.interp(even_E, merged_Im[:, 0], merged_Im[:, 1],
						   left=0, right=0)
	del even_E  # maximise available memory, we can recalculate even_E when needed.
	logger.debug("First FFT (Part 2/4)")
	temp_wave = 2 / math.pi * scipy.fftpack.fft(even_Im, n=2*N_pow2).imag
	del even_Im
	temp_wave[N_pow2:-1] = -temp_wave[N_pow2:-1]
	logger.debug("Second FFT (Part 3/4)")
	temp_wave = scipy.fftpack.ifft(temp_wave)
	temp_wave = temp_wave[0:N_E]
	temp_wave = math.pi*temp_wave.real + relativistic_correction
	logger.debug("Reinterpolate (Part 4/4)")
	even_E = numpy.arange(0, merged_Im[-1, 0], FFT_step)
	KK_Re = numpy.interp(merged_Im[:, 0], even_E, temp_wave)
	logger.debug("Done!")
	return KK_Re


def KK_PP(merged_Im, relativistic_correction):
	"""Calculate Kramers-Kronig transform with
	"Piecewise Polynomial" algorithm.

	Parameters
	----------
	merged_Im : two dimensional array_like
		TODO: explain with more detail
	relativistic_correction : float
		The relativistic correction to the Kramers-Kronig transform.
		You can calculate the value using the
		`calc_relativistic_correction` function.

	Returns
	-------
	TODO

	"""
	logger = logging.getLogger(__name__)
	logger.info("Calculate Kramers-Kronig transform (PP)")
	len_E = len(merged_Im[:, 0])
	X1 = merged_Im[0:-1, 0]
	X2 = merged_Im[1:, 0]
	Y1 = merged_Im[0:-1, 1]
	Y2 = merged_Im[1:, 1]
	M = (Y2-Y1)/(X2-X1)
	E = numpy.tile(merged_Im[:, 0], (len_E-1, 1)).T
	# Find areas between data points, assuming linear interpolation of Im data
	Symb_1 = Y1*(X2-X1)+M*(0.5*(X2**2-X1**2)-(X1-E)*(X2-X1))+E*(Y1-M*(X1-E))*numpy.log(numpy.absolute((X2-E)/(X1-E)))
	Symb_2 = Y1*(X2-X1)+M*(0.5*(X2**2-X1**2)-(X1+E)*(X2-X1))-E*(Y1-M*(X1+E))*numpy.log(numpy.absolute((X2+E)/(X1+E)))
	Symb_1[~numpy.isfinite(Symb_1)] = 0  # Ignore singularities for now
	Symb_A = numpy.sum(Symb_2-Symb_1, axis=1)  # Sum areas for approximate integral
	del X1, X2, Y1, Y2, M, E, Symb_1, Symb_2
	# Patch singularities by integrating across two intervals at once, avoiding evaluation at the singularity.
	# Note we will not calculate this at the end-points and assume it is zero (due to symmetry)
	X1 = merged_Im[0:-2, 0]
	XE = merged_Im[1:-1, 0]
	X2 = merged_Im[2:, 0]
	Y1 = merged_Im[0:-2, 1]
	YE = merged_Im[1:-1, 1]
	Y2 = merged_Im[2:, 1]
	M1 = (YE-Y1)/(XE-X1)
	M2 = (Y2-YE)/(X2-XE)
	Symb_singularities = numpy.zeros(len_E)
	Symb_singularities[1:-1] = 0.5*M2*(X2**2-XE**2)-0.5*M1*(X1**2-XE**2)+YE*((X2-X1)+XE*numpy.log(numpy.absolute((X2-XE)/(X1-XE))))
	# Finish things off
	KK_Re = (Symb_A - Symb_singularities) / (math.pi * merged_Im[:, 0]) + relativistic_correction
	logger.debug("Done!")
	return KK_Re


def KK_PP_BL(merged_Im, relativistic_correction, BL_coefficients, BL_range):
	"""Calculate Kramers-Kronig transform with "Piecewise Polynomial"
	algorithm plus the Biggs and Lighthill extended data.

	Parameters
	----------
	merged_Im : two dimensional array_like
		TODO: explain with more detail
	relativistic_correction : float
		The relativistic correction to the Kramers-Kronig transform.
		You can calculate the value using the
		`calc_relativistic_correction` function.
	BL_coefficients : ?
		TODO: explain
	BL_range : ?
		TODO: explain

	Returns
	-------
	TODO

	"""
	logger = logging.getLogger(__name__)
	logger.info("Calculate Kramers-Kronig transform (PP) plus BL data")
	len_E = len(merged_Im[:, 0])
	M = (merged_Im[1:, 1] - merged_Im[0:-1, 1]) / (merged_Im[1:, 0] - merged_Im[0:-1, 0])
	B = merged_Im[0:-1, 1] - M*merged_Im[0:-1, 0]
	E = merged_Im[:, 0]
	Full_coeffs = numpy.zeros((len_E-1, 5))
	Full_coeffs[:, 0] = M
	Full_coeffs[:, 1] = B
	# B&L extension
	C = BL_coefficients
	X = BL_range
	E = E[0:-1]
	for i in range(len(X)-1):
		Y1 = coeffs_to_ASF(X[i]-0.2, Full_coeffs[-1, :])
		Y2 = coeffs_to_ASF(X[i]+0.2, C[i, :])
		M = (Y2-Y1)/0.4
		B = Y1-M*(X[i]-0.2)
		E = numpy.append(E, [X[i]-0.2, X[i]+0.2])
		Full_coeffs = numpy.append(Full_coeffs, [[M, B, 0, 0, 0]], axis=0)
		Full_coeffs = numpy.append(Full_coeffs, [C[i, :]], axis=0)
	E = numpy.append(E, X[-1])
	X1 = E[0:-1]
	X2 = E[1:]
	E = numpy.tile(E, (len(E)-1, 1)).T
	Full_coeffs = Full_coeffs.T
	Ident = numpy.identity(len(E))  # Use this to annul illegal operations
	temp = (1-(Ident[:, 1:]+Ident[:, 0:-1]))
	Symb_1 = (1-(Ident[:, 1:]+Ident[:, 0:-1]))*((Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)+(Full_coeffs[0, :]*E**2+Full_coeffs[1, :]*E+Full_coeffs[2, :]+Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2-E+Ident[:, 1:])/(X1-E+Ident[:, 0:-1])))-(Full_coeffs[3, :]/E+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))+Full_coeffs[4, :]/E*(X2**-1-X1**-1))
	Symb_2 =                                   (-Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)+(Full_coeffs[0, :]*E**2-Full_coeffs[1, :]*E+Full_coeffs[2, :]-Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2+E)/(X1+E)))                          +(Full_coeffs[3, :]/E-Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))-Full_coeffs[4, :]/E*(X2**-1-X1**-1)
	Symb_B = numpy.sum(Symb_2-Symb_1, axis=1)  # Sum areas for approximate integral
	# Patch singularities
	X1 = E[0:-2, 0]
	XE = E[1:-1, 0]
	X2 = E[2:, 0]
	C1 = Full_coeffs[:, 0:-1]
	C2 = Full_coeffs[:, 1:]
	Symb_singularities = numpy.zeros(len(E))
	Symb_singularities[1:-1] = (C2[0, :]*XE**2+C2[1, :]*XE+C2[2, :]+C2[3, :]*XE**-1+C2[4, :]*XE**-2)*numpy.log(numpy.absolute((X2-XE)/(X1-XE)))+(C2[0, :]*XE+C2[1, :])*(X2-XE)+0.5*C2[0, :]*(X2**2-XE**2)-(C2[3, :]*XE**-1+C2[4, :]*XE**-2)*numpy.log(numpy.absolute(X2/XE))+C2[4, :]*XE**-1*(X2**-2-XE**-2)
	Symb_singularities[1:-1] = Symb_singularities[1:-1]+(C1[0, :]*XE+C1[1, :])*(XE-X1)+0.5*C1[0, :]*(XE**2-X1**2)-(C1[3, :]*XE**-1+C1[4, :]*XE**-2)*numpy.log(numpy.absolute(XE/X1))+C1[4, :]*XE**-1*(XE**-2-X1**-2)
	# Finish things off
	cut = 2 * (len(BL_range) - 1)  # remove calculated values at energies higher than 30 keV
	KK_Re = (Symb_B[:-cut]-Symb_singularities[:-cut]) / (math.pi*E[:-cut, 0]) + relativistic_correction
	logger.debug("Done!")
	return KK_Re
