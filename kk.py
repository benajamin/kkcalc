#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""This module implements the Kramers-Kronig transformation."""

import logging, sys
logger = logging.getLogger(__name__)
if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG)
	logging.StreamHandler(stream=sys.stdout)

import math
import numpy
import os
import data


def calc_relativistic_correction(stoichiometry):
	"""Calculate the relativistic correction to the Kramers-Kronig transform.

	Parameters:
	-----------
	stoichiometry : array of integer/float pairs
		Each pair in the list consists of an atomic number and the relative proportion of that element.

	Returns
	-------
	This function returns a ``float`` holding the relativistic
	corection to the Kramers-Kronig transform.

	"""
	correction = 0
	for z, n in stoichiometry:
		correction += (z - (z/82.5)**2.37) * n
	return correction


def KK_PP(Eval_Energy, Energy, imaginary_spectrum, relativistic_correction):
	"""Calculate Kramers-Kronig transform with "Piecewise Polynomial"
	algorithm plus the Biggs and Lighthill extended data.

	Parameters
	----------
	Eval_Energy : numpy vector of `float`
		Set of photon energies describing points at which to evaluate the real spectrum
	Energy : numpy vector of `float`
		Set of photon energies describing intervals for which each row of `imaginary_spectrum` is valid
	imaginary_spectrum : two-dimensional `numpy.array` of `float`
		The array consists of five columns of polynomial coefficients: A_1, A_0, A_-1, A_-2, A_-3
	relativistic_correction : float
		The relativistic correction to the Kramers-Kronig transform.
		You can calculate the value using the `calc_relativistic_correction` function.

	Returns
	-------
	This function returns the real part of the scattering factors evaluated at photon energies specified by Eval_Energy.

	"""
	logger = logging.getLogger(__name__)
	logger.info("Calculate Kramers-Kronig transform using piecewise-polynomial algorithm")
	X1 = Energy[0:-1]
	X2 = Energy[1:]
	E = numpy.tile(Eval_Energy, (len(Energy)-1, 1)).T
	Full_coeffs = imaginary_spectrum.T
	Symb_1 = (( Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)-(Full_coeffs[3, :]/E+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))+Full_coeffs[4, :]/E*(X2**-1-X1**-1))
	Symb_2 = ((-Full_coeffs[0, :]*E+Full_coeffs[1, :])*(X2-X1)+0.5*Full_coeffs[0, :]*(X2**2-X1**2)+(Full_coeffs[3, :]/E-Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute(X2/X1))-Full_coeffs[4, :]/E*(X2**-1-X1**-1))+(Full_coeffs[0, :]*E**2-Full_coeffs[1, :]*E+Full_coeffs[2, :]-Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2+E)/(X1+E)))
	Symb_3 = (1-1*((X2==E)|(X1==E)))*(Full_coeffs[0, :]*E**2+Full_coeffs[1, :]*E+Full_coeffs[2, :]+Full_coeffs[3, :]*E**-1+Full_coeffs[4, :]*E**-2)*numpy.log(numpy.absolute((X2-E+1*(X2==E))/(X1-E+1*(X1==E))))
	Symb_B = numpy.sum(Symb_2 - Symb_1 - Symb_3, axis=1)  # Sum areas for approximate integral
	# Patch singularities
	hits = Energy[1:-1]==E[:,0:-1]
	E_hits = numpy.insert(numpy.any(hits, axis=0),[0,0],False)
	Eval_hits = numpy.any(hits, axis=1)
	X1 = Energy[E_hits[2:]]
	XE = Energy[E_hits[1:]]
	X2 = Energy[E_hits]
	C1 = Full_coeffs[:, E_hits[2:]]
	C2 = Full_coeffs[:, E_hits[1:]]
	Symb_singularities = numpy.zeros(len(Eval_Energy))
	Symb_singularities[Eval_hits] = (C2[0, :]*XE**2+C2[1, :]*XE+C2[2, :]+C2[3, :]*XE**-1+C2[4, :]*XE**-2)*numpy.log(numpy.absolute((X2-XE)/(X1-XE)))
	# Finish things off
	KK_Re = (Symb_B-Symb_singularities) / (math.pi*Eval_Energy) + relativistic_correction
	logger.debug("Done!")
	return KK_Re

def kk_calculate_real(NearEdgeDataFile, ChemicalFormula, load_options=None, input_data_type=None, merge_points=None, add_background=False, fix_distortions=False):
	"""Do all data loading and processing and then calculate the kramers-Kronig transform.
	Parameters
	----------
	NearEdgeDataFile : string
		Path to file containg near-edge data
	ChemicalFormula : string
		A standard chemical formula string consisting of element symbols, numbers and parentheses.
	merge_points : list or tuple pair of `float` values, or None
		The photon energy values (low, high) at which the near-edge and scattering factor data values
		are set equal so as to ensure continuity of the merged data set.

	Returns
	-------
	This function returns a numpy array with columns consisting of the photon energy, the real and the imaginary parts of the scattering factors.
	"""
	Stoichiometry = data.ParseChemicalFormula(ChemicalFormula)
	Relativistic_Correction = calc_relativistic_correction(Stoichiometry)
	ASF_E, ASF_Data = data.calculate_asf(Stoichiometry)
	NearEdge_Data = data.convert_data(data.load_data(NearEdgeDataFile, load_options),FromType=input_data_type,ToType='asf')
	Full_E, Imaginary_Spectrum = data.merge_spectra(NearEdge_Data, ASF_E, ASF_Data, merge_points=merge_points, add_background=add_background, fix_distortions=fix_distortions)
	Real_Spectrum = KK_PP(Full_E, Full_E, Imaginary_Spectrum, Relativistic_Correction)
	
	Imaginary_Spectrum_Values = data.coeffs_to_ASF(Full_E, numpy.vstack((Imaginary_Spectrum,Imaginary_Spectrum[-1])))
	return numpy.vstack((Full_E,Real_Spectrum,Imaginary_Spectrum_Values)).T, Imaginary_Spectrum

if __name__ == '__main__':
	#use argparse here to get command line arguments
	#process arguments and pass to a pythonic function
	
	#I will abuse this section of code for initial testing
	#Output, Im = kk_calculate_real('data/Xy_norm_bgsub.txt', 'C10SH14', input_data_type='NEXAFS')
	Output, Im = kk_calculate_real('data/LaAlO3_Exp.csv', 'LaAlO3', input_data_type='NEXAFS', fix_distortions=True)
	#Output, Im = kk_calculate_real('data/As.xmu.csv', 'GaAs', input_data_type='NEXAFS', fix_distortions=True)
	
	Stoichiometry = data.ParseChemicalFormula('LaAlO3')
	#Stoichiometry = data.ParseChemicalFormula('GaAs')
	Relativistic_Correction = calc_relativistic_correction(Stoichiometry)
	ASF_E, ASF_Data = data.calculate_asf(Stoichiometry)
	ASF_Data3 = data.coeffs_to_linear(ASF_E, ASF_Data, 0.1)
	ASF_Data2 = data.coeffs_to_ASF(ASF_E, numpy.vstack((ASF_Data,ASF_Data[-1])))
	
	#Test_E = (Output[1:,0]+Output[0:-1,0])*0.5
	#Test_E = numpy.linspace(41257.87,41259.87,num=21)
	#Real_Spectrum2 = KK_PP(Test_E, Output[:,0], Im, Relativistic_Correction)
	
	
	import matplotlib
	matplotlib.use('WXAgg')
	import pylab
	
	pylab.figure()
	pylab.plot(Output[:,0],Output[:,1],'xg-',Output[:,0],Output[:,2],'xb-')
	pylab.plot(ASF_E,ASF_Data2,'+r')
	#pylab.plot(ASF_E,ASF_Data22,'xr')
	pylab.plot(ASF_Data3[0],ASF_Data3[1],'r-')
	#pylab.plot(Test_E,Real_Spectrum2,'*y')
	pylab.xscale('log')
	pylab.show()
	
	
	