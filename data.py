#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""This module implements data formats and conversions."""

import logging
import json, numpy


# Constants
CLASSICAL_ELECTRON_RADIUS = 2.81794029957951365441605230194258e-15  # meters
PLANCKS_CONSTANT = 4.1356673310e-15  # eV*seconds
SPEED_OF_LIGHT = 2.99792458e8  # meters per second
AVOGADRO_CONSTANT = 6.02214129e23  # no unit


def convert_Beta_to_ASF(raw_data, density):
	"""Convert Beta data (index of refraction) to atomic scattering
	factors (ASF).

	The Beta value is the imaginary part of the index of refraction.
	This represents the absorption.

	Parameters
	----------
	raw_data : two-dimensional `numpy.array` of `float`
		The array consists of two columns: Energy and magnitude.

	Returns
	-------
	The function returns a `numpy.array` of atomic scattering factors.
	They are made up of the energy and the magnitude of the imaginary
	part of the atomic scattering factors.

	"""
	raw_Im = raw_data[:, :2].copy()
	raw_Im[:, 1] = density*AVOGADRO_CONSTANT*2*numpy.pi*raw_data[:, 0]**2*raw_data[:, 1]/(self.MolecularMass*CLASSICAL_ELECTRON_RADIUS*(PLANCKS_CONSTANT*SPEED_OF_LIGHT)**2)
	return raw_Im


def convert_BL_to_ASF(E, coeffs, Atomic_mass):
	"""Convert Biggs and Lighthill (BL) to atomic scattering factors (ASF).

	Biggs and Lighthill offers photoelectric cross-section (PECS) with the
	sum of AnE^-n for n=1-4  {E in keV and PECS in cm^2/g}.

	Henke scattering factors related by
	f2 = PECS*E/(2*r0*h*c)  {E in eV and PECS in cm^2/atom}.

	Parameters
	----------
	E : float
		Energy in eV
	coeffs : array of float
		Polynomial coefficients describing the function of the
		scattering factors
	Atomic_mass : float
		Atomic mass of the element (in units = g/mol)

	Returns
	-------
	The function returns the magnitude of the imaginary
	part of the atomic scattering factors at energy `E`.

	"""
	return (coeffs[0] + coeffs[1]/(E*0.001) + coeffs[2]/((E*0.001)**2) + coeffs[3]/((E*0.001)**3))*Atomic_mass/(2*AVOGADRO_CONSTANT*CLASSICAL_ELECTRON_RADIUS*PLANCKS_CONSTANT*SPEED_OF_LIGHT)*0.1


def convert_NEXAFS_to_ASF(raw_data):
	"""Convert NEXAFS photoabsorption data to atomic scattering factors (ASF).

	Parameters
	----------
	raw_data : two-dimensional `numpy.array` of `float`
		The array consists of two columns: Energy and magnitude.

	Returns
	-------
	The function returns a `numpy.array` of atomic scattering factors.
	They are made up of the energy and the magnitude of the imaginary
	part of the atomic scattering factors.

	"""
	raw_Im = raw_data[:, :2].copy()
	raw_Im[:, 1] = raw_data[:, 0]*raw_data[:, 1]/(2*CLASSICAL_ELECTRON_RADIUS*PLANCKS_CONSTANT*SPEED_OF_LIGHT)
	return raw_Im

def load_Element_Database():
	with open('ASF.json','r') as f:
		Element_Database = json.load(f)
	for Z in range(1,93):
		Element_Database[str(Z)]['E'] = numpy.array(Element_Database[str(Z)]['E'])
		Element_Database[str(Z)]['Im'] = numpy.array(Element_Database[str(Z)]['Im'])
		Element_Database[str(Z)]['Re'] = numpy.array(Element_Database[str(Z)]['Re'])
	return Element_Database
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	