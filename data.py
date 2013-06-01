#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module implements data formats and conversions."""

import logging


# Constants
CLASSICAL_ELECTRON_RADIUS = 2.81794029957951365441605230194258e-15  # meters
PLANCKS_CONSTANT = 4.1356673310e-15  # eV*seconds
SPEED_OF_LIGHT = 2.99792458e8  # meters per second
AVOGADRO_CONSTANT = 6.02214129e23  # no unit


def convert_Beta_to_ASF(raw_data, density):
	"""Convert Beta data to atomic scattering factors (ASF).

	Parameters
	----------
	raw_data : two-dimensional array of `float`
		TODO: explain

	Returns
	-------
	TODO: explain

	"""
	raw_Im = raw_data[:, :2].copy()
	raw_Im[:, 1] = density*AVOGADRO_CONSTANT*2*numpy.pi*raw_data[:, 0]**2*raw_data[:, 1]/(self.MolecularMass*CLASSICAL_ELECTRON_RADIUS*(PLANCKS_CONSTANT*SPEED_OF_LIGHT)**2)
	return raw_Im


def convert_BL_to_ASF(E, coeffs, Atomic_mass):
	"""Convert Biggs and Lighthill (BL) to atomic scattering factors (ASF).

	Biggs and Lighthill offers photoelectric cross-section with the
	sum of AnE^-n for n=1-4  {E in keV and PECS in cm^2/g}.

	Henke scattering factors related by
	f2 = PECS*E/(2*r0*h*c)  {E in eV and PECS in cm^2/atom}.

	Parameters
	----------
	E : float
		Energy in eV
	coeffs : array of float
		TODO: explain
	Atomic_mass : float
		Atomic mass of the element (in units = g/mol)

	Returns
	-------
	TODO: explain

	"""
	return (coeffs[0] + coeffs[1]/(E*0.001) + coeffs[2]/((E*0.001)**2) + coeffs[3]/((E*0.001)**3))*Atomic_mass/(2*AVOGADRO_CONSTANT*CLASSICAL_ELECTRON_RADIUS*PLANCKS_CONSTANT*SPEED_OF_LIGHT)*0.1


def convert_NEXAFS_to_ASF(raw_data):
	"""Convert NEXAFS photoabsorption data to atomic scattering factors (ASF).

	Parameters
	----------
	raw_data : two-dimensional array of `float`
		TODO: explain

	Returns
	-------
	TODO: explain

	"""
	raw_Im = raw_data[:, :2].copy()
	raw_Im[:, 1] = raw_data[:, 0]*raw_data[:, 1]/(2*CLASSICAL_ELECTRON_RADIUS*PLANCKS_CONSTANT*SPEED_OF_LIGHT)
	return raw_Im
