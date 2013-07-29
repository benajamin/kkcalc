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
import re


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
	"""Loads atomic scattering factor database from a json file.
	
	The database has been previously created by PackRawData.py
	
	Parameters
	----------
	None
	
	Returns
	-------
	The function returns a dictionary of elements, each consisting of a dictionary
	of data types

	"""
	with open('ASF.json','r') as f:
		Element_Database = json.load(f)
	for Z in range(1,93):
		Element_Database[str(Z)]['E'] = numpy.array(Element_Database[str(Z)]['E'])
		Element_Database[str(Z)]['Im'] = numpy.array(Element_Database[str(Z)]['Im'])
		Element_Database[str(Z)]['Re'] = numpy.array(Element_Database[str(Z)]['Re'])
	return Element_Database


def ConvertElementSymbol(SymbolString):
	"""Replace list of elemental symbols with the corresponding atomic numbers.

	Parameters
	----------
	SymbolString : String representing an elemental symbol

	Returns
	-------
	The function returns an integer atomic number corresponding to the input symbol.
	Zero is returned when the string is not recognised.
	"""
	Elements = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U']
	try:
		AtomicNumber = Elements.index(SymbolString)+1
	except ValueError:
		print "\""+SymbolString+"\" is not a known element!"
		AtomicNumber = 0
	return AtomicNumber

def ParseChemicalFormula(Formula):
	"""Parse a chemical formula string to obtain a stoichiometry.

	Parameters
	----------
	Formula : string consisting of element symbols, numbers and parentheses

	Returns
	-------
	The function returns a list of elemental symbol,number pairs
	"""
	Stoichiometry = []
	m=re.search('((?P<Element>[A-Z][a-z]?)|\((?P<Paren>.*)\))(?P<Number>\d*(\.\d+)?)(?P<Remainder>.*)',Formula)
	#print m.group('Element'),m.group('Number'),m.group('Paren'),m.group('Remainder')
	if len(m.group('Number')) is not 0:
		Number = float(m.group('Number'))
	else:
		Number = 1.0
	if m.group('Element') is not None:
		Z = ConvertElementSymbol(m.group('Element'))
		if Z is not 0:
			Stoichiometry.append([Z,Number])
	elif len(m.group('Paren')) > 0:
		Stoichiometry +=[[x[0],x[1]*Number] for x in ParseChemicalFormula(m.group('Paren'))]
	if len(m.group('Remainder')) is not 0:
		Stoichiometry += ParseChemicalFormula(m.group('Remainder'))
	return Stoichiometry
	
	
	
	
	
	
	
	