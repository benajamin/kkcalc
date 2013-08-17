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
logger = logging.getLogger(__name__)
import json, os
import numpy, re


# Constants
CLASSICAL_ELECTRON_RADIUS = 2.81794029957951365441605230194258e-15  # meters
PLANCKS_CONSTANT = 4.1356673310e-15  # eV*seconds
SPEED_OF_LIGHT = 2.99792458e8  # meters per second
AVOGADRO_CONSTANT = 6.02214129e23  # no unit
LIST_OF_ELEMENTS = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U']

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
	logger.info("Loading element database")
	with open('ASF.json','r') as f:
		Element_Database = json.load(f)
	for Z in range(1,93):
		Element_Database[str(Z)]['E'] = numpy.array(Element_Database[str(Z)]['E'])
		Element_Database[str(Z)]['Im'] = numpy.array(Element_Database[str(Z)]['Im'])
		Element_Database[str(Z)]['Re'] = numpy.array(Element_Database[str(Z)]['Re'])
	return Element_Database

ELEMENT_DATABASE = load_Element_Database() #spectral data, plus atomic masses

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
	The function returns the magnitude of the real or imaginary
	part of the atomic scattering factors at energy (or energies) `E`.

	"""
	if len(coeffs.shape) == 1:
		return coeffs[0]*E + coeffs[1] + coeffs[2]/E + coeffs[3]/(E**2) + coeffs[4]/(E**3)
	else:
		return coeffs[:,0]*E + coeffs[:,1] + coeffs[:,2]/E + coeffs[:,3]/(E**2) + coeffs[:,4]/(E**3)

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
	try:
		AtomicNumber = LIST_OF_ELEMENTS.index(SymbolString)+1
	except ValueError:
		logger.warn("\""+SymbolString+"\" is not a known element!")
		AtomicNumber = 0
	return AtomicNumber

def ParseChemicalFormula(Formula,recursion_flag=False):
	"""Parse a chemical formula string to obtain a stoichiometry.

	Parameters
	----------
	Formula : string consisting of element symbols, numbers and parentheses

	Returns
	-------
	The function returns a list of elemental symbol,number pairs
	"""
	if not recursion_flag:
		logger.info("Parsing '"+Formula+"' as a chemical formula")
	Stoichiometry = []
	m=re.search('((?P<Element>[A-Z][a-z]?)|\((?P<Paren>.*)\))(?P<Number>\d*(\.\d+)?)(?P<Remainder>.*)',Formula)
	if len(m.group('Number')) is not 0:
		Number = float(m.group('Number'))
	else:
		Number = 1.0
	if m.group('Element') is not None:
		Z = ConvertElementSymbol(m.group('Element'))
		if Z is not 0:
			Stoichiometry.append([Z,Number])
	elif len(m.group('Paren')) > 0:
		Stoichiometry +=[[x[0],x[1]*Number] for x in ParseChemicalFormula(m.group('Paren'),recursion_flag=True)]
	if len(m.group('Remainder')) is not 0:
		Stoichiometry += ParseChemicalFormula(m.group('Remainder'),recursion_flag=True)
	return Stoichiometry
	
	
def load_data(filename, load_options):
	"""Read a standard ASCII file and return a list of lists of floats.
	
	Parameters
	----------
	filename : string path to data file
	load_options : dictionary of optional settings

	Returns
	-------
	The function returns a numpy array with two columns: Photon energy and Imaginary scattering factor values
"""
	logger.info("Load data from file")
	data = []
	if os.path.isfile(filename):
		for line in open(filename):
			try:
				data.append([float(f) for f in line.split()])
			except ValueError:
				pass
		data = numpy.array(data)
	else:
		logger.error(filename+" is not a valid file name.")
	if len(data)==0:
		logger.error("no data found in "+filename)
		return None
	else:
		if isinstance(load_options,dict) and load_options.has_key('E_column'):
			E_column = int(load_options['E_column'])
		else:
			E_column = 0
		if isinstance(load_options,dict) and load_options.has_key('data_column'):
			data_column = int(load_options['data_column'])
		else:
			data_column = data.shape[1]-1
		if isinstance(load_options,dict) and load_options.has_key('convert_from'):
			data_type = load_options['convert_from']
		else:
			data_type = 'photoabsorption'
		if data_type == 'photoabsorption':
			converted_data = convert_NEXAFS_to_ASF(data[:,[E_column,data_column]])
		if data_type == 'beta':
			if isinstance(load_options,dict) and load_options.has_key('density'):
				density = load_options['density']
			else:
				density = 1.0
				logger.warn("Assuming material density of 1 gcm^-3 for conversion from Beta to ASF.")
			converted_data = convert_Beta_to_ASF(data[:,[E_column,data_column]],density)
		if data_type == 'asf':
			converted_data = data[:,[E_column,data_column]]
		return converted_data
	
	
def calculate_asf(Stoichiometry):
	logger.info("Calculate material scattering factor data from the given stoichiometry")
	if len(Stoichiometry) is 0:
		logger.error("No elements described by input.")
		return None
	else:
		# get unique energy points
		total_E = numpy.array([])
		for element,n in Stoichiometry:
			total_E = numpy.concatenate((total_E, ELEMENT_DATABASE[str(element)]['E']))
		total_E = numpy.unique(total_E)
		# add weighted asf data sets for KK calculation
		total_Im_coeffs = numpy.zeros((len(total_E)-1, 5))
		counters = numpy.zeros((len(Stoichiometry)))
		for i,E in enumerate(total_E[1:]):
			sum_Im_coeffs = 0
			for j in range(len(counters)):
				sum_Im_coeffs += Stoichiometry[j][1]*ELEMENT_DATABASE[str(Stoichiometry[j][0])]['Im'][counters[j],:]
				counters[j] += ELEMENT_DATABASE[str(Stoichiometry[j][0])]['E'][counters[j]+1] == E
			total_Im_coeffs[i,:] = sum_Im_coeffs
		return total_E, total_Im_coeffs

def merge_spectra(NearEdge_Data, ASF_E, ASF_Data, merge_points=None, add_background=False, fix_distortions=False):
	logger.info("Merge near-edge data with wide-range scattering factor data")
	if merge_points is None:
		merge_points = NearEdge_Data[[0,-1],0]
	NearEdge_merge_ind = [numpy.where(NearEdge_Data[:,0] > merge_points[0])[0][0], numpy.where(NearEdge_Data[:,0] < merge_points[1])[0][-1]]
	asf_merge_ind = [numpy.where(ASF_E > merge_points[0])[0][0]-1, numpy.where(ASF_E > merge_points[1])[0][0]-1]
	NearEdge_merge_values = numpy.interp(merge_points, NearEdge_Data[:,0], NearEdge_Data[:,1])
	NearEdge_Data = numpy.vstack(((merge_points[0],NearEdge_merge_values[0]),NearEdge_Data[NearEdge_merge_ind[0]:NearEdge_merge_ind[1],:],(merge_points[1],NearEdge_merge_values[1])))
	asf_merge_values = [coeffs_to_ASF(merge_points[0], ASF_Data[asf_merge_ind[0],:]), coeffs_to_ASF(merge_points[1], ASF_Data[asf_merge_ind[1],:])]
	if add_background:
		logger.error("Not implemented!")
		#get pre-edge region
		#extrapolate background
	else:
		scale = (asf_merge_values[1]-asf_merge_values[0])/(NearEdge_merge_values[1]-NearEdge_merge_values[0])
		NearEdge_Data[:,1] = ((NearEdge_Data[:, 1]-NearEdge_merge_values[0])*scale)+asf_merge_values[0]
	if fix_distortions:
		logger.error("Not implemented!")
		#use fitting function to estimate curvature
	
	#Now convert point values to coefficients
	NearEdge_Coeffs = numpy.zeros((len(NearEdge_Data),5))
	M = (NearEdge_Data[1:,1]-NearEdge_Data[:-1,1])/(NearEdge_Data[1:,0]-NearEdge_Data[:-1,0])
	NearEdge_Coeffs[:-1,0] = M
	NearEdge_Coeffs[:-1,1] = NearEdge_Data[:-1,1]-M*NearEdge_Data[:-1,0]
	NearEdge_Coeffs[-1,:] = ASF_Data[asf_merge_ind[1],:]
	#Paste data sections together
	Full_E = numpy.hstack((ASF_E[0:asf_merge_ind[0]],NearEdge_Data[:,0],ASF_E[asf_merge_ind[1]+1:]))
	Full_Coeffs = numpy.vstack((ASF_Data[0:asf_merge_ind[0],:],NearEdge_Coeffs,ASF_Data[asf_merge_ind[1]+1:,:]))
	
	return Full_E, Full_Coeffs
	
	
	
	
	
	