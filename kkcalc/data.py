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
import json, os, warnings
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
	with open(os.path.join(os.path.dirname(__file__),'ASF.json'),'r') as f:
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
		if len(E) == coeffs.shape[0]+1:
			coeffs = numpy.vstack((coeffs,coeffs[-1,:])) #Use last defined polynomial to calculate the ASF values at both second last and the last energies.
		return coeffs[:,0]*E + coeffs[:,1] + coeffs[:,2]/E + coeffs[:,3]/(E**2) + coeffs[:,4]/(E**3)

def convert_Beta_to_ASF(raw_data, density=None, formula_mass=None, stoichiometry=None, number_density=None, reverse=False):
	"""Convert Beta data (index of refraction) to atomic scattering
	factors (ASF).

	The Beta value is the imaginary part of the index of refraction.
	This represents the absorption.

	Parameters
	----------
	raw_data : two-dimensional `numpy.array` of `float`
		The array consists of two columns: Energy and magnitude.
	density : float
		material density in grams per millilitre.
	formula_mass : float
		atomic mass sum of the materials chemical formula (molecular mass).
	stoichiometry : a list of elemental symbol,number pairs
		description of the combination of elements composing the material.
	number_density : float
		material density in atoms per millilitre.
	reverse : boolean
		flag to indicate the reverse conversion.
	
	Note that the density in millilitres makes for a factor of a million in the equations.

	Returns
	-------
	The function returns a `numpy.array` of atomic scattering factors.
	They are made up of the energy and the magnitude of the imaginary
	part of the atomic scattering factors.

	"""
	if number_density is None:
		if density is None and formula_mass is None and stoichiometry is None:
			logger.info("Material number density required for conversions involving Beta.")
		if density is None:
			logger.info("Assuming a material mass density of 1 g/ml.")
			density = 1.0
		if formula_mass is None:
			if stoichiometry is not None:
				formula_mass = calculate_FormulaMass(stoichiometry)
			else:
				logger.info("Assuming a formula mass 100.")
				formula_mass = 100.
		number_density = density*AVOGADRO_CONSTANT/formula_mass
	if reverse: #ASF to Beta
		raw_Im = raw_data[:, :2].copy()
		raw_Im[:, 1] = 1000000.*raw_data[:, 1]*(number_density*CLASSICAL_ELECTRON_RADIUS*(PLANCKS_CONSTANT*SPEED_OF_LIGHT)**2)/(2*numpy.pi*raw_data[:, 0]**2)
	else: #Beta to ASF
		raw_Im = raw_data[:, :2].copy()
		raw_Im[:, 1] = 0.000001*raw_data[:, 1]*2*numpy.pi*raw_data[:, 0]**2/(number_density*CLASSICAL_ELECTRON_RADIUS*(PLANCKS_CONSTANT*SPEED_OF_LIGHT)**2)
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


def convert_NEXAFS_to_ASF(raw_data, reverse=False):
	"""Convert NEXAFS photoabsorption data to atomic scattering factors (ASF).

	Parameters
	----------
	raw_data : two-dimensional `numpy.array` of `float`
		The array consists of two columns: Energy and magnitude.
	reverse : boolean
		flag to indicate the reverse conversion

	Returns
	-------
	The function returns a `numpy.array` of atomic scattering factors.
	They are made up of the energy and the magnitude of the imaginary
	part of the atomic scattering factors.

	"""
	if reverse:
		raw_Im = raw_data[:, :2].copy()
		raw_Im[:, 1] = (2*CLASSICAL_ELECTRON_RADIUS*PLANCKS_CONSTANT*SPEED_OF_LIGHT)*raw_data[:, 1]/raw_data[:, 0]
	else:
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
	
	TODO: return understood formula (i.e. input strign with errors removed)
	"""
	if not recursion_flag:
		logger.info("Parsing '"+Formula+"' as a chemical formula")
	Stoichiometry = []
	m=re.search('((?P<Element>[A-Z][a-z]?)|\((?P<Paren>.*)\))(?P<Number>\d*(\.\d+)?)(?P<Remainder>.*)',Formula)
	if len(m.group('Number')) != 0:
		Number = float(m.group('Number'))
	else:
		Number = 1.0
	if m.group('Element') is not None:
		Z = ConvertElementSymbol(m.group('Element'))
		if Z != 0:
			Stoichiometry.append([Z,Number])
	elif len(m.group('Paren')) > 0:
		Stoichiometry +=[[x[0],x[1]*Number] for x in ParseChemicalFormula(m.group('Paren'),recursion_flag=True)]
	if len(m.group('Remainder')) != 0:
		Stoichiometry += ParseChemicalFormula(m.group('Remainder'),recursion_flag=True)
	return Stoichiometry
	
def calculate_FormulaMass(Stoichiometry):
	"""Sum atomic masses
	Parameters
	-------
	Stoichiometry : a list of elemental symbol,number pairs

	Returns
	-------
	The function returns the fomula mass as a float
	"""
	FormulaMass = 0
	for element, number in Stoichiometry:
		FormulaMass += number*ELEMENT_DATABASE[str(element)]['mass']
	return FormulaMass
	
def load_data(filename, load_options=None):
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
		if isinstance(load_options,dict) and 'E_column' in load_options:
			E_column = int(load_options['E_column'])
		else:
			E_column = 0
		if isinstance(load_options,dict) and 'data_column' in load_options:
			data_column = int(load_options['data_column'])
		else:
			data_column = data.shape[1]-1   
		return data[:,[E_column,data_column]]
	
def export_data(filename, data, header_info=None,convert_to=None):
	"""Write spectral data into an ASCII file.
	
	Parameters
	----------
	filename : string path to data file
	data : array of floats, arranged in columns, with energy in the first.
	save_options : dictionary of optional settings

	Returns
	-------
	None
	"""
	logger.info("Export data to file")
	if data is not None:
		column_headings = '# E(eV)\tf1\tf2\n'
		data_type = "Scattering factors"
		if convert_to == 'photoabsorption':
			data = convert_data(data[:,[0,2]],'ASF',convert_to, Density=header_info['Density'])
			column_headings = '# E(eV)\tphotoabsorption\n'
			data_type = "Photoabsorption spectrum"
		elif convert_to == 'refractive_index':
			data[:,1] = convert_data(data[:,[0,1]],'ASF',convert_to, Density=header_info['Density'], Formula_Mass=header_info["Formula Mass"])[:,1]
			data[:,2] = convert_data(data[:,[0,2]],'ASF',convert_to, Density=header_info['Density'], Formula_Mass=header_info["Formula Mass"])[:,1]
			column_headings = '# E(eV)\tDelta\tBeta\n'
			data_type = "Refractive index"
		with open(filename, "w") as outfile:
			if header_info is not None:
				for key in header_info:
					outfile.write('# '+key+' = '+str(header_info[key])+'\n')
			outfile.write(column_headings)
			numpy.savetxt(outfile,data,fmt="%7g",delimiter="\t")
		logger.info(data_type+" for "+header_info["Molecular Formula"]+" saved to "+filename)
	else:
		logger.info("Nothing to save.")









def convert_data(Data, FromType, ToType, Density=None, Formula_Mass=None):
	"""Switchyard function for converting between data types.
	
	Parameters
	----------
	Data : a numpy array with two columns: Photon energy and absorption data values
	FromType : string describing input data type
	ToType : string describing desired output data type

	Returns
	-------
	The function returns a numpy array with two columns: Photon energy and absorption data values
	"""
	logger.info("Convert data from "+FromType+" to "+ToType+".")
	if FromType.lower() in ['photoabsorption', 'nexafs', 'xanes']:
		converted_data = convert_NEXAFS_to_ASF(Data)
	elif FromType.lower() in ['beta']:
		converted_data = convert_Beta_to_ASF(Data, density=Density)
	else: #if already ASF
		converted_data = Data.copy()
	#data should now be in terms of ASF
	if ToType.lower() in ['photoabsorption', 'nexafs', 'xanes']:
		converted_data = convert_NEXAFS_to_ASF(converted_data, reverse=True)
	elif ToType.lower() in ['refractive_index','beta']:
		converted_data = convert_Beta_to_ASF(converted_data, density=Density, formula_mass=Formula_Mass, reverse=True)
	return converted_data
	
	
def calculate_asf(Stoichiometry):
	"""Sum scattering factor data for a given chemical stoichiometry.
	
	Parameters
	----------
	Stoichiometry : a list of elemental symbol,number pairs
	
	Returns
	-------
	total_E: 1D numpy array listing the starting photon energies of the segments that the spectrum is broken up into.
	total_Im_coeffs: nx5 numpy array in which each row lists the polynomial coefficients describing the shape of the spectrum in that segment.
	"""
	logger.info("Calculate material scattering factor data from the given stoichiometry")
	if len(Stoichiometry) == 0:
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
		counters = numpy.zeros((len(Stoichiometry)),dtype=int)
		for i,E in enumerate(total_E[1:]):
			sum_Im_coeffs = 0
			for j in range(len(counters)):
				sum_Im_coeffs += Stoichiometry[j][1]*ELEMENT_DATABASE[str(Stoichiometry[j][0])]['Im'][counters[j],:]
				counters[j] += ELEMENT_DATABASE[str(Stoichiometry[j][0])]['E'][counters[j]+1] == E
			total_Im_coeffs[i,:] = sum_Im_coeffs
		return total_E, total_Im_coeffs

def merge_spectra(NearEdge_Data, ASF_E, ASF_Data, merge_points=None, add_background=False, fix_distortions=False, plotting_extras=False):
	"""Normalise the user-provided, near-edge data to the scattering factor data and combine them.
	
	Parameters
	----------
	NearEdge_Data : a numpy array with two columns: Photon energy and absorption data values.
	ASF_E : 1D numpy array listing the starting photon energies for each spectrum segment.
	ASF_Data: nx5 numpy array in which each row lists the polynomial coefficients describing the shape of the spectrum in that segment.
	merge_points : a pair of photon energies indicating the data range of the NearEdge_Data to be used.
	add_background : boolean flag for adding a background function to the provided near-edge data.
	fix_distortions : boolean flag for removing erroneous slope from the provided near-edge data.
	plotting_extras : boolean flag for providing plottable feedback data that isn't required for actual calculations.

	Returns
	-------
	Full_E : 1D numpy array listing the updated starting photon energies for each spectrum segment.
	Full_Coeffs : nx5 numpy array in which each row lists the polynomial coefficients describing the shape of the spectrum in that segment.
	plot_scaled_NearEdge_Data (optional) : an updated verion of NearEdge_Data that has been scaled to match the scattering factor data.
	splice_points (optional) : a list of two pairs of photon energy-magnitude values to represent the merge_points in a plot.
	"""
	logger.info("Merge near-edge data with wide-range scattering factor data")
	if merge_points is None:
		merge_points = NearEdge_Data[[0,-1],0]
	NearEdge_merge_ind = [numpy.where(NearEdge_Data[:,0] > merge_points[0])[0][0], numpy.where(NearEdge_Data[:,0] < merge_points[1])[0][-1]]
	NearEdge_merge_values = numpy.interp(merge_points, NearEdge_Data[:,0], NearEdge_Data[:,1])
	if ASF_Data is not None:
		asf_merge_ind = [numpy.where(ASF_E > merge_points[0])[0][0]-1, numpy.where(ASF_E > merge_points[1])[0][0]-1]
		asf_merge_values = [coeffs_to_ASF(merge_points[0], ASF_Data[asf_merge_ind[0],:]), coeffs_to_ASF(merge_points[1], ASF_Data[asf_merge_ind[1],:])]
		if add_background:
			logger.info("Add background")
			logger.error("Not implemented!")
			#get pre-edge region
			#extrapolate background
			scale = (asf_merge_values[1]-asf_merge_values[0])/(NearEdge_merge_values[1]-NearEdge_merge_values[0])
			scaled_NearEdge_Data = numpy.vstack((NearEdge_Data[:,0],((NearEdge_Data[:, 1]-NearEdge_merge_values[0])*scale)+asf_merge_values[0])).T
		else:# don't add background
			scale = (asf_merge_values[1]-asf_merge_values[0])/(NearEdge_merge_values[1]-NearEdge_merge_values[0])
			scaled_NearEdge_Data = numpy.vstack((NearEdge_Data[:,0],((NearEdge_Data[:, 1]-NearEdge_merge_values[0])*scale)+asf_merge_values[0])).T
		try:
			import scipy.optimize
			SCIPY_FLAG = True
		except ImportError:
			SCIPY_FLAG = False
			logger.info('Failed to import the scipy.optimize module - disabling the \'fix distortions\' option.')
		if SCIPY_FLAG and fix_distortions:
			logger.info("Fix distortions")
			import scipy.optimize
			ASF_fitY = 0.0*NearEdge_Data[:, 0]
			for i, E in enumerate(NearEdge_Data[:,0]):
				ASF_fitY[i] = coeffs_to_ASF(E, ASF_Data[numpy.where(ASF_E > E)[0][0]-1])
			fitfunc = lambda p, x, y, asf_mv, asf: ((y-p*x)-(y[0]-p*x[0]))/((y[-1]-p*x[-1])-(y[0]-p*x[0]))*(asf_mv[1]-asf_mv[0])+asf_mv[0] - asf
			p0 = -(NearEdge_merge_values[1]-NearEdge_merge_values[0])/((asf_merge_values[1]-asf_merge_values[0])*NearEdge_Data[0,0])
			logger.debug("Fix distortions - start fit with p0 ="+str(p0))
			p1, success = scipy.optimize.leastsq(fitfunc, p0, args=(NearEdge_Data[:, 0], NearEdge_Data[:, 1], asf_merge_values, ASF_fitY))
			logger.debug("Fix distortions - complete fit with p1 ="+str(p1[0]))
			NearEdge_fitY = asf_merge_values[0]+((NearEdge_Data[:,1]-p1[0]*NearEdge_Data[:,0])-(NearEdge_Data[0,1]-p1[0]*NearEdge_Data[0,0]))*(asf_merge_values[1]-asf_merge_values[0])/((NearEdge_Data[-1,1]-p1[0]*NearEdge_Data[-1,0])-(NearEdge_Data[0,1]-p1[0]*NearEdge_Data[0,0]))
			scaled_NearEdge_Data = numpy.vstack((NearEdge_Data[:,0],NearEdge_fitY)).T
	else:
		scaled_NearEdge_Data = NearEdge_Data.copy()
		asf_merge_values = NearEdge_merge_values.copy()
	plot_scaled_NearEdge_Data = scaled_NearEdge_Data.copy()
	scaled_NearEdge_Data = numpy.vstack(((merge_points[0],asf_merge_values[0]),scaled_NearEdge_Data[NearEdge_merge_ind[0]:NearEdge_merge_ind[1]+1,:],(merge_points[1],asf_merge_values[1])))

	#Now convert point values to coefficients
	NearEdge_Coeffs = numpy.zeros((len(scaled_NearEdge_Data),5))
	M = (scaled_NearEdge_Data[1:,1]-scaled_NearEdge_Data[:-1,1])/(scaled_NearEdge_Data[1:,0]-scaled_NearEdge_Data[:-1,0])
	NearEdge_Coeffs[:-1,0] = M
	NearEdge_Coeffs[:-1,1] = scaled_NearEdge_Data[:-1,1]-M*scaled_NearEdge_Data[:-1,0]
	if ASF_Data is None:
		NearEdge_Coeffs = NearEdge_Coeffs[0:-1,:]
		#NearEdge_Coeffs[-1,:] = NearEdge_Coeffs[-2,:]
		#Paste data sections together
		Full_E = scaled_NearEdge_Data[:,0]
		Full_Coeffs = NearEdge_Coeffs
		splice_points = [0,len(Full_E)-1]
	else:
		NearEdge_Coeffs[-1,:] = ASF_Data[asf_merge_ind[1],:]
		#Paste data sections together
		Full_E = numpy.hstack((ASF_E[0:asf_merge_ind[0]+1],scaled_NearEdge_Data[:,0],ASF_E[asf_merge_ind[1]+1:]))
		Full_Coeffs = numpy.vstack((ASF_Data[0:asf_merge_ind[0]+1,:],NearEdge_Coeffs,ASF_Data[asf_merge_ind[1]+1:,:]))
		splice_points = [asf_merge_ind[0]+1, asf_merge_ind[0]+len(scaled_NearEdge_Data[:,0])]
	if plotting_extras:
		return Full_E, Full_Coeffs, plot_scaled_NearEdge_Data, splice_points
	else:
		return Full_E, Full_Coeffs
	
	
def coeffs_to_linear(E, coeffs, threshold):
	"""Convert a curved data set (described by polynomial coefficients) to an approximate data set composed of linear segments.
	   This should be useful for plotting
	
	Parameters
	----------
	E : 1D numpy array listing the starting photon energies for each spectrum segment.
	coeffs: nx5 numpy array in which each row lists the polynomial coefficients describing the shape of the spectrum in that segment.
	threshold: scalar value indicating the maximum amount of error to be tolerated.
	
	Returns
	-------
	linear_E: 1D numpy array listing the starting photon energies of the segments that the spectrum is broken up into.
	linear_Vals: 1D numpy array listing the  intensity values corresponding to the energies listed in linear_E.
	"""
	logger.info("Linearise data for plotting.")
	curvature = 2*coeffs[:,2]/(E[0:-1]**3) + 6*coeffs[:,3]/(E[0:-1]**4) + 12*coeffs[:,4]/(E[0:-1]**5)
	linear_E = numpy.array([])
	linear_Vals = numpy.array([])
	for i in range(len(E)-1):
		linear_E = numpy.append(linear_E, E[i])
		linear_Vals = numpy.append(linear_Vals, coeffs_to_ASF(E[i],coeffs[i]))
		if not numpy.all(coeffs[i,2::] == [0,0,0]):
			new_points = split_segment_recursively(E[i], E[i+1], coeffs[i], threshold)
			#print new_points
			if new_points is not None:
				for p in new_points:
					linear_E = numpy.append(linear_E, p[0])
					linear_Vals = numpy.append(linear_Vals, p[1])
	linear_E = numpy.append(linear_E, E[-1])
	linear_Vals = numpy.append(linear_Vals, coeffs_to_ASF(E[-1],coeffs[-1]))
	return linear_E, linear_Vals

def split_segment_recursively(E1, E2, coeffs, threshold):
	"""Recursively use split_segment() to divide a spectrum segment until all points along the piecewise linear representation
	is within 'threshold' of the curve defined by 'coeffs'.
	
	Parameters
	----------
	E1 : Energy value at the start of the segment.
	E2 : Energy value at the end of the segment.
	coeffs: 1x5 numpy array listing the polynomial coefficients describing the shape of the spectrum in that segment.
	threshold: scalar value indicating the maximum amount of error to be tolerated.
	
	Returns
	-------
	The function returns a list of energy, magnitude tuples representing an approximation of the curve described by 'coeffs'.
	"""
	new_point = split_segment(E1, E2, coeffs)
	if new_point is not None and new_point[2] > threshold:
		#print new_point
		point_list = []
		new_point_2 = split_segment_recursively(E1, new_point[0], coeffs, threshold)#recursively check first half
		if new_point_2 is not None:
			point_list = point_list + new_point_2
		point_list.append(new_point[0:2])
		new_point_3 = split_segment_recursively(new_point[0], E2, coeffs, threshold)#recursively check second half
		if new_point_3 is not None:
			point_list = point_list + new_point_3
		return point_list
	else:
		return None

def split_segment(E1, E2, coeffs):
	"""Calculate which point in a spectrum segment shows the greatest error when represented by a straight line.
	
	Parameters
	----------
	E1 : Energy value at the start of the segment.
	E2 : Energy value at the end of the segment.
	coeffs: 1x5 numpy array listing the polynomial coefficients describing the shape of the spectrum in that segment.
	
	Returns
	-------
	The function returns the energy, magnitude and error value corresponding to the point at which the greatest error occurs.
	"""
	vals = [coeffs_to_ASF(E1, coeffs), coeffs_to_ASF(E2, coeffs)]
	m = (vals[1] - vals[0])/(E2 - E1)
	a = m-coeffs[0]
	b = 0.0
	c = coeffs[2]
	d = 2*coeffs[3]
	e = 3*coeffs[4]
	
	# see http://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots
	D0 = c**2-3*b*d+12*a*e
	D1 = 2*c**3-9*b*c*d+27*b**2*e+27*a*d**2-72*a*c*e
	p = (8*a*c-3*b**2)/(8.0*a**2)
	q = (b**3-4*a*b*c+8*a**2*d)/(8.0*a**3)
	Q_cubed = (D1+(D1**2-4*D0**3)**0.5)*0.5
	Q = numpy.sign(Q_cubed)*abs(Q_cubed)**(1/3.0)
	S = 0.5*(-3*p/2.0+1/(3.0*a)*(Q+D0/Q))**0.5
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		#This line will generate warnings when failing to calculate complex roots properly, but we only care about real roots anyway. 
		roots = numpy.array([-b/(4.0*a)-S+0.5*(-4*S**2-2*p+q/S)**0.5, -b/(4.0*a)-S-0.5*(-4*S**2-2*p+q/S)**0.5, -b/(4.0*a)+S+0.5*(-4*S**2-2*p-q/S)**0.5, -b/(4.0*a)+S-0.5*(-4*S**2-2*p-q/S)**0.5])
	good_roots = []
	for r in roots:
		if numpy.isfinite(r) & (E1<r<E2):
			good_roots.append(r)
	if len(good_roots) == 0:
		return None
	elif len(good_roots) == 1:
		new_E = good_roots[0]
		distance = abs(vals[0]+(new_E-E1)*m - coeffs_to_ASF(new_E, coeffs))
	else:
		distance = abs(vals[0]+(good_roots-E1)*m - coeffs_to_ASF(good_roots, coeffs))
		best_root_ind = numpy.argmin(distance)
		new_E = good_roots[best_root_ind]
		distance = distance[best_root_ind]
	return new_E, coeffs_to_ASF(new_E, coeffs), distance


	
