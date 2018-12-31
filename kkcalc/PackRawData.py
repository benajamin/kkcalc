#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""
This module takes data from different sources and packages it in a form to be used by and
distributed with the Kramers-Kronig Calculator software package.


Workflow to accomodate:
1. Read data from .nff files and BL files.
2. Combine Henke and BL data sets.
3. Convert to useful format for internal use.
4. Write to json file for distribution.
5. Load data for use by KKcalc.
6. Combine data for different elements as selected by user.
 a) figure out energy values (abscissa).
 b) add coefficients/intensity values in selected proportions
7. Provide combined data in required formats.
 a) list of tuples for plotting.
 b) list of energy ranges, each with corresponding list of polynomial coefficients, (i.e. piecewise polynomial format) for PP KK calculation.

Items 1-4 not usually performed by users. Items 5-7 must be integrated into KKcalc program.

I think that PP format is ideal representation and so will use that as the stored database format. 
Plotable list of tuples can be easily generated from PP energy ranges (maybe with added option to add extra points).


"""

import os, os.path
import scipy, scipy.io, scipy.interpolate
import numpy, math, json

BASEDIR = os.path.dirname(os.path.realpath(__file__))
classical_electron_radius = 2.81794029957951365441605230194258e-15# meters
Plancks_constant = 4.1356673310e-15 # eV*seconds
speed_of_light = 2.99792458e8 # meters per second
Avogadro_constant = 6.02214129e23

Elements_DATA = [line.strip("\r\n").split() for line in open(os.path.join(BASEDIR, 'asf', 'elements.dat'))]
Database = dict()

#################################################################################################################
def LoadData(filename):
	""" Read a standard ascii file and return a list of lists of floats"""
	data = []
	if os.path.isfile(filename):
		for line in open(filename):
			try:
				data.append([float(f) for f in line.split()])
			except ValueError:
				pass
		data = numpy.array(data)
	else:
		print("Error:", filename, "is not a valid file name.")
	if len(data)==0:
		print("Error: no data found in", filename)
		return []
	else:
		return data

def parse_BL_file():
	continue_norm = True # Normalise the Biggs and Lighthill data as the published scattering factors do, rather than as Henke et al says.
	BLfile = {}
	for line in open('original_biggs_file.dat'):
		try:
			values = [float(f) for f in line.split()]
			if values[3] > 10:
				Norm_value = 0 #will calculate actual normalisation value later
				if not continue_norm and values[2] > 10 and values[2] not in [20, 100, 500, 100000]:
					Norm_value = 1
				elif not continue_norm and values[0] == 42 and values[2] > 10 and values[2] not in [100, 500, 100000]:#Mo needs special handling
					#print "Mo seen at", values[0], values[2]
					Norm_value = 1
				values.append(Norm_value)
				if values[2] not in [0.01, 0.1, 0.8, 4, 20, 100, 500, 100000] or (values[0] == 42 and values[2] == 20):
					values.append(1)#this is an absorption edge!
				else:
					values.append(0)#this is not an absorption edge
				BLfile[int(values[0])].append(values)
		except ValueError:
			pass
		except IndexError:
			pass
		except KeyError:
			BLfile[int(values[0])] = [values]
	for elem,coeffs in list(BLfile.items()):
		BLfile[elem] = numpy.array(coeffs)[:,2:]
	return BLfile

def BL_to_ASF(E,coeffs,Atomic_mass):
	"""Biggs and Lighthill offers photoelectric cross-section with the sum of AnE^-n for n=1-4 {E in keV and PECS in cm^2/g}.
			Henke scattering factors related by f2 = PECS*E/(2*r0*h*c)  {E in eV and PECS in cm^2/atom}."""
	return (coeffs[0] + coeffs[1]/(E*0.001) + coeffs[2]/((E*0.001)**2) + coeffs[3]/((E*0.001)**3))*Atomic_mass/(2*Avogadro_constant*classical_electron_radius*Plancks_constant*speed_of_light)*0.1

def Coeffs_to_ASF(E,coeffs):
	"""Calculate Henke scattering factors from polynomial coefficients. {E in eV and PECS in cm^2/atom}."""
	return coeffs[0]*E + coeffs[1] + coeffs[2]/E + coeffs[3]/(E**2) + coeffs[4]/(E**3)
###########################################################################################################
BL_data = parse_BL_file()

#for z, symbol, name, atomic_mass, Henke_file in [Elements_DATA[0]]:
for z, symbol, name, atomic_mass, Henke_file in Elements_DATA:
	print(z, symbol, name, atomic_mass, Henke_file)
	#Get basic metadata
	Element_Database = dict()
	Element_Database['mass'] = float(atomic_mass)
	Element_Database['name'] = name
	Element_Database['symbol'] = symbol
	
	#Get basic data
	print("Load nff data from:", os.path.join(BASEDIR, 'asf', Henke_file))
	asf_RawData = LoadData(os.path.join(BASEDIR, 'asf', Henke_file))
	if min(asf_RawData[1:-1,0]-asf_RawData[0:-2,0])<0:
		print("Warning! Energies in ", Henke_file, "are not in ascending order! (Sorting now..)")
		asf_RawData.sort()
	#print BL_data[int(z)]
	
	#Convert and normalise BL data
	#get normalisation values
	ASF_norm = scipy.interpolate.splev(10000,scipy.interpolate.splrep(asf_RawData[:,0],asf_RawData[:,2],k=1),der=0)
	BL_norm = BL_to_ASF(10000,BL_data[int(z)][0][3:7],float(atomic_mass))
	#print "Norms:", ASF_norm, BL_norm, BL_norm/ASF_norm
	
	temp_E = []
	BL_coefficients = []
	for line in BL_data[int(z)]:
		if float(line[1]) >= 30:
			temp_E.append(float(line[0]))
			BL_coefficients.append(line[2:7]/BL_norm*ASF_norm*[0,1,1000,1000000,1000000000]*float(atomic_mass)/(2*Avogadro_constant*classical_electron_radius*Plancks_constant*speed_of_light)*0.1)
	#store for use in calculation
	C = numpy.array(BL_coefficients)
	#(insert 30000.1 here to use linear section from 30000 to 30000.2 to ensure continuity between data sets)
	X = numpy.array([30.0001]+temp_E[1:])*1000
	
	#Express asf data in PP
	M = (asf_RawData[1:,2]-asf_RawData[0:-1,2])/(asf_RawData[1:,0]-asf_RawData[0:-1,0])
	B = asf_RawData[0:-1,2]-M*asf_RawData[0:-1,0]
	E = asf_RawData[:,0]
	Full_coeffs = numpy.zeros((len(asf_RawData[:,0])-1,5))
	Full_coeffs[:,0] = M
	Full_coeffs[:,1] = B
	#Append B&L data and make sure it is continuous
	E = E[0:-1]
	for i in range(len(X)-1):
		Y1 = Coeffs_to_ASF(X[i]-0.1,Full_coeffs[-1,:])
		Y2 = Coeffs_to_ASF(X[i]+0.1,C[i,:])
		M = (Y2-Y1)/0.2
		B = Y1-M*(X[i]-0.1)
		E = numpy.append(E,[X[i]-0.1,X[i]+0.1])
		Full_coeffs = numpy.append(Full_coeffs,[[M,B,0,0,0]],axis=0)
		Full_coeffs = numpy.append(Full_coeffs,[C[i,:]],axis=0)
	E = numpy.append(E,X[-1])
	
	#convert numpy arrays to nested lists to enable json serialisation with the default converter.
	Element_Database['E'] = E.tolist()
	Element_Database['Im'] = Full_coeffs.tolist()
	Element_Database['Re'] = asf_RawData[:,1].tolist()
	Database[int(z)] = Element_Database

with open('ASF.json','w') as f:
	json.dump(Database,f,indent=1)
