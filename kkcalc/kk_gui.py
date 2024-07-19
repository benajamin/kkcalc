#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""This module implements a GUI using the wxPython toolkit."""
import logging
logger = logging.getLogger(__name__)
if __name__ == '__main__':
	import sys
	logging.basicConfig(level=logging.DEBUG)
	logging.StreamHandler(stream=sys.stdout)

import wx
import wx.lib.plot as plot
import numpy
import os
from . import kk, data
from .__init__ import __version__ as version

try:
	import scipy.optimize
	SCIPY_FLAG = True
except ImportError:
	SCIPY_FLAG = False
	logger.info('Failed to import the scipy.optimize module - disabling the \'fix distortions\' checkbox.')


class MyFrame(wx.Frame):

	def __init__(self):
		wx.Frame.__init__(self, None, wx.ID_ANY, f"Kramers-Kronig Calculator {version}", size=(500, 800))

		# Initialise variables
		self.dirname = ''
		self.raw_file = None
		self.total_asf = None
		self.total_Im_coeffs = None
		self.merged_Im = None
		self.nexafs_CutOut = []
		self.MolecularMass = 1
		self.asf_bg = None
		#New set of variables to initialise. All those above might want to be removed.
		self.ChemicalFormula = None
		self.Stoichiometry = None
		self.Relativistic_Correction = None
		self.NearEdgeData = None
		self.splice_ind = None
		self.ASF_E = None
		self.ASF_Data = None
		self.Full_E = None
		self.Imaginary_Spectrum = None
		self.KK_Real_Spectrum = None


		# Setting up the menus.
		filemenu = wx.Menu()
		filemenu.Append(wx.ID_OPEN, "L&oad", " Load photoabsorption data from file")
		filemenu.AppendSeparator()
		filemenu.Append(wx.ID_SAVE, "&Save", " Export results to file")
		exportmenu = wx.Menu()
		exportmenu.Append(201,"Photoabsorption", " Export X-ray absorption data")
		exportmenu.Append(202,"Refractive Index", " Export beta and delta")
		filemenu.AppendSubMenu(exportmenu,"Export")   # Adding the "exportmenu" to the filemenu
		filemenu.AppendSeparator()
		filemenu.Append(wx.ID_EXIT, "E&xit", " Terminate the program")
		helpmenu = wx.Menu()
		helpmenu.Append(wx.ID_HELP, "&Help", " How to use this program")
		helpmenu.Append(205, "&Article", " Publication describing how this calculation works")
		helpmenu.AppendSeparator()
		helpmenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")  # Adding the "filemenu" to the MenuBar
		menuBar.Append(helpmenu, "&Help")  # Adding the "helpmenu" to the MenuBar
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		self.Bind(wx.EVT_MENU, self.OnOpen, id=wx.ID_OPEN)
		self.Bind(wx.EVT_MENU, self.OnSave, id=wx.ID_SAVE)
		self.Bind(wx.EVT_MENU, self.OnSave, id=201)   # will set convert_to="photoabsorption" when ID is recognised
		self.Bind(wx.EVT_MENU, self.OnSave, id=202)   # will set convert_to="refractive_index" when ID is recognised
		self.Bind(wx.EVT_MENU, self.OnExit, id=wx.ID_EXIT)
		self.Bind(wx.EVT_MENU, self.OnAbout, id=wx.ID_ABOUT)
		self.Bind(wx.EVT_MENU, self.OnArticle, id=205)
		self.Bind(wx.EVT_MENU, self.OnHelp, id=wx.ID_HELP)


		Sizer1 = wx.BoxSizer(wx.HORIZONTAL)  # create outer sizer
		SizerL = wx.BoxSizer(wx.VERTICAL)  # create left-hand sizer for controls
		SizerR = wx.BoxSizer(wx.VERTICAL)  # create right-hand sizer for plots
		############################Data box
		DataBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Near-Edge Data"), wx.VERTICAL)
		self.FileText = wx.StaticText(self, -1, "File: (None)")
		DataBox.Add(self.FileText, 1, wx.GROW)
		DataTypeLabel = wx.StaticText(self, -1, "Data Type: ")
		self.DataTypeCombo = wx.ComboBox(self, -1, value='Photoabsorption', style=wx.CB_READONLY)
		self.DataTypeCombo.Append('Photoabsorption')
		self.DataTypeCombo.Append('Beta')
		self.DataTypeCombo.Append('Scattering Factor')
		self.DataTypeCombo.Bind(wx.EVT_COMBOBOX, self.MergeAdd_check)
		DataTypeSizer = wx.BoxSizer(wx.HORIZONTAL)
		DataTypeSizer.Add(DataTypeLabel)
		DataTypeSizer.Add(self.DataTypeCombo, 2, wx.GROW)
		DataBox.Add(DataTypeSizer, 1, wx.GROW)
		SpliceSizer = wx.BoxSizer(wx.HORIZONTAL)
		self.SpliceText1 = wx.TextCtrl(self, -1, "Start", style=wx.TE_PROCESS_ENTER)
		self.SpliceText1.Bind(wx.EVT_KILL_FOCUS, self.Splice_Text_check)
		self.SpliceText1.Bind(wx.EVT_TEXT_ENTER, self.Splice_Text_check)
		SpliceSizer.Add(self.SpliceText1, 1)
		self.SpliceText2 = wx.TextCtrl(self, -1, "End", style=wx.TE_PROCESS_ENTER)
		self.SpliceText2.Bind(wx.EVT_KILL_FOCUS, self.Splice_Text_check)
		self.SpliceText2.Bind(wx.EVT_TEXT_ENTER, self.Splice_Text_check)
		SpliceSizer.Add(self.SpliceText2, 1)
		DataBox.Add(SpliceSizer, 1, wx.GROW)

		# Background_CloseSizer = wx.BoxSizer(wx.HORIZONTAL)
#		self.InvertDataCheckBox = wx.CheckBox(self, -1, "Invert Data")
#		self.InvertDataCheckBox.Bind(wx.EVT_CHECKBOX, self.Splice_Text_check)
#		DataBox.Add(self.InvertDataCheckBox, 0)
		self.AddBackgroundCheckBox = wx.CheckBox(self, -1, "Add background")
		self.AddBackgroundCheckBox.Bind(wx.EVT_CHECKBOX, self.Splice_Text_check)
		self.AddBackgroundCheckBox.Disable()
		self.AddBackgroundCheckBox.SetToolTip(wx.ToolTip("Not implemented"))
		DataBox.Add(self.AddBackgroundCheckBox, 0)
		self.FixDistortionsCheckBox = wx.CheckBox(self, -1, "Fix distortions")
		self.FixDistortionsCheckBox.Bind(wx.EVT_CHECKBOX, self.Splice_Text_check)
		if not SCIPY_FLAG:
			self.FixDistortionsCheckBox.Disable()
			self.FixDistortionsCheckBox.SetToolTip(wx.ToolTip("Install the SciPy module to use this feature"))
		DataBox.Add(self.FixDistortionsCheckBox, 0)
		# Background_CloseSizer.Add(self.AddBackgroundCheckBox, 0)
		# self.AddBackgroundCheckBox.Bind(wx.EVT_CHECKBOX, self.MergeAdd_check)
		# Background_CloseSizer.AddStretchSpacer(1)
		# self.CloseFile = wx.Button(self, -1, "X", style= wx.BU_EXACTFIT)
		# Background_CloseSizer.Add(self.CloseFile, 0)
		# DataBox.Add(Background_CloseSizer, 1, wx.GROW)

		############################Material box
		self.MaterialBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Material"), wx.VERTICAL)
		DensitySizer = wx.BoxSizer(wx.HORIZONTAL)
		DensitySizer.Add(wx.StaticText(self, -1, "Density: "))
		self.DensityText = wx.TextCtrl(self, -1, "1", style=wx.TE_PROCESS_ENTER)
		self.DensityText.Bind(wx.EVT_KILL_FOCUS, self.Splice_Text_check)
		self.DensityText.Bind(wx.EVT_TEXT_ENTER, self.Splice_Text_check)
		DensitySizer.Add(self.DensityText, 1)
		DensitySizer.Add(wx.StaticText(self, -1, " g/ml"))
		self.MaterialBox.Add(DensitySizer, 0)
		
		StoichiometrySizer = wx.BoxSizer(wx.HORIZONTAL)
		StoichiometrySizer.Add(wx.StaticText(self, -1, "Stoichiometry: "))
		self.StoichiometryText = wx.TextCtrl(self, -1, "", style=wx.TE_PROCESS_ENTER)
		self.StoichiometryText.Bind(wx.EVT_KILL_FOCUS, self.Stoichiometry_Text_check)
		self.StoichiometryText.Bind(wx.EVT_TEXT_ENTER, self.Stoichiometry_Text_check)
		StoichiometrySizer.Add(self.StoichiometryText, 1)
		self.MaterialBox.Add(StoichiometrySizer, 0)

		############################Calc box
		CalcBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Calculation"), wx.VERTICAL)
		CalcButton = wx.Button(self, -1, "Calculate")
		CalcBox.Add(CalcButton, 1, wx.GROW)
		CalcButton.Bind(wx.EVT_BUTTON, self.calculate)



		SizerL.Add(DataBox, 0, wx.GROW)
		SizerL.Add(self.MaterialBox, 1, wx.GROW)
		SizerL.AddStretchSpacer(1)
		SizerL.Add(CalcBox, 0, wx.GROW)




		self.PlotAxes = plot.PlotCanvas(self)

		SizerR.Add(self.PlotAxes, 1, wx.GROW)
		#SizerR.Add(self.Rplot, 1, wx.GROW)
		# enable the zoom feature (drag a box around area of interest)
		self.PlotAxes.enableZoom = True
		#self.Rplot.SetEnableZoom(True)


		Sizer1.Add(SizerL, 1, wx.GROW)
		Sizer1.Add(SizerR, 3, wx.GROW)
		self.SetAutoLayout(True)
		self.SetSizer(Sizer1)			   # add outer sizer to frame
		self.Fit()

		self.Show(True)
		self.plot_data()
		#self.Test()

	def Test(self):
		"""Convenience function for repetitive testing"""
		self.filename = "NC-Xy_norm_bgsub.txt"
		self.dirname = "data"
		self.FileText.SetLabel("File: "+self.filename)
		#self.raw_file = self.LoadData(os.path.join(self.dirname, self.filename))
		self.AddBackgroundCheckBox.SetValue(True)
		self.combine_data()
		self.PP_AlgorithmRadio.SetValue(True)
		self.plot_data()




	def OnAbout(self, e):
		d = wx.MessageDialog(self, " A utility for calculating the real part of soft X-ray spectra written by Dr. Benjamin Watts at the Paul Scherrer Institute. If you make use of this utility, please consider citing the associated article:\n\nBenjamin Watts, ""Calculation of the Kramers-Kronig transform of X-ray spectra by a piecewise Laurent polynomial method,"" Opt. Express 22, 23628-23639 (2014)", "About KKcalc", wx.OK)
		# Create a message dialog box
		d.ShowModal() # Shows it
		d.Destroy() # finally destroy it when finished.

	def OnExit(self, e):
		self.Close(True)  # Close the frame.

	def OnOpen(self, e):
		"""Load data from a file."""
		success = False
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.FD_OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			success = True
			self.dirname, self.filename = os.path.split(dlg.GetPath())
		dlg.Destroy()
		if success:
			self.FileText.SetLabel("File: "+self.filename)
			self.raw_file = data.load_data(os.path.join(self.dirname, self.filename))
			self.combine_data()
			self.plot_data()

	def OnHelp(self, e):
		logger.info("Opening web browser for help files.")
		import webbrowser
		webbrowser.open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"README.rst"))
		#webbrowser.open("https://doi.org/10.1364/OE.22.023628")#README.rst")

	def OnArticle(self, e):
		logger.info("Opening web browser for article.")
		import webbrowser
		webbrowser.open("https://doi.org/10.1364/OE.22.023628")

	def OnSave(self, e):
		"""Write data to file."""
		convert_to = None
		if e.Id == 201:
			convert_to = "photoabsorption"
		elif e.Id == 202:
			convert_to = "refractive_index"
		logger.info("Save")
		fd = wx.FileDialog(self, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
		if fd.ShowModal()==wx.ID_OK:
			metadata = {"Density": float(self.DensityText.GetValue()), "Molecular Formula":self.StoichiometryText.GetValue(),"Formula Mass":data.calculate_FormulaMass(self.Stoichiometry)}
			data.export_data(fd.GetPath(), numpy.transpose(numpy.vstack((self.Full_E,self.KK_Real_Spectrum,data.coeffs_to_ASF(self.Full_E,self.Imaginary_Spectrum)))), header_info=metadata, convert_to=convert_to)

	def combine_data(self):
		"""Combine users near-edge data with extended spectrum data."""
		self.Full_E = None
		self.Imaginary_Spectrum = None
		if self.raw_file is not None:
			logger.info("Convert to scattering factors")
			self.NearEdgeData = data.convert_data(self.raw_file,self.DataTypeCombo.GetValue(),'ASF')
#			if self.InvertDataCheckBox.GetValue():
#				self.NearEdgeData[:,1] =  numpy.abs(self.NearEdgeData[:,1] - 2*numpy.mean(self.NearEdgeData[:,1]))
		logger.info("Combine Data")
		# Get splice points
		splice_eV = numpy.array([10.0, 30000.0])  # Henke limits
		if self.SpliceText1.GetValue() == "Start":
			if self.raw_file is not None:
				splice_eV[0] = self.NearEdgeData[0, 0]
		else:
			splice_eV[0] = float(self.SpliceText1.GetValue())
		if self.SpliceText2.GetValue() == "End":
			if self.raw_file is not None:
				splice_eV[1] = self.NearEdgeData[-1, 0]
		else:
			splice_eV[1] = float(self.SpliceText2.GetValue())
		if self.raw_file is not None and self.ASF_Data is None:
			self.Full_E, self.Imaginary_Spectrum, self.NearEdgeData, self.splice_ind  = data.merge_spectra(self.NearEdgeData, self.ASF_E, self.ASF_Data, merge_points=splice_eV, add_background=self.AddBackgroundCheckBox.GetValue(), plotting_extras=True)

		elif self.raw_file is None and self.ASF_Data is not None:
			self.Full_E = self.ASF_E
			self.Imaginary_Spectrum = self.ASF_Data

		elif self.raw_file is not None and self.ASF_Data is not None:
			
			self.Full_E, self.Imaginary_Spectrum, self.NearEdgeData, self.splice_ind  = data.merge_spectra(self.NearEdgeData, self.ASF_E, self.ASF_Data, merge_points=splice_eV, add_background=self.AddBackgroundCheckBox.GetValue(), fix_distortions=self.FixDistortionsCheckBox.GetValue(), plotting_extras=True)

			### get start and end Y values from nexafs and asf data
			##splice_nexafs_Im = numpy.interp(splice_eV, raw_Im[:, 0], raw_Im[:, 1])
			###splice_asf_Im = numpy.interp(splice_eV, self.total_asf[:, 0], self.total_asf[:, 2])
			##splice_asf_Im = (data.coeffs_to_ASF(splice_eV[0],self.total_Im_coeffs[numpy.where(self.total_E<splice_eV[0])[0][-1]]),data.coeffs_to_ASF(splice_eV[1],self.total_Im_coeffs[numpy.where(self.total_E<splice_eV[1])[0][-1]]))
			##cut_boolean = (splice_eV[0]<raw_Im[:, 0]) == (raw_Im[:, 0]<splice_eV[1])
			### Merge Y values
			##if not self.AddBackgroundCheckBox.GetValue():
				##logger.info("Merge data sets")
				##scale = (splice_asf_Im[1]-splice_asf_Im[0])/(splice_nexafs_Im[1]-splice_nexafs_Im[0])
				##scaled_nexafs_Im = ((raw_Im[:, 1]-splice_nexafs_Im[0])*scale)+splice_asf_Im[0]
				##self.asf_bg = None  # We won't be using this variable this time
			##else:
				##logger.info("Add data sets (this will currently only work at energies below 30 keV)")
				### Set up background function
				### We trust this point to be just before the absorption edge
				##trusted_ind = max(0, numpy.where(self.total_asf[:, 0]>splice_eV[0])[0][0]-1)
				##Log_total_asf = numpy.log(self.total_asf[:, 2])
				### Lets trust the 5 points before our trusted point and make an initial guess at the background function
				##p = numpy.polyfit(self.total_asf[(trusted_ind-5):trusted_ind, 0], Log_total_asf[(trusted_ind-5):trusted_ind], 1)
				### Now lets look for the points up util the absorption edge
				##p_vals = numpy.exp(numpy.polyval(p, self.total_asf[(trusted_ind-5):-1, 0]))
				##p_err = max(p_vals[0:5]-self.total_asf[(trusted_ind-5):trusted_ind, 2])
				##edge_ind = numpy.where(self.total_asf[trusted_ind:-1, 2]-p_vals[4:-1]>p_err*10)
				##if len(edge_ind[0])!=0:
					##edge_ind = edge_ind[0][0]
				##else:
					##edge_ind = trusted_ind
				### Redo background using the 5 points before the background point
				##p = numpy.polyfit(self.total_asf[(trusted_ind+edge_ind-5):trusted_ind+edge_ind, 0], Log_total_asf[(trusted_ind+edge_ind-5):trusted_ind+edge_ind], 1)
				##asf_bg = numpy.exp(numpy.polyval(p, raw_Im[:, 0]))
				##logger.info("Background defined as: y=exp(%(p1)ex %(p0)+e)" % {"p1":p[1], "p0":p[0]})
				### Apply background function
				##scale = (splice_asf_Im[1]-numpy.exp(numpy.polyval(p, splice_eV[1])))/splice_nexafs_Im[1]
				##scaled_nexafs_Im = raw_Im[:, 1]*scale+asf_bg
				### store background data for plotting
				##cut_boolean_wide = numpy.roll(cut_boolean, -1) + numpy.roll(cut_boolean, 1)
				##self.asf_bg = [[trusted_ind+edge_ind-5, trusted_ind+edge_ind], numpy.vstack((raw_Im[cut_boolean_wide, 0], asf_bg[cut_boolean_wide])).T]
			
			##nexafs_cut = numpy.vstack((raw_Im[cut_boolean, 0], scaled_nexafs_Im[cut_boolean])).T
			####Merge point-wise data sets together
			##asf_cut_high = self.total_asf[self.total_asf[:, 0]>splice_eV[1], :]
			##asf_cut_low = self.total_asf[self.total_asf[:, 0]<splice_eV[0], :]
			##self.merged_Im = numpy.vstack((asf_cut_low[:, [0, 2]], (splice_eV[0], splice_asf_Im[0]), nexafs_cut, (splice_eV[1], splice_asf_Im[1]), asf_cut_high[:, [0, 2]]))
			
			####Merge coeff data together
			##coeffs_cut_high = self.total_Im_coeffs[self.total_E[:-1]>splice_eV[1],:]
			##coeffs_cut_low = self.total_Im_coeffs[self.total_E[:-1]<splice_eV[0],:]
			###convert points to coeffs
			##nexafs_coeffs_cut = numpy.zeros((len(nexafs_cut)+1,5))
			##Y = numpy.append(numpy.insert(nexafs_cut[:,1],0,splice_asf_Im[0]),splice_asf_Im[1])
			##nexafs_E = numpy.append(numpy.insert(nexafs_cut[:,0],0,splice_eV[0]),splice_eV[1])
			##M = (Y[1:]-Y[:-1])/(nexafs_E[1:]-nexafs_E[:-1])
			##nexafs_coeffs_cut[:,0] = M
			##nexafs_coeffs_cut[:,1] = Y[:-1]-M*nexafs_E[:-1]
			###assemble merged coeffs and energy values
			##self.merged_Im_coeffs = numpy.vstack((coeffs_cut_low, nexafs_coeffs_cut, self.total_Im_coeffs[-coeffs_cut_high.shape[0]-2,:], coeffs_cut_high))
			##self.merged_E = numpy.concatenate((self.total_E[self.total_E<splice_eV[0]], nexafs_E, self.total_E[self.total_E>splice_eV[1]]))
			### Extras for plotting
			##self.splice_ind = (len(asf_cut_low[:, 0]), -len(asf_cut_high[:, 0]))
			##cut_boolean = (splice_eV[0]<=raw_Im[:, 0]) != (raw_Im[:, 0]<=splice_eV[1])
			##self.nexafs_CutOut = numpy.vstack((raw_Im[cut_boolean, 0], scaled_nexafs_Im[cut_boolean])).T
		### Previous calculation of f_1 is no longer matching displayed f_2 data
		##self.KK_Real_Spectrum = None

	def plot_data(self):
		"""Plot data.
		
		Parameters:
		-----------
		self.Full_E : vector of floats 
			photon energies at which the real and imaginary scattering factor data will be plotted.
		self.Imaginary_Spectrum : Array of float 
			polynomial coefficients that can be evaluated to give the values of the imaginary scattering factors.
		self.KK_Real_Spectrum : vector of float 
			the values of the real scattering factors.

		Returns
		-------
		The GUI is updated, but nothing is returned.
		"""
		logger.info("plotting data")
		# List of things to plot
		plotlist = []
		# get initial guess at X limits
		X_min = 0
		X_max = 30000
		Y_max = 1
		Y_min = 0
		if self.NearEdgeData is not None:
			X_min = self.NearEdgeData[0, 0]
			X_max = self.NearEdgeData[-1, 0]
		if self.SpliceText1.GetValue() != "Start":
			X_min = float(self.SpliceText1.GetValue())
		if self.SpliceText2.GetValue() != "End":
			X_max = float(self.SpliceText2.GetValue())
		if self.Imaginary_Spectrum is not None:
			if self.Stoichiometry is not None:
				scale = sum([Z*count for Z, count in self.Stoichiometry])
			else:
				scale = 1.
			Im_energies, Im_values = data.coeffs_to_linear(self.Full_E, self.Imaginary_Spectrum, 0.001*scale)
			plotlist.append(plot.PolyLine(list(zip(Im_energies,Im_values)), colour='black', width=1))
			
			# get Y limits
			if self.splice_ind is None:
				Y_max = max(Im_values)
				Y_min = min(Im_values)
			else:
				Y_max = max(Im_values[self.splice_ind[0]:self.splice_ind[1]])
				Y_min = min(Im_values[self.splice_ind[0]:self.splice_ind[1]])
		if self.NearEdgeData is not None:
			Y_max = max(self.NearEdgeData[:,1])
			Y_min = min(self.NearEdgeData[:,1])
			plotlist.append(plot.PolyMarker(list(zip(self.NearEdgeData[:,0], self.NearEdgeData[:,1])), colour='blue', marker='plus', size=1))
			
		if self.splice_ind is not None:
			splice_values = data.coeffs_to_ASF(self.Full_E[self.splice_ind], self.Imaginary_Spectrum[[self.splice_ind[0],min(self.splice_ind[1],self.Imaginary_Spectrum.shape[0]-1)]])
			plotlist.append(plot.PolyMarker(list(zip(self.Full_E[self.splice_ind], splice_values)), colour='red', marker='cross', size=1))


		if self.raw_file is not None and self.Imaginary_Spectrum is None:
			logger.info("plot raw data only")
			plotlist.append(plot.PolyLine(self.NearEdgeData, colour='blue', width=1))  # User data
			if self.asf_bg is not None:
				plotlist.append(plot.PolyMarker(self.total_asf[self.asf_bg[0][0]:self.asf_bg[0][1], [0, 2]], colour='red', marker='cross', size=1))
				plotlist.append(plot.PolyLine(self.asf_bg[1], colour='red', width=1))

			# Real part
			#plotlist.append(plot.PolyLine(self.total_asf[:, [0, 1]], colour='black', width=1))
		if self.KK_Real_Spectrum is not None:
			if self.splice_ind is None:
				Y_max = max(self.KK_Real_Spectrum)
				Y_min = min(self.KK_Real_Spectrum)
			else:
				Y_max = max(Y_max, max(self.KK_Real_Spectrum[self.splice_ind[0]:self.splice_ind[1]]))
				Y_min = min(Y_min, min(self.KK_Real_Spectrum[self.splice_ind[0]:self.splice_ind[1]]))
			plotlist.append(plot.PolyLine(list(zip(self.Full_E, self.KK_Real_Spectrum)), colour='green', width=1))

		# Expand plotting limits for prettiness
		window_width = X_max-X_min
		X_max = X_max+window_width*0.1
		X_min = max(X_min-window_width*0.1, 0)
		window_Im_height = Y_max-Y_min
		window_Re_height = Y_max-Y_min
		Y_max = Y_max+window_Im_height*0.1
		Y_min = Y_min-window_Im_height*0.1
		Y_max = Y_max+window_Re_height*0.1
		Y_min = Y_min-window_Re_height*0.1
		# set up text, axis and draw
		#print plotlist
		#print X_min, X_max, Y_min, Y_max
		self.PlotAxes.Draw(plot.PlotGraphics(plotlist, '', 'Energy (eV)', 'Magnitude'), xAxis=(X_min, X_max), yAxis=(0, Y_max))
		#print "Plotlist =", len(plotlist)

	def Splice_Text_check(self, evt):
		self.combine_data()
		self.plot_data()

	def MergeAdd_check(self, evt):
		self.combine_data()
		self.plot_data()

	def Stoichiometry_Text_check(self, evt):
		if len(self.StoichiometryText.GetValue()) == 0:
			self.ChemicalFormula = None
			self.Stoichiometry = None
			self.Relativistic_Correction = None
			self.ASF_E = None
			self.ASF_Data = None
		else:
			self.ChemicalFormula = self.StoichiometryText.GetValue()
			self.Stoichiometry = data.ParseChemicalFormula(self.ChemicalFormula)
			self.Relativistic_Correction = kk.calc_relativistic_correction(self.Stoichiometry)
			self.ASF_E, self.ASF_Data = data.calculate_asf(self.Stoichiometry)
		self.combine_data()
		self.plot_data()



	def calculate(self, button):
		"""Calculate Button."""
		logger.debug("Calculate button")
		if self.Imaginary_Spectrum is not None:
			logger.info("Calculate Kramers-Kronig transform (PP)")
			self.KK_Real_Spectrum = kk.KK_PP(self.Full_E, self.Full_E, self.Imaginary_Spectrum, self.Relativistic_Correction)
			logger.info("Done!")
			self.plot_data()

def start_wx():
	app = wx.App()
	f = MyFrame()
	app.SetTopWindow(f)
	app.MainLoop()

if __name__ == '__main__':
	start_wx()
