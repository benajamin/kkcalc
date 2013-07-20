#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

"""This module implements a GUI using the wxPython toolkit."""

from wx.lib.fancytext import StaticFancyText
import math
import numpy
import os
import scipy
import scipy.fftpack
import scipy.interpolate
import scipy.io
import time  # only for profiling
import wx
import wx.lib.plot as plot

import data
import kk


class SaveFrame(wx.Frame):

	def __init__(self, Parent):
		print Parent.MolecularFormula, Parent.MolecularMass
		wx.Frame.__init__(self, None, wx.ID_ANY, "Export Data", size=(500, 800))
		SizerV = wx.BoxSizer(wx.VERTICAL)
		MetadataBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Metadata"), wx.VERTICAL)
		SizerV.Add(MetadataBox, 3, wx.GROW)

		SizerH0 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH0.Add(wx.StaticText(self, -1, "Acronym:"))
		SizerH0.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH0, 0, wx.GROW)

		SizerH1 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH1.Add(wx.StaticText(self, -1, "Name:"))
		SizerH1.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH1, 0, wx.GROW)

		SizerH2 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH2.Add(wx.StaticText(self, -1, "Formula:"))
		SizerH2.Add((1, 1), 1, wx.GROW)
		SizerH2.Add(wx.StaticText(self, -1, Parent.MolecularFormula, style=wx.ALIGN_RIGHT))
		MetadataBox.Add(SizerH2, 0, wx.GROW)

		SizerH3 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH3.Add(wx.StaticText(self, -1, "Formula Mass:"))
		SizerH3.Add((1, 1), 1, wx.GROW)
		SizerH3.Add(wx.StaticText(self, -1, str(Parent.MolecularMass), style=wx.ALIGN_RIGHT))
		MetadataBox.Add(SizerH3, 0, wx.GROW)

		SizerH4 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH4.Add(wx.StaticText(self, -1, "Density:"))
		self.DensityText = wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER)
		SizerH4.Add(self.DensityText, 1, wx.GROW)
		SizerH4.Add(StaticFancyText(self, -1, "(g/cm<sup>3</sup>)"))
		MetadataBox.Add(SizerH4, 0, wx.GROW)

		SizerH5 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH5.Add(wx.StaticText(self, -1, "Reference:"))
		SizerH5.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH5, 0, wx.GROW)

		SizerH6 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH6.Add(wx.StaticText(self, -1, "Source:"))
		SizerH6.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH6, 0, wx.GROW)

		SizerH7 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH7.Add(wx.StaticText(self, -1, "Acquisition:"))
		SizerH7.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH7, 0, wx.GROW)

		SizerH8 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH8.Add(wx.StaticText(self, -1, "Abs. Edge (Resolution):"))
		SizerH8.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		MetadataBox.Add(SizerH8, 0, wx.GROW)

		SizerH9 = wx.BoxSizer(wx.HORIZONTAL)
		SizerH9.Add(wx.StaticText(self, -1, "Incidence Angle:"))
		SizerH9.Add(wx.TextCtrl(self, -1, "", style=wx.ALIGN_RIGHT|wx.TE_PROCESS_ENTER), 1, wx.GROW)
		SizerH9.Add(wx.StaticText(self, -1, "(\xb0)".decode('ISO8859-1')))
		MetadataBox.Add(SizerH9, 0, wx.GROW)

		FormatBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Format"), wx.VERTICAL)
		FormatBox.Add(wx.RadioButton(self, -1, "Scattering Factor", style=wx.RB_GROUP))
		self.CRIRadio = wx.RadioButton(self, -1, "Complex Refractive Index")
		FormatBox.Add(self.CRIRadio)
		FormatBox.Add(wx.RadioButton(self, -1, "Photoabsorption Coefficient"))
		SizerV.Add(FormatBox, 1, wx.GROW)

		SaveButton = wx.Button(self, -1, "Save")
		SaveButton.Bind(wx.EVT_BUTTON, self.Save)
		SizerV.Add(SaveButton, 0, wx.CENTER)

		self.SetSizer(SizerV)			   # add outer sizer to frame
		self.Fit()

		self.Show(True)

	def Save(self, evt):
		print "do something"
		print self.CRIRadio.GetValue(), self.DensityText.GetValue()

		if self.CRIRadio.GetValue():
			num_value = None
			try:
				num_value = float(self.DensityText.GetValue())
			except ValueError:
				print "'", self.DensityText.GetValue(), "' is not a number!"
			if num_value is not None:
			  print "can safely save"
		# Get values

		# Convert values

		self.Close(True)  # Close the frame.

		# choose filename

		# write data



class MyFrame(wx.Frame):

	def __init__(self):
		wx.Frame.__init__(self, None, wx.ID_ANY, "Kramers-Kronig Calculator", size=(500, 800))

		# Initialise variables
		self.dirname = ''
		self.raw_file = None
		self.total_asf = None
		self.total_Im_coeffs = None
		self.merged_Im = None
		self.nexafs_CutOut = []
		self.KK_Re = None
		self.MolecularMass = 1
		self.asf_bg = None
		# Get data about elements
		self.Element_Database = data.load_Element_Database()
#		self.Elements = [line.strip("\r\n").split() for line in open(os.path.join(os.getcwd(), 'asf', 'elements.dat'))]

		# Setting up the menus.
		filemenu = wx.Menu()
		filemenu.Append(wx.ID_OPEN, "L&oad", " Load photoabsorption data from file")
		filemenu.AppendSeparator()
		filemenu.Append(wx.ID_SAVE, "&Save", " Export results to file")
		filemenu.AppendSeparator()
		filemenu.Append(wx.ID_EXIT, "E&xit", " Terminate the program")
		helpmenu = wx.Menu()
		helpmenu.Append(wx.ID_HELP, "&Help", " How to use this program")
		helpmenu.AppendSeparator()
		helpmenu.Append(wx.ID_ABOUT, "&About", " Information about this program")
		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")  # Adding the "filemenu" to the MenuBar
		menuBar.Append(helpmenu, "&Help")  # Adding the "helpmenu" to the MenuBar
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		wx.EVT_MENU(self, wx.ID_OPEN, self.OnOpen)
		wx.EVT_MENU(self, wx.ID_SAVE, self.OnSave)
		wx.EVT_MENU(self, wx.ID_EXIT, self.OnExit)
		wx.EVT_MENU(self, wx.ID_ABOUT, self.OnAbout)
		wx.EVT_MENU(self, wx.ID_HELP, self.OnHelp)


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
		self.AddBackgroundCheckBox = wx.CheckBox(self, -1, "Add background")
		self.AddBackgroundCheckBox.Bind(wx.EVT_CHECKBOX, self.Splice_Text_check)
		DataBox.Add(self.AddBackgroundCheckBox, 0)
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
		self.MaterialBox.Add(wx.StaticText(self, -1, "Stoichiometry: "))
		self.MaterialBox.contents = []
		self.add_element()
		self.MaterialBox.AddStretchSpacer(1)

		############################Calc box
		CalcBox = wx.StaticBoxSizer(wx.StaticBox(self, label="Calculation"), wx.VERTICAL)
		self.PP_AlgorithmRadio = wx.RadioButton(self, -1, "Piecewise-polynomial", style=wx.RB_GROUP)
		self.FFT_AlgorithmRadio = wx.RadioButton(self, -1, "FFT-based")
		CalcBox.Add(self.PP_AlgorithmRadio, 1, wx.GROW)
		CalcBox.Add(self.FFT_AlgorithmRadio, 1, wx.GROW)
		CalcButton = wx.Button(self, -1, "Calculate")
		CalcBox.Add(CalcButton, 1, wx.GROW)
		CalcButton.Bind(wx.EVT_BUTTON, self.calculate)



		SizerL.Add(DataBox, 0, wx.GROW)
		SizerL.Add(self.MaterialBox, 1, wx.GROW)
		SizerL.Add(CalcBox, 0, wx.GROW)




		self.Iplot = plot.PlotCanvas(self)
		self.Rplot = plot.PlotCanvas(self)

		SizerR.Add(self.Iplot, 1, wx.GROW)
		SizerR.Add(self.Rplot, 1, wx.GROW)
		# enable the zoom feature (drag a box around area of interest)
		self.Iplot.SetEnableZoom(True)
		self.Rplot.SetEnableZoom(True)


		Sizer1.Add(SizerL, 1, wx.GROW)
		Sizer1.Add(SizerR, 3, wx.GROW)
		self.SetAutoLayout(True)
		self.SetSizer(Sizer1)			   # add outer sizer to frame
		self.Fit()

		self.Show(True)
		self.plot_data()
		self.Test()

	def Test(self):
		"""Convenience function for repetitive testing"""
		self.filename = "NC-Xy_norm_bgsub.txt"
		self.dirname = "data"
		self.FileText.SetLabel("File: "+self.filename)
		self.raw_file = self.LoadData(os.path.join(self.dirname, self.filename))
		self.AddBackgroundCheckBox.SetValue(True)
		self.combine_data()
		self.PP_AlgorithmRadio.SetValue(True)
		self.plot_data()




	def OnAbout(self, e):
		d = wx.MessageDialog(self, " A utility for calculating the real part of soft X-ray spectra.\nWritten by Dr. Benjamin Watts at the Paul Scherrer Institut", "About KKcalc", wx.OK)
		# Create a message dialog box
		d.ShowModal() # Shows it
		d.Destroy() # finally destroy it when finished.

	def OnExit(self, e):
		self.Close(True)  # Close the frame.

	def OnOpen(self, e):
		"""Load data from a file."""
		success = False
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			success = True
			self.dirname, self.filename = os.path.split(dlg.GetPath())
		dlg.Destroy()
		if success:
			self.FileText.SetLabel("File: "+self.filename)
			self.raw_file = self.LoadData(os.path.join(self.dirname, self.filename))

			self.combine_data()
			self.plot_data()

	def ConvertData(self, raw_data):
		if len(raw_data) != 0:
			data_type = self.DataTypeCombo.GetSelection()
			if data_type == 1:  # Beta
				print "Convert from Beta (type=", data_type, ") to Scattering Factors."
				density = float(self.DensityText.GetValue())
				raw_Im = data.convert_Beta_to_ASF(raw_data, density)
			elif data_type == 2:  # Scattering factor
				print "Data is already in terms of Scattering Factors (type=", data_type, ")."
			else:  # (data_type == 0 or -1) Assume Photoabsorption data
				print "Convert NEXAFS photoabsorption data to Scattering Factors."
				raw_Im = data.convert_NEXAFS_to_ASF(raw_data)
			return raw_Im

	def OnHelp(self, e):
		print "Opening web browser for help files."
		import webbrowser
		webbrowser.open("KKcalc.html")

	def OnSave(self, e):
		"""Write data to file."""
		print "Save"
		f = SaveFrame(self)
		print "done!"
		if self.KK_Re is not None:
			fd = wx.FileDialog(self, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
			if fd.ShowModal()==wx.ID_OK:
				outfile = open(os.path.join(fd.GetDirectory(), fd.GetFilename()), "w")
				outfile.write('Scattering factors for '+self.MolecularFormula+'\n')
				outfile.write('E(eV)\tf1\tf2\n')
				for i in xrange(len(self.merged_Im[:, 0])):
	#				outfile.write("{0}\t{1}\t{2}\n".format(self.merged_Im[i, 0], self.KK_Re[i], self.merged_Im[i, 1]))  # Python 3.0 style
					outfile.write("%(E)#7g\t%(Re)#7g\t%(Im)#7g\n"%{'E':self.merged_Im[i, 0], 'Re':self.KK_Re[i], 'Im':self.merged_Im[i, 1]})  # old formatting style
				outfile.close()
			print "Scattering factors for", self.MolecularFormula, "saved to ", fd.GetFilename()
		else:
			print "Nothing to save."

	def LoadData(self, filename):
		"""Read a standard ASCII file and return a list of lists of floats."""
		data = []
		if os.path.isfile(filename):
			for line in open(filename):
				try:
					data.append([float(f) for f in line.split()])
				except ValueError:
					pass
			data = numpy.array(data)
		else:
			print "Error:", filename, "is not a valid file name."
		if len(data)==0:
			print "Error: no data found in", filename
			return []
		else:
			return data

	def combine_data(self):
		"""Combine users near-edge data with extended spectrum data."""
		#TODO do all data scaling using the coeff extended spectrum data (not with the point-wise extended spectrum data)
		if self.raw_file is not None:
			print "Convert to scattering factors"
			raw_Im = self.ConvertData(self.raw_file)
		print "Combine Data"
		# Get splice points
		splice_eV = numpy.array([10, 30000])  # Henke limits
		if self.SpliceText1.GetValue() == "Start":
			if self.raw_file is not None:
				splice_eV[0] = raw_Im[0, 0]
		else:
			splice_eV[0] = float(self.SpliceText1.GetValue())
		if self.SpliceText2.GetValue() == "End":
			if self.raw_file is not None:
				splice_eV[1] = raw_Im[-1, 0]
		else:
			splice_eV[1] = float(self.SpliceText2.GetValue())

		if self.raw_file is not None and self.total_asf is None:
			splice_nexafs_Im = numpy.interp(splice_eV, raw_Im[:, 0], raw_Im[:, 1])
			cut_boolean = (splice_eV[0]<raw_Im[:, 0]) == (raw_Im[:, 0]<splice_eV[1])
			nexafs_cut = raw_Im[cut_boolean]
			self.merged_Im = numpy.vstack(((splice_eV[0], splice_nexafs_Im[0]), nexafs_cut, (splice_eV[1], splice_nexafs_Im[1])))
			# Extras for plotting
			self.splice_ind = (0, -1)
			cut_boolean = (splice_eV[0]<=raw_Im[:, 0]) != (raw_Im[:, 0]<=splice_eV[1])
			self.nexafs_CutOut = raw_Im[cut_boolean]
			self.asf_bg = None  # We won't be using this variable this time

		elif self.raw_file is None and self.total_asf is not None:
			self.merged_Im = self.total_asf[:, [0, 2]]
			# Extras for plotting
			self.splice_ind = (0, -1)
			self.nexafs_CutOut = None
			self.asf_bg = None  # We won't be using this variable this time

		elif self.raw_file is not None and self.total_asf is not None:
			# get start and end Y values from nexafs and asf data
			splice_nexafs_Im = numpy.interp(splice_eV, raw_Im[:, 0], raw_Im[:, 1])
			splice_asf_Im = numpy.interp(splice_eV, self.total_asf[:, 0], self.total_asf[:, 2])
			cut_boolean = (splice_eV[0]<raw_Im[:, 0]) == (raw_Im[:, 0]<splice_eV[1])
			# Merge Y values
			if not self.AddBackgroundCheckBox.GetValue():
				print "Merge data sets"
				scale = (splice_asf_Im[1]-splice_asf_Im[0])/(splice_nexafs_Im[1]-splice_nexafs_Im[0])
				scaled_nexafs_Im = ((raw_Im[:, 1]-splice_nexafs_Im[0])*scale)+splice_asf_Im[0]
				self.asf_bg = None  # We won't be using this variable this time
			else:
				print "Add data sets"
				# Set up background function
				# We trust this point to be just before the absorption edge
				trusted_ind = max(0, numpy.where(self.total_asf[:, 0]>splice_eV[0])[0][0]-1)
				Log_total_asf = numpy.log(self.total_asf[:, 2])
				# Lets trust the 5 points before our trusted point and make an initial guess at the background function
				p = numpy.polyfit(self.total_asf[(trusted_ind-5):trusted_ind, 0], Log_total_asf[(trusted_ind-5):trusted_ind], 1)
				# Now lets look for the points up util the absorption edge
				p_vals = numpy.exp(numpy.polyval(p, self.total_asf[(trusted_ind-5):-1, 0]))
				p_err = max(p_vals[0:5]-self.total_asf[(trusted_ind-5):trusted_ind, 2])
				edge_ind = numpy.where(self.total_asf[trusted_ind:-1, 2]-p_vals[4:-1]>p_err*10)
				if len(edge_ind[0])!=0:
					edge_ind = edge_ind[0][0]
				else:
					edge_ind = trusted_ind
				# Redo background using the 5 points before the background point
				p = numpy.polyfit(self.total_asf[(trusted_ind+edge_ind-5):trusted_ind+edge_ind, 0], Log_total_asf[(trusted_ind+edge_ind-5):trusted_ind+edge_ind], 1)
				asf_bg = numpy.exp(numpy.polyval(p, raw_Im[:, 0]))
				print "Background defined as: y=exp(%(p1)ex %(p0)+e)" % {"p1":p[1], "p0":p[0]}
				# Apply background function
				scale = (splice_asf_Im[1]-numpy.exp(numpy.polyval(p, splice_eV[1])))/splice_nexafs_Im[1]
				scaled_nexafs_Im = raw_Im[:, 1]*scale+asf_bg
				# store background data for plotting
				cut_boolean_wide = numpy.roll(cut_boolean, -1) + numpy.roll(cut_boolean, 1)
				self.asf_bg = [[trusted_ind+edge_ind-5, trusted_ind+edge_ind], numpy.vstack((raw_Im[cut_boolean_wide, 0], asf_bg[cut_boolean_wide])).T]
			
			nexafs_cut = numpy.vstack((raw_Im[cut_boolean, 0], scaled_nexafs_Im[cut_boolean])).T
			##Merge point-wise data sets together
			asf_cut_high = self.total_asf[self.total_asf[:, 0]>splice_eV[1], :]
			asf_cut_low = self.total_asf[self.total_asf[:, 0]<splice_eV[0], :]
			self.merged_Im = numpy.vstack((asf_cut_low[:, [0, 2]], (splice_eV[0], splice_asf_Im[0]), nexafs_cut, (splice_eV[1], splice_asf_Im[1]), asf_cut_high[:, [0, 2]]))
			
			##Merge coeff data together
			coeffs_cut_high = self.total_Im_coeffs[self.total_E[:-1]>splice_eV[1],:]
			coeffs_cut_low = self.total_Im_coeffs[self.total_E[:-1]<splice_eV[0],:]
			#convert points to coeffs
			nexafs_coeffs_cut = numpy.zeros((len(nexafs_cut)+1,5))
			Y = numpy.append(numpy.insert(nexafs_cut[:,1],0,splice_asf_Im[0]),splice_asf_Im[1])
			nexafs_E = numpy.append(numpy.insert(nexafs_cut[:,0],0,splice_eV[0]),splice_eV[1])
			M = (Y[1:]-Y[:-1])/(nexafs_E[1:]-nexafs_E[:-1])
			nexafs_coeffs_cut[:,0] = M
			nexafs_coeffs_cut[:,1] = Y[:-1]-M*nexafs_E[:-1]
			#assemble merged coeffs and energy values
			self.merged_Im_coeffs = numpy.vstack((coeffs_cut_low, nexafs_coeffs_cut, self.total_Im_coeffs[-coeffs_cut_high.shape[0]-2,:], coeffs_cut_high))
			self.merged_E = numpy.concatenate((self.total_E[self.total_E<splice_eV[0]], nexafs_E, self.total_E[self.total_E>splice_eV[1]]))
			# Extras for plotting
			self.splice_ind = (len(asf_cut_low[:, 0]), -len(asf_cut_high[:, 0]))
			cut_boolean = (splice_eV[0]<=raw_Im[:, 0]) != (raw_Im[:, 0]<=splice_eV[1])
			self.nexafs_CutOut = numpy.vstack((raw_Im[cut_boolean, 0], scaled_nexafs_Im[cut_boolean])).T
		# Previous calculation of f_1 is no longer matching displayed f_2 data
		self.KK_Re = None

	def plot_data(self):
		"""Plot data."""
		print "plotting data"
		# List of things to plot
		plotlist_Im = []
		plotlist_Re = []
		# get initial guess at X limits
		X_min = 0
		X_max = 30000
		Y_Im_max = 1
		Y_Im_min = 0
		Y_Re_max = 1
		Y_Re_min = 0
		if self.raw_file is not None:
			X_min = self.raw_file[0, 0]
			X_max = self.raw_file[-1, 0]
		if self.SpliceText1.GetValue() != "Start":
			X_min = float(self.SpliceText1.GetValue())
		if self.SpliceText2.GetValue() != "End":
			X_max = float(self.SpliceText2.GetValue())
		if self.raw_file is not None:
			print "plot raw data only"
			# get Y limits
			Y_Im_max = max(self.merged_Im[self.splice_ind[0]:self.splice_ind[1], 1])
			Y_Im_min = min(self.merged_Im[self.splice_ind[0]:self.splice_ind[1], 1])
			plotlist_Im.append(plot.PolyLine(self.merged_Im[self.splice_ind[0]:self.splice_ind[1], :], colour='blue', width=1))  # User data
			if len(self.nexafs_CutOut)!=0:
				plotlist_Im.append(plot.PolyMarker(self.nexafs_CutOut, colour='blue', marker='cross', size=1))
			if self.asf_bg is not None:
				plotlist_Im.append(plot.PolyMarker(self.total_asf[self.asf_bg[0][0]:self.asf_bg[0][1], [0, 2]], colour='red', marker='cross', size=1))
				plotlist_Im.append(plot.PolyLine(self.asf_bg[1], colour='red', width=1))

		if self.total_asf is not None:
			if self.raw_file is not None:
				print "plot everything"
				# get Y limits
				Y_Im_max = max(self.merged_Im[self.splice_ind[0]:self.splice_ind[1], 1])
				Y_Im_min = min(self.merged_Im[self.splice_ind[0]:self.splice_ind[1], 1])
				plotlist_Im.append(plot.PolyLine(self.merged_Im[0:(self.splice_ind[0]+1), :], colour='black', width=1))  # Low energy Henke data
				plotlist_Im.append(plot.PolyLine(self.merged_Im[(self.splice_ind[1]-1):-1, :], colour='black', width=1))  # High energy Henke data
			else:
				print "plot asf data only"
				Y_Im_max = max(self.total_asf[(X_min<=self.total_asf[:, 0])==(self.total_asf[:, 0]<=X_max), 2])
				Y_Im_min = min(self.total_asf[(X_min<=self.total_asf[:, 0])==(self.total_asf[:, 0]<=X_max), 2])
				plotlist_Im.append(plot.PolyLine(self.total_asf[:, [0, 2]], colour='black', width=1))


			# Real axes
			total_asf_CUT = (X_min<=self.total_asf[:, 0])==(self.total_asf[:, 0]<=X_max)
			Y_Re_max = max(self.total_asf[total_asf_CUT, 1])
			rel_corr = kk.calc_relativistic_correction(self.Z, self.stoichiometry)
			Y_Re_min = min(self.total_asf[total_asf_CUT & (self.total_asf[:, 1]>-100*rel_corr), 1])
			plotlist_Re.append(plot.PolyLine(self.total_asf[:, [0, 1]], colour='black', width=1))
		if self.KK_Re is not None:
			Y_Re_max = max(Y_Re_max, max(self.KK_Re[self.splice_ind[0]:self.splice_ind[1]]))
			Y_Re_min = min(Y_Re_min, min(self.KK_Re[self.splice_ind[0]:self.splice_ind[1]]))
			plotlist_Re.append(plot.PolyLine(numpy.vstack((self.merged_Im[:, 0], self.KK_Re)).T, colour='green', width=1))

		# Expand plotting limits for prettiness
		window_width = X_max-X_min
		X_max = min(X_max+window_width*0.1, 30000)
		X_min = max(X_min-window_width*0.1, 0)
		window_Im_height = Y_Im_max-Y_Im_min
		window_Re_height = Y_Re_max-Y_Re_min
		Y_Im_max = Y_Im_max+window_Im_height*0.1
		Y_Im_min = Y_Im_min-window_Im_height*0.1
		Y_Re_max = Y_Re_max+window_Re_height*0.1
		Y_Re_min = Y_Re_min-window_Re_height*0.1
		# set up text, axis and draw
		self.Iplot.Draw(plot.PlotGraphics(plotlist_Im, '', '', 'Imaginary'), xAxis=(X_min, X_max), yAxis=(0, Y_Im_max))
		self.Rplot.Draw(plot.PlotGraphics(plotlist_Re, '', 'Energy (eV)', 'Real'), xAxis=(X_min, X_max), yAxis=(Y_Re_min, Y_Re_max))

	def select_element(self, evt):
		"""Select an element."""
		cb = evt.GetEventObject()
		datalist = []
		n = 0
		N = None
		for box in self.MaterialBox.contents:
			datalist.append(box[1].GetSelection())
			if box[1]==cb:
				N = n
			n = n+1
		#print 'datalist =', datalist
		try:  # check corresponding textctrl and put in a 1 if empty
			num_value = float(self.MaterialBox.contents[N][2].GetValue())
		except ValueError:
			self.MaterialBox.contents[N][2].SetValue("1")

		if evt.GetSelection() is 0:  # might need to remove a box
			if 0 in datalist and datalist.index(0) != len(datalist)-1:
				# remove extra boxes
				dead_box = datalist.index(0)
				print dead_box
				# Remove and destroy combobox
#				self.MaterialBox.contents[dead_box][0].Clear(True)
				self.MaterialBox.contents[dead_box][0].Detach(0)
				self.MaterialBox.contents[dead_box][1].Destroy()
				# Remove and destroy textbox
				self.MaterialBox.contents[dead_box][0].Detach(0)
				self.MaterialBox.contents[dead_box][2].Destroy()
				# Remove (no need to destroy) ElementSizer
				self.MaterialBox.Remove(dead_box+2)
				# delete reference to items in saved list
				self.MaterialBox.contents.remove(self.MaterialBox.contents[dead_box])
				# adjust GUI layout
				self.Layout()
		else:
			if 0 not in datalist:
				# add spare box
				#print 'Add box'
				self.add_element()
				self.Layout()
				if self.MaterialBox.GetSize()[1] < self.MaterialBox.CalcMin()[1]:
					self.Fit()
		self.calc_asfdata()

	def add_element(self):
		"""Add element GUI items."""
		# make GUI objects
		element_ComboBox = wx.ComboBox(self, -1, value='', style=wx.CB_READONLY)
		self.populate_elements(element_ComboBox)
		element_Text = wx.TextCtrl(self, -1, '', style=wx.TE_PROCESS_ENTER)
		# put into a sizer
		ElementSizer = wx.BoxSizer(wx.HORIZONTAL)
		ElementSizer.contents = (element_ComboBox, element_Text)
		ElementSizer.Add(element_ComboBox, 1)#, wx.Grow)
		ElementSizer.Add(element_Text, 1)#, wx.Grow)
		#print "Insert in materialsbox at:", max(1, len(self.MaterialBox.GetChildren())-2)
		self.MaterialBox.Insert(max(2, len(self.MaterialBox.GetChildren())-1), ElementSizer, 0, wx.GROW)
		# various bookkeeping
		element_ComboBox.Bind(wx.EVT_COMBOBOX, self.select_element)
#		element_Text.Bind(wx.EVT_TEXT, self.element_Text_check)
		element_Text.Bind(wx.EVT_KILL_FOCUS, self.element_Text_check)
		element_Text.Bind(wx.EVT_TEXT_ENTER, self.element_Text_check)
		self.MaterialBox.contents.append([ElementSizer, element_ComboBox, element_Text])
		#print "there are now", len(self.MaterialBox.GetChildren())-3, "element boxes"

	def Splice_Text_check(self, evt):
		self.combine_data()
		self.plot_data()

	def MergeAdd_check(self, evt):
		self.combine_data()
		self.plot_data()

	def element_Text_check(self, evt):
		tb = evt.GetEventObject()
		num_value = None
		try:
			num_value = float(tb.GetValue())
		except ValueError:
			print "'", tb.GetValue(), "' is not a number!"
		if num_value != None and numpy.isfinite(num_value):
			if num_value < 0:
				print tb.GetValue(), " is not a positive number!"
			else:
				print "Execute with", num_value
				self.calc_asfdata()

	def populate_elements(self, cb):
		"""Append elements."""
		cb.Append('')
		for element in range(1,93):
			cb.Append(self.Element_Database[str(element)]['symbol'])

	def calc_asfdata(self):
		"""Calculate atomic scattering factors."""
		print "Calculate total asf data"
		# start from clean slate
		self.total_asf_Re = None
		self.total_asf_Im = None
		# get element list
		self.stoichiometry = []
		self.Z = []
		self.MolecularFormula = ""
		self.MolecularMass = 0
		for i in range(len(self.MaterialBox.contents)-1):
			try:  # sanitise inputs
				self.stoichiometry.append(float(self.MaterialBox.contents[i][2].GetValue()))
				self.Z.append(int(self.MaterialBox.contents[i][1].GetSelection()))  # Need Z for relativistic correction
				# Construct molecular formula for savefile header
				num = self.stoichiometry[-1]
				if num==round(num):num=int(num)
				if num!=0:
					self.MolecularFormula = self.MolecularFormula+self.MaterialBox.contents[i][1].GetValue()+str(num)
					self.MolecularMass = self.MolecularMass+num*float(self.Element_Database[str(self.MaterialBox.contents[i][1].GetSelection())]['mass'])
					#print self.MolecularMass
			except ValueError:
				pass
		if len(self.stoichiometry)!=0:
			# get unique energy points
			temp_E = numpy.array([])
			for element in self.Z:
				temp_E = numpy.concatenate((temp_E, self.Element_Database[str(element)]['E']))
			temp_E = numpy.unique(temp_E)
			# add weighted asf data sets for KK calculation
			self.total_Im_coeffs = numpy.zeros((len(temp_E)-1, 5))
			self.total_E = temp_E
			counters = numpy.zeros((len(self.Z)))
			for i,E in enumerate(temp_E[1:]):
				#print "\nE =", E, "\t counters =", counters, "\t compare = ",
				sum_Im_coeffs = 0
				for j in range(len(counters)):
					sum_Im_coeffs += self.stoichiometry[j]*self.Element_Database[str(self.Z[j])]['Im'][counters[j],:]
					counters[j] += self.Element_Database[str(self.Z[j])]['E'][counters[j]+1] == E
				self.total_Im_coeffs[i,:] = sum_Im_coeffs
			# add weighted data sets for plotting
			Henke_end_index = numpy.where(temp_E==30000)[0][0]+1
			self.total_asf = numpy.zeros((Henke_end_index, 3))
			self.total_asf[:, 0] = temp_E[:Henke_end_index]
			counters = numpy.zeros((len(self.Z)))
			for i in range(Henke_end_index):
				for j in range(len(self.Z)):
					if self.Element_Database[str(self.Z[j])]['E'][counters[j]] == temp_E[i]:#keeping pace
						self.total_asf[i, 1] += self.stoichiometry[j]*self.Element_Database[str(self.Z[j])]['Re'][counters[j]]
						counters[j] += 1
					else: #need to interpolate extra point
						# Use method of Y = y0 + (y1-y0)*[(X-x0)/(x1-x0)]
						X_fraction = (temp_E[i] - self.Element_Database[str(self.Z[j])]['E'][counters[j]-1])/(self.Element_Database[str(self.Z[j])]['E'][counters[j]] - self.Element_Database[str(self.Z[j])]['E'][counters[j]-1])
						self.total_asf[i, 1] += self.stoichiometry[j]*(self.Element_Database[str(self.Z[j])]['Re'][counters[j]-1]+X_fraction*(self.Element_Database[str(self.Z[j])]['Re'][counters[j]]-self.Element_Database[str(self.Z[j])]['Re'][counters[j]-1]))
			self.total_asf[:, 2] = self.Coeffs_to_ASF(temp_E[:Henke_end_index],self.total_Im_coeffs[:Henke_end_index,:])
		# resplice with raw nexafs data, if any is loaded
		self.combine_data()
		# plot asf data
		self.plot_data()

	def calculate(self, button):
		"""Calculate Button."""
		print "Calculate button"
		if self.merged_Im is not None:
			tic = time.time()
			if self.FFT_AlgorithmRadio.GetValue():
				self.KK_FFT()
			else:
				# self.KK_PP()
				self.KK_PP_BL()
			print "Completed in ", round(time.time()-tic, 3), "seconds."
			self.plot_data()

	def KK_FFT(self):
		"""Calculate Kramers-Kronig transform with FFT algorithm."""
		print "Calculate Kramers-Kronig transform (FFT)"
		rel_corr = kk.calc_relativistic_correction(self.Z, self.stoichiometry)
		self.KK_Re = kk.KK_FFT(self.merged_Im, rel_corr)
		print "Done!"

	def KK_PP(self):
		"""Calculate Kramers-Kronig transform with
		"Piecewise Polynomial" algorithm.

		"""
		print "Calculate Kramers-Kronig transform (PP)"
		rel_corr = kk.calc_relativistic_correction(self.Z, self.stoichiometry)
		self.KK_Re = kk.KK_PP(self.merged_Im, rel_corr)
		print "Done!"

	def Coeffs_to_ASF(self, E, coeffs):
		"""Calculate Henke scattering factors from polynomial coefficients.

		{E in eV and PECS in cm^2/atom}.

		"""
		return coeffs[:,0]*E + coeffs[:,1] + coeffs[:,2]/E + coeffs[:,3]/(E**2) + coeffs[:,4]/(E**3)

	def KK_PP_BL(self):
		"""Calculate Kramers-Kronig transform with "Piecewise Polynomial"
		algorithm plus the Biggs and Lighthill extended data.

		"""
		print "Calculate Kramers-Kronig transform (PP) plus BL data"
		rel_corr = kk.calc_relativistic_correction(self.Z, self.stoichiometry)
		self.KK_Re = kk.KK_PP_BL(self.merged_Im, rel_corr,
								 self.BL_coefficients, self.BL_range)
		print "Done!"


if __name__ == '__main__':
	app = wx.App()
	f = MyFrame()
	app.SetTopWindow(f)
	app.MainLoop()
