#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the Kramers-Kronig Calculator software package.
#
# Copyright (c) 2013 Benjamin Watts, Daniel J. Lauk
#
# The software is licensed under the terms of the zlib/libpng license.
# For details see LICENSE.txt

import kkcalc as kk
import numpy as np
import matplotlib.pyplot as plt

filename = 'Fe.cf'  # full path and filename
chemical_formula = 'Fe'

x_min = 690
x_max = 740

output = kk.kk_calculate_real(filename,
                              chemical_formula,
                              load_options=None,
                              input_data_type='Beta',
                              merge_points=[x_min, x_max],
                              add_background=False,
                              fix_distortions=False,
                              curve_tolerance=None,
                              curve_recursion=50)

Stoichiometry = kk.data.ParseChemicalFormula(chemical_formula)
ASF_E, ASF_Data = kk.data.calculate_asf(Stoichiometry)
ASF_Data2 = kk.data.coeffs_to_ASF(ASF_E, np.vstack((ASF_Data, ASF_Data[-1])))

plt.figure()
plt.plot(output[:, 0], output[:, 1], label='f1 kkcalc')
plt.plot(output[:, 0], output[:, 2], label='f2 kkcalc')
plt.plot(ASF_E, ASF_Data2, label='Henke f2')
plt.legend()
plt.xlim(x_min, x_max)
plt.title('{:d} eV - {:d} eV'.format(x_min, x_max))
plt.xlabel('Energy [eV]')
plt.ylabel('f1, f2')
plt.show()
