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

filename = ''  # full path and filename
chemical_formular = ''

output = kk.kk_calculate_real(filename,
                              chemical_formular,
                              load_options=None,
                              input_data_type='Beta',
                              merge_points=None,
                              add_background=False,
                              fix_distortions=False,
                              curve_tolerance=None,
                              curve_recursion=50)

plt.figure()
plt.plot(output[:, 0], output[:, 1], label='f1')
plt.plot(output[:, 0], output[:, 2], label='f2')
plt.legend()
plt.show()