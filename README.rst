======
KKcalc
======


Introduction
============

KKcalc is a program to calculate the real part of the refractive index of a material in the X-ray region of the spectrum from measured absorption data via the Kramers-Kronig transform. KKcalc boasts the following features:

- Easily extend measured spectra with published scattering factors.
- Efficient peicewise polynomial (direct integration) and FFT-based algorithms for Kramers-Kronig transform.
- Automatic selection of appropriate FFT sampling rate.
- Automatic calculation of the relativistic correction.
- Python source code is easy to read, so you can assess its correctness yourself.
- Freely available on a wide range of platforms.
- Graphical interface for easy use.

For support, please contact benjaminDOTwattsATgmailDOTcom.

The major emphasis of this program is to provide **correct** results. I would be glad to be informed of any problems.


Installation
============

Simply download the software (e.g. as a zip file) and run `KKcalc.py`. The zip file contains::

    KKcalc.py
    LICENSE.txt
    README.rst
    asf/
        elements.dat
        h.nff
        he.nff
        ...
        u.nff

You will also need to install the following dependencies:

- Python_
- Numpy_ (Numerical Python)
- Scipy_ (Scientific Python)
- wxPython_

.. _Python: http://www.python.org/
.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://scipy.org/
.. _wxPython: http://wxpython.org/


Basic Usage
===========

1. Load NEXAFS data from file by selecting: File>Load. Data is converted from photoabsorption to atomic scattering factors and plotted upon successful loading.
2. Enter molecular formula in the "Material" section.
3. Click on the "Calculate" button. The new real scattering factors are plotted upon completion of the Kramers-Kronig transform calculation (should only take a few seconds if the computer has sufficient power).
4. Save data (as molecular scattering factors) by selecting: File>Save.


Details
=======


"Near-Edge Data" Section
------------------------

Here we define the imaginary spectrum to be inputted into the Kramers-Kronig transform. The "File:" label simply displays the currently loaded data file. To load a new set of near-edge absorption data, select: *File>Load*. Data files must be in two columns, the first giving the photon energy in electron volts and the second giving the absorption spectrum data.

The two text-boxes in the next row, initially containing "Start" and "End" control the energy points between which the loaded data will be used. Valid entries in these controls include energy values and the strings "Start" or "End", which are understood by the program to correspond to the lowest and highest energy points, respectively, in the loaded data set.

If the "Extend with Henke data" checkbox is checked, then the atomic scattering factors published by Henke et al. [HEN1993]_ will be calculated from the information given in the "Material" section and used to extend the loaded data to the range 10 eV to 30,000 eV. Splicing together the data loaded from file and the Henke scattering factor data can be done in two ways:

1. The "Merge" selection will simply rescale the loaded data to fit the Henke scattering factor data.
2. The "Add" selection will calculate an appropriate background function and add this to the loaded data. Use this if the loaded data has previously had a background subtracted.


"Material" Section
------------------

In this section, we define the elemental composition of the material in question. Elements can be chosen from the combo-boxes on the left and the corresponding proportion can be entered in the text-box to its right. New sets of controls will appear (and disappear) as needed. The data entered here will be used to calculate the appropriate Henke scattering factors to extend the imaginary spectrum as well as the relativistic correction to the Kramers-Kronig transform.


"Calculation" Section
---------------------

Currently, two Kramers-Kronig algorithms are implemented: one FFT-based and the other using a piecewise polynomial, direct integration scheme. These algorithms should produce identical output. The two algorithms also utilise the computers resources differently, with the FFT-based algorithm being very RAM intensive while the piecewise polynomial algorithm is more processor limited, though much more efficient. No "step size" or "quality" parameters are required from the user because the most appropriate choice is made automatically in the case of the FFT-based algorithm and is unnecessary in the case of the piecewise polynomial algorithm.

The FFT-based algorithm is described by Bruzzoni et al. [BRU2002]_ and requires evenly spaced samples across the entire imaginary spectrum. To ensure that no information is lost during sampling, a step size equal to half of the smallest step in the loaded data (i.e. double the Nyquist frequency) is utilised. The outputted real spectrum is then re-interpolated to the same energy values as the input data.

The "Piecewise Polynomial" algorithm performs direct integration of the area between each data-point by assuming the imaginary spectrum is a straight line in each interval (i.e. linear interpolation). With this assumption, the symbolic form of the Kramers-Kronig transform integral is precisely known and can be fully written symbolically (albeit piecewise). This form is then trivial to symbolically integrate in a piecewise fashion everywhere except at the singularity, which is avoided by integrating across two intervals at once (terms referencing the singularity cancel out). The only assumption of this method is that the interpolation is a good representation of the imaginary spectrum, all remaining steps are exact to machine precision. This algorithm is very efficient because it doesn't require equally spaced steps, which would correspond to a very large number of samples over the full energy range of the spectrum.

This program uses linear interpolation (or assumes a piecewise-linear form in the case of the piecewise polynomial algorithm) only. The appropriateness of this varies according to the curvature of the imaginary spectrum and the spacing of the recorded data points. However, since it is common practice to choose, and often to vary also, the energy spacing of measured data points to reflect the curvature of the expected imaginary spectrum, then the piecewise-linear assumption will tend to be as valid as the user cares for it to be. It would be trivial to implement higher order interpolation (and assume a higher order polynomial for the piecewise polynomial algorithm), but a polynomial order of one was chosen here because it was deemed good enough, while higher order interpolation showed a great deal of oscillation near the sharp absorption edges of the Henke data used to extend the energy range of the measured imaginary spectrum. It is possible to apply an adaptive interpolation algorithm for a smoothly varying interpolated function, while eliminating oscillation at the Henke data absorption edges and this may become a future enhancement.

The calculation of the relativistic correction deserves some mention too, since I have seen a number of programs not calculating it correctly. Information on the types and number of atoms present are taken from the "Material" box and the equation :math:`Z - (\frac{Z}{82.5})^{2.37}` (as described by Henke et al. [HEN1993]_) is applied to each atom separately and the individual corrections then summed. The same code calculates the relativistic correction regardless of which Kramers-Kronig algorithm is being used.

The absorption contributions above 30 keV are not included in the published Henke data and thus are not considered in the calculation of the real spectrum. The effect of this omission will be negligible for light elements and most significant for heavy elements. The real spectrum published by Henke et al. includes absorption contributions up to 500 keV and so can be considered more accurate than that calculated by this program. One can observe the difference by comparing the black line of the real Henke data with the green line plotted for the real spectrum calculated by this program from the imaginary Henke data.


Plots
-----

The plot area contains two plot, the imaginary spectrum at the top and the real spectrum below. Loaded imaginary, near-edge data is plotted in blue. Blue crosses mark loaded data points that have been excluded from input into the Kramers-Kronig transform calculation. The background function (if selected) is plotted in red, with red crosses marking the scattering factor data-points used in its construction. Henke scattering factor data (both real and imaginary) is plotted in black. The newly calculated real spectrum is plotted in green


References
----------

.. [BRU2002] P. Bruzzoni, R.M. Carranza, J.R. Collet Lacoste, and E.A. Crespo
  "Kramers-Kronig transforms calculation with a fast convolution algorithm"
  *Electrochimica Acta* **48** (2002) 341-347.

.. [HEN1993] B.L. Henke, E.M. Gullikson, and J.C. Davis
  "X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92"
  *Atomic Data and Nuclear Data Tables* **54**\ (2) (1993) 181-342.

.. [MOH2008] N. Mohankumar, and A. Natarajan
  "On the numerical solution of Cauchy singular integral equations in neutron transport"
  *Annals of Nuclear Energy* **35**\ (10) (2008) 1800-1804.


Known Issues
============

- Critical error when selecting blank entry of a combo-box in the "Material" section. This is a known bug in wxPython and I hope it'll be fixed in future releases. The error is non-fatal and can be safely ignored.
- "Out of Memory Error" when calculating Kramers-Kronig transform with a very large number of FFT-sampled points (very small step size in the loaded imaginary spectrum). Most users will not see this if they have a reasonable amount of RAM. There are ways to alleive this somewhat with more complex code that I won't bother to explore until someone tells me that it is a real problem. If you see this error, try the "double exponential" algorithm. The stubborn could also add more physical RAM or increase the virtual memory (increase size of swap partition in Linux or increase the size of the pagefile in Windows).
- Loading and saving is too simplistic. I'll work on it. Suggestions and code submissions welcome!
- Higher order interpolation could be made to work without oscillations at the Henke data absorption edges. I'll get to it eventually. Suggestions and code submissions welcome!
