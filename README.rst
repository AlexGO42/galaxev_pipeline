galaxev_pipeline
================

Python and C code for creating synthetic images using the GALAXEV
population synthesis code (Bruzual & Charlot 2003).

Brief description
-----------------

This code has three main stages:

- ``stellar_photometrics.py`` : For a given snapshot and a list of broadband
  filters, this calculates the apparent magnitude of a solar mass star
  as a function of metallicity and stellar age, which is stored as a 2D
  array. This includes effects from unresolved dust attenuation
  (e.g. Charlot & Fall 2000; Nelson et al. 2018, Section 3.2).
- The table created in the previous step is used to calculate the
  apparent magnitude of every stellar particle in a given snapshot.
- Once the magnitudes have been calculated, the code generates synthetic
  images at a given pixel scale, implementing an adaptive smoothing length
  (e.g. Torrey et al. 2015).

Requires numpy, matplotlib, scipy, h5py, astropy.

Compilation
-----------

Compile the external C module as:

.. code:: bash

    gcc -o adaptive_smoothing.so -shared -fPIC -O3 adaptive_smoothing.c

How to use
----------

First stage:
- Go to the BC03 ``src`` directory.
- Switch to C shell:
$ csh
- Convert relevant model files to ASCII:
$ ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_hr_stelib_m22_chab_ssp.ised
$ ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_hr_stelib_m32_chab_ssp.ised
...
$ ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_hr_stelib_m82_chab_ssp.ised
- Optionally, rest-frame magnitudes can be calculated as follows:
- Append filter names to ``filters.log`` and filter response functions to
  ``filterfrm.res`` (wavelengths are in angstroms).
- Specify filter names of interest in ``RF_COLORS.filters``
- Create binary filter file:
$ ./build_filterbin

Author
------

Vicente Rodriguez-Gomez (vrg [at] jhu.edu)

Licensing
---------

Licensed under a 3-Clause BSD License.
