galaxev_pipeline
================

Python and C code for creating synthetic images using the GALAXEV
population synthesis code (Bruzual & Charlot 2003).

Brief description
-----------------

This code has two main stages:

- ``stellar_photometrics.py`` : For a given snapshot and a list of
  broadband filters, this program calculates the apparent magnitude of a
  single stellar population (SSP) (normalized to a mass of 1 Msun)
  as a function of metallicity and stellar age, which is stored as a
  2D array. Optionally, this calculation can include effects from an
  unresolved dust distribution using a Charlot & Fall (2000) model
  (e.g., Nelson et al. 2018, Section 3.2).
- ``create_images.py`` : Once the magnitude tables have been calculated,
  this program generates synthetic images for each broadband filter at a
  given pixel scale, implementing an adaptive smoothing length equal to
  the Nth nearest neighbor (usually N=16, e.g. Torrey et al. 2015).

Requires numpy, matplotlib, scipy, h5py, astropy.

How to use
----------

** First stage (apparent magnitudes): **

- Download GALAXEV code from http://www.bruzual.org/~gbruzual/bc03/Updated_version_2013/
  and compile it using C Shell (csh) and g77.
- After a successful installation, go to the BC03 ``src`` directory.
- Switch to C shell:

.. code:: bash

    csh

- Convert relevant model files to ASCII, e.g.:

.. code:: bash

    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m22_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m32_chab_ssp.ised
    ...
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m82_chab_ssp.ised

- Run ``stellar_photometrics.py`` with the appropriate input parameters.
  For convenience, an example batch script has been provided for this purpose:

.. code:: bash

    bash stellar_photometrics.sh

** Second stage (creating images): **

- Compile the external C module as:

.. code:: bash

    gcc -o adaptive_smoothing.so -shared -fPIC -O3 adaptive_smoothing.c

- The images are created by running ``create_images.py`` with
  appropriate input parameters. For convenience, an example batch script
  is also provided for this purpose:

.. code:: bash

    sbatch create_images.sbatch

Notes
-----

- By default, each image is a square with 20 stellar half-mass radii
  on each side (``num_rhalfs = 10``). In general, this is enough to
  include low surface brightness features in the outskirts and/or to
  reliably fit Sersic profiles.
- For reasons dating back to past projects (from which the current
  code has evolved), the code internally works in units of the stellar
  half-mass radius.
- As part of the input parameters, the code requires the angular momentum
  vectors of all galaxies in the snapshot of interest. This allows the
  code to (optionally) create face-on and edge-on projections. However,
  if the angular momentum vectors are not available, the code can still
  work with minimal modification (e.g., setting ``jvec = None``, etc.).

Author
------

Vicente Rodriguez-Gomez (vrg [at] jhu.edu)

Licensing
---------

Licensed under a 3-Clause BSD License.
