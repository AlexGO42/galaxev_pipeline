================
galaxev_pipeline
================

Python and C code for creating synthetic images using the GALAXEV
population synthesis code
(`Bruzual & Charlot 2003 <https://ui.adsabs.harvard.edu/abs/2003MNRAS.344.1000B>`_).

Overview
========

This repository contains the code used to generate the synthetic
images of Illustris and IllustrisTNG galaxies presented in
`Rodriguez-Gomez et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4140R>`_.
This so-called "GALAXEV pipeline" complements the "SKIRT pipeline"
also described in Rodriguez-Gomez et al. (2019), which additionally
takes into account dust radiative transfer using SKIRT and employs
MAPPINGS-III models for starbursting regions. However,
both pipelines produce essentially identical images for galaxies
with low fractions of star-forming gas (fig. 2 from Rodriguez-Gomez
et al. 2019), and the GALAXEV pipeline has the advantage of being
substantially faster and easier to use. As a low-cost alternative
to full dust radiative transfer, the current pipeline can optionally
include the effects of a spatially unresolved dust distribution
according to the model of
`Charlot & Fall (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...539..718C>`_.

The main idea of this code is the following. Every stellar particle in the
simulation can be considered to be a single stellar population (SSP),
whose spectral energy distribution (SED) can be modeled with GALAXEV
(for performance reasons, we do not actually deal with the spectra of
individual stellar particles, but the basic idea is similar).
So, for a given set of broadband filters, we can calculate the magnitude
(in the frame of the observer) of each stellar particle and its associated
flux. We then add the fluxes from all the stellar particles to create a final
image (with a given pixel scale), but instead of treating the particles as
point-like, their light contribution is convolved with a kernel with a
spatially varying smoothing scale, similar to the methodology described in
`Torrey et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015MNRAS.447.2753T>`_.
This procedure yields "idealized" images (i.e. as they would be observed
by a perfect telescope with a point-like PSF and zero background noise),
which can be further processed to create "realistic" images that can be
directly compared to observations (see Rodriguez-Gomez et al. 2019,
for more details).

Setup
=====

Requirements
------------

**Python packages:** numpy, scipy, h5py, astropy, illustris_python

**Compilers:** C (for the adaptive smoothing module), Fortran (for GALAXEV)

**Optional:** mpi4py (for parallel execution), astroquery
(for retrieving filter curves)

Obtaining GALAXEV spectra
-------------------------

The first step is to download the GALAXEV source code, which
can be done with the `csh` shell as follows (see also section 3 of the
`GALAXEV manual <http://www.bruzual.org/bc03/doc/bc03.pdf>`_):

.. code:: bash

    csh
    wget http://www.bruzual.org/bc03/src.tgz
    tar -zxf src.tgz

We must also download some GALAXEV models. For our purposes, the
low-resolution version (BaSeL spectral library) of the models computed
with the Padova 1994 evolutionary tracks and a Chabrier (2003)
initial mass function are appropriate:

.. code:: bash

    wget http://www.bruzual.org/bc03/Updated_version_2013/bc03.models.padova_1994_chabrier_imf.tar.gz
    tar -zxf bc03.models.padova_1994_chabrier_imf.tar.gz

A ``bc03`` folder should have appeared, which contains a ``src``
subdirectory with the source files as well as a
``Padova1994/chabrier`` subdirectory with the model data.

In order to save disk space, GALAXEV spectra are usually stored as
binary files, which we need to convert to text (ASCII). To do this,
we first compile GALAXEV (this requires a Fortran compiler):

.. code:: bash

    cd /path/to/bc03/src
    make all

We then convert the relevant model files to ASCII using the
``ascii_ised`` program:

.. code:: bash

    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m22_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m32_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m42_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m52_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m62_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m72_chab_ssp.ised
    ./ascii_ised /path/to/bc03/Padova1994/chabrier/bc2003_lr_BaSeL_m82_chab_ssp.ised

This produces similarly named files but with an ``*.ised_ASCII`` extension,
which can be read as simple text. Finally, exit the csh shell:

.. code:: bash

    exit

Compiling the adaptive smoothing module
---------------------------------------

An external C module for speeding up the adaptive smoothing
calculations needs to be compiled as follows:

.. code:: bash

    gcc -o adaptive_smoothing.so -shared -fPIC -O3 adaptive_smoothing.c

How to use
==========

Now that we have the GALAXEV models, we can define the instrument
for which we wish to generate synthetic images. In particular,
we need to provide a set of broadband filters and an
appropriate pixel scale. We will also need to set the redshift
at which the sources are assumed to be observed, which is
typically (but not necessarily) equal to the redshift of the
simulation snapshot.

For concreteness, let us create synthetic images of galaxies from
snapshot 91 (z = 0.0994) of the TNG50 simulation with settings that
mimic the Hyper Suprime-Cam (HSC) on the Subaru Telescope. To this
end, let us create a directory to store all the relevant data:

.. code:: bash

    mkdir hsc

Obtaining filter curves
-----------------------

Filter transmission curves for numerous instruments can be obtained from the
`SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_
or from the instrument/survey websites. In principle, this part of
the process could be automatized with
`Astroquery <https://astroquery.readthedocs.io/en/latest/svo_fps/svo_fps.html>`_,
which would be very nice, but at the time of this writing there seem to be
`connectivity issues <https://github.com/astropy/astroquery/issues/2508>`_,
so instead we manually download the Subaru/HSC g,r,i,z,Y filter curves
(in ASCII format) from the SVO Filter Profile Service and and place them in
a new `filter_curves` subdirectory:

.. code:: bash

    cd hsc
    mkdir filter_curves  # store filter curves here

This text files contain two columns, where the first one should consist
of the wavelengths in angstroms.

We also create a text file called ``filters.txt`` where we specify the
filenames of the filter curves (I renamed them slightly), one per line:

.. code:: bash

    hsc_g
    hsc_r
    hsc_i
    hsc_z
    hsc_y

Note that these transmission curves already include the contribution from the
instrument and atmosphere (Filter + Instrument + Atmosphere).

Calculating magnitudes
----------------------

Although GALAXEV provides Fortran programs for calculating the magnitudes
of SSPs for a given filter response function, we instead perform these
computations in Python so that they integrate better with the rest of
the code.

This is done by the ``stellar_photometrics.py`` program, which
calculates the apparent (observer-frame) magnitude of an SSP
(normalized to a mass of 1 Msun) as a function of metallicity and
stellar age, which is stored as a grid (2D array) that can be
interpolated to obtain a magnitude for any metallicity and stellar age.
Optionally, this calculation can include effects from an
unresolved dust distribution using a Charlot & Fall (2000) model.

For convenience, an example batch script has been provided in order to
run ``stellar_photometrics.py`` with the appropriate input parameters:

.. code:: bash

    bash stellar_photometrics.sh

This will generate a file called ``stellar_photometrics_091.hdf5``
(and/or ``stellar_photometrics_cf00_091.hdf5``, if the dust prescription
from Charlot & Fall 2000 is included) that contains a 2D array
with the observer-frame magnitudes for each filter.

---------- Work in progress... ----------

Creating the images
-------------------

Now that we have precalculated all the magnitudes,


** Second stage (creating images): **


- ``create_images.py`` : Once the magnitude tables have been calculated,
  this program generates synthetic images for each broadband filter at a
  given pixel scale, implementing an adaptive smoothing length equal to
  the Nth nearest neighbor (usually N=16, e.g. Torrey et al. 2015).



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

Vicente Rodriguez-Gomez (vrodgom.astro@gmail.com)

Citing
------

This code is fully described in
`Rodriguez-Gomez et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4140R>`_.

Licensing
---------

Licensed under a 3-Clause BSD License.
