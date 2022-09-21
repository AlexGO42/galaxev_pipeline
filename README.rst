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
This procedure yields "idealized" images (i.e., as they would be observed
by a perfect telescope with a point-like PSF and negligible noise),
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
snapshot 78 (*z* = 0.2977) of the TNG100 simulation with settings that
mimic the Hyper Suprime-Cam (HSC) on the Subaru Telescope. We will
create images for the the HSC g,r,i,z,Y filters. To this
end, let us create a directory to store all the relevant data:

.. code:: bash

    mkdir hsc

In the following sections, some shell scripts are provided as examples
(also for the current HSC example). In general, however, the programs
``stellar_photometrics.py`` and ``create_images.py`` can be run for any
combination of filters, pixel scales, and redshifts, and accept various
command-line parameters to customize the synthetic images.

Obtaining filter curves
-----------------------

If the filters of interest (e.g. HSC g,r,i,z,Y) are listed in the
`SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_,
then they can be retrieved automatically via 
`Astroquery <https://astroquery.readthedocs.io/en/latest/svo_fps/svo_fps.html>`_.
For convenience, a batch script (please modify as needed) is provided
to carry out this task:

.. code:: bash

    bash get_filter_curves.sh

This writes the filter IDs to a file ``hsc/filters.txt`` and stores
the filter curves in the folder ``hsc/filter_curves``. Each
filter curve file is in ASCII format and includes two columns:
the wavelength in angstroms (AA) and the transmission values.
Alternatively, if the desired filters are not found in the
SVO Filter Profile Service, they can be included manually (without
Astroquery) by listing the filter names in the text file
``hsc/filters.txt`` and including their transmission curves in the folder
``hsc/filter_curves``, following the same convention (ASCII format,
two columns, wavelength in angstroms).

Note that the transmission curves used in this example (HSC g,r,i,z,Y)
already include the contribution from the instrument and atmosphere
(Filter + Instrument + Atmosphere).

Calculating magnitudes
----------------------

The ``stellar_photometrics.py`` program
calculates the apparent (observer-frame) magnitude of an SSP
(normalized to a mass of 1 Msun) as a function of metallicity and
stellar age, which is stored as a grid (2D array) that can be
interpolated to obtain a magnitude for any metallicity and stellar age.
Optionally, this calculation can include effects from an
unresolved dust distribution using a Charlot & Fall (2000) model.

For convenience, an example batch script (please modify as needed)
has been provided to guide the user on how to run
``stellar_photometrics.py`` with appropriate input parameters:

.. code:: bash

    bash stellar_photometrics.sh

This generates a file called ``stellar_photometrics_078.hdf5``
(or ``stellar_photometrics_cf00_078.hdf5``, if the dust prescription
from Charlot & Fall 2000 is included) that contains a 2D array
with the apparent (observer-frame) magnitudes for each filter, assuming
that the observed source is located at a redshift *z* = 0.2977
(snapshot 78 in IllustrisTNG).

Creating the images
-------------------

Now that we have the precalculated magnitude tables, the program
``create_images.py`` can be used to generate synthetic images of
simulated galaxies for the chosen broadband filters, implementing an
adaptive smoothing scheme for the stellar particles with a smoothing
scale given by the distance to the Nth nearest neighbor (usually N=32).

The image generation stage requires knowing the pixel scale of the
instrument (in arcsec) and setting a few other parameters, such as the
field of view, the projection (e.g. face-on, edge-on, or aligned with the
axes of the simulation box), and whether or not to include neighboring
galaxies (within the same parent halo).

For convenience, an example batch script is once again provided for
this purpose, which can be copied and modified as needed:

.. code:: bash

    bash create_images.sh

Note that ``create_images.py`` is able to create many images in parallel
using MPI, and that the script ``create_images.sh`` can be easily modified
for submission to a job scheduler (e.g. SLURM) in a computer cluster.

Running the above script (after setting the correct directories, etc.)
creates *idealized* synthetic images for the most massive galaxies
(Mstar > 10^12 Msun, including the intracluster light) from snapshot 78
(z ~ 0.3) of TNG100, showing stars from neighboring galaxies as well,
with a fixed image size of 224x224 pixels and HSC settings
(pixel scale and filters). The resulting idealized images have units of "maggies"
(following SDSS nomenclature), which means that the zero-point is
exactly zero and that magnitudes can be calculated simply as
MAG = -2.5 * log10(DATA). The output is stored in FITS files, each
with a main HDU object in which the different "layers" correspond
to the filters of interest. The figure below shows RGB composite images
(for the HSC *i,r,g* bands, respectively), generated with the example
script ``extra/view_rgb_composites.py``, for the first nine objects
considered in this example.

Applying realism
================

The goal of the GALAXEV pipeline is to create "idealized" synthetic images,
i.e., as seen by an instrument with a point-like point spread function (PSF)
and without any sources of noise (infinite signal-to-noise ratio). These
idealized images can then be further processed in order to make them
"realistic", which facilitates robust comparisons to observations, as done in 
`Rodriguez-Gomez et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4140R>`_
and other works.

Since the "realism" stage largely depends on the set of observations
that one wishes to compare against, it would be difficult to implement
(and maintain) a unified solution that works for all instruments and
surveys. Therefore, the realism stage is ultimately left to the end user.
However, an example is provided in ``extra/apply_realism.py`` that
applies a moderate degree of realism to the HSC *i*-band images of the
*z* ~ 0.3 objects considered in the example above.

Briefly, the procedure in ``extra/apply_realism.py`` consists in
(i) convolving with a PSF, (ii) applying shot (Poisson) noise, and
(iii) adding uniform Gaussian background noise, adopting settings
(seeing, sensitivity, zero-point, etc.) that match those of
real HSC *i*-band images. In step (i), the PSF is assumed to be
a simple 2D Gaussian, which might not be accurate for some
applications. In step (iii), instead of adding uniform
background noise, the idealized images could be inserted into
real HSC backgrounds, which would achieve a higher level of
realism, but is not done here for the sake of simplicity.
The figure below shows the results of applying steps (i), (ii),
and (iii) to the HSC *i*-band images of the *z* ~ 0.3 objects
from the previous section.

Author
======

Vicente Rodriguez-Gomez (vrodgom.astro@gmail.com)

Citing
======

This code is fully described in
`Rodriguez-Gomez et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4140R>`_.

Licensing
=========

Licensed under a 3-Clause BSD License.
