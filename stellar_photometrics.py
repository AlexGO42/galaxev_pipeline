"""
Calculate magnitudes and other photometric properties for a set of
broadband filters based on models from GALAXEV.
"""
# Author: Vicente Rodriguez-Gomez <vrodgom.astro@gmail.com>
# Licensed under a 3-Clause BSD License.
import numpy as np
import h5py
import sys
import os
import time
import scipy.interpolate as ip
import scipy.integrate as it

import astropy.constants as const
import cosmology as cosmo

# Constants
h = const.h.value  # J*s
c = const.c.value  # m/s


def read_bc03(bc03_model_dir, high_resolution=False):
    """
    Read single stellar population (SSP) model data from GALAXEV.
    Most parameters are hardcoded and might only work with the
    2013 version of GALAXEV.
    """
    if high_resolution:
        num_wavelengths = 6917
        model_filenames = [
            'bc2003_hr_stelib_m22_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m32_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m42_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m52_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m62_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m72_chab_ssp.ised_ASCII',
            'bc2003_hr_stelib_m82_chab_ssp.ised_ASCII',
        ]
    else:
        num_wavelengths = 1238
        model_filenames = [
            'bc2003_lr_BaSeL_m22_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m32_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m42_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m52_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m62_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m72_chab_ssp.ised_ASCII',
            'bc2003_lr_BaSeL_m82_chab_ssp.ised_ASCII',
        ]
    metallicities = np.array([
        0.0001,
        0.0004,
        0.004,
        0.008,
        0.02,
        0.05,
        0.1
    ], dtype=np.float64)
    num_metallicities = len(metallicities)  # 7
    num_stellar_ages = 221
    datacube = np.zeros(
        (num_metallicities, num_stellar_ages, num_wavelengths), dtype=np.float64)

    # Read and parse ASCII file
    print('Reading BC03 ASCII data...')
    start = time.time()
    for k in range(num_metallicities):
        with open('%s/%s' % (bc03_model_dir, model_filenames[k]), 'r') as f:
            # The first line has the stellar ages
            line = f.readline()
            words = line.split()
            stellar_ages = np.array(words[1:], dtype=np.float64)
            assert len(stellar_ages) == int(words[0]) == num_stellar_ages
            # Next 5 lines are not needed:
            for i in range(5):
                line = f.readline()
            # The next line has the wavelengths:
            line = f.readline()
            words = line.split()
            wavelengths = np.array(words[1:], dtype=np.float64)
            assert len(wavelengths) == int(words[0]) == num_wavelengths
            # The next 221 lines are spectra
            for j in range(num_stellar_ages):
                line = f.readline()
                cur_sed = np.array(line.split()[1:], dtype=np.float64)
                # There are some extra data points in each line (at least
                # when using the high-resolution stelib data), so we
                # truncate to the supposed number of wavelengths:
                datacube[k, j, :] = cur_sed[:num_wavelengths]
            # There are a few more lines after this (with 221 values per line),
            # which we don't need.
    print('Time: %g s.' % (time.time() - start))

    return datacube, metallicities, stellar_ages, wavelengths


def apply_cf00(datacube, stellar_ages, wavelengths):
    t_BC = 1e7  # yr
    tau_BC = 1.0
    tau_ISM = 0.3
    n = -0.7
    lambda_0 = 5500.0  # angstroms

    tau_cube = np.zeros_like(datacube)
    tau_cube[:, stellar_ages <= t_BC, :] = tau_BC
    tau_cube[:, stellar_ages >  t_BC, :] = tau_ISM

    datacube *= np.exp(-tau_cube*(wavelengths[np.newaxis,np.newaxis,:]/lambda_0)**n)

    return datacube


if __name__ == '__main__':

    try:
        suite = sys.argv[1]
        writedir = sys.argv[2]
        bc03_model_dir = sys.argv[3]
        use_cf00 = bool(int(sys.argv[4]))
        snapnum = int(sys.argv[5])
        use_z = float(sys.argv[6])
        mock_set = sys.argv[7]  # 'hsc', etc.
    except:
        print('Arguments: suite writedir bc03_model_dir',
              'use_cf00 snapnum use_z mock_set')
        sys.exit()

    # Some additional directories and filenames
    suitedir = '%s/%s' % (writedir, suite)
    filter_dir = '%s/filter_curves' % (writedir,)
    filename_filters = '%s/filters.txt' % (writedir,)

    # If True, use high resolution data (at 3 angstrom intervals) in the
    # wavelength range from 3200 to 9500 angstroms:
    high_resolution = False

    # Make sure that write directories exist
    if not os.path.lexists(suitedir):
        os.makedirs(suitedir)

    # Get luminosity distance (if redshift is too small, use 10 Mpc)
    if use_z < 2.5e-3:
        use_z = 0.0
        d_L = 10.0 * 3.086e22  # 10 Mpc in meters
        print('WARNING: Assuming that source is at 10 Mpc...')
    else:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_L = cosmo.luminosity_distance(use_z, params)  # meters

    # Read BC03 model data
    datacube, metallicities, stellar_ages, wavelengths = read_bc03(
        bc03_model_dir, high_resolution)
    num_stellar_ages = len(stellar_ages)
    num_metallicities = len(metallicities)

    # Apply Charlot & Fall (2000) model
    if use_cf00:
        datacube = apply_cf00(datacube, stellar_ages, wavelengths)

    # Shift rest-frame wavelengths to observer-frame; convert to meters
    wavelengths *= (1.0 + use_z) * 1e-10  # m

    # AB magnitude system in wavelength units
    FAB_nu = 3631.0 * 1e-26  # W/m^2/Hz
    FAB_lambda = FAB_nu * c / (wavelengths)**2  # W/m^2/m

    # Convert rest-frame spectra from Lsun/angstrom to W/m
    datacube *= 3.826e26 * 1e10

    # Convert luminosity to flux
    datacube *= 1.0 / (4.0 * np.pi * d_L**2)  # W/m^2/m

    # Read filter names
    with open(filename_filters, 'r') as f:
        filter_names = list(map(lambda s: s.strip('\n'), f.readlines()))

    # Write output to this file
    if use_cf00:
        filename = '%s/stellar_photometrics_cf00_%03d.hdf5' % (suitedir, snapnum)
    else:
        filename = '%s/stellar_photometrics_%03d.hdf5' % (suitedir, snapnum)
    f = h5py.File(filename, 'w')
    f.create_dataset('metallicities', data=metallicities)
    f.create_dataset('stellar_ages', data=stellar_ages)

    # Iterate over filters
    for filter_name in filter_names:
        # Store magnitudes here:
        magnitudes = np.zeros((num_metallicities, num_stellar_ages), dtype=np.float32)

        # Read filter response function
        filter_data = np.loadtxt('%s/%s' % (filter_dir, filter_name))
        filter_lambda, filter_response = filter_data.T
        filter_lambda *= 1e-10  # to meters
        filter_interp = ip.interp1d(filter_lambda, filter_response, bounds_error=False, fill_value=0.0)

        # Apply eq. (8) from the BC03 manual to calculate the apparent
        # magnitude of the integrated photon flux collected by a detector
        # with filter response R(lambda).
        R = filter_interp(wavelengths)
        denominator = it.trapz(
            FAB_lambda * wavelengths * R, x=wavelengths)

        # For now, iterate for every Z, age combination:
        for k in range(num_metallicities):
            for j in range(num_stellar_ages):
                F_lambda = datacube[k, j, :]
                # (1+z) factor comes from the stretching of dlambda
                # (spreading of photons in wavelength)
                numerator = it.trapz(
                    F_lambda / (1.0 + use_z) * wavelengths * R, x=wavelengths)
                magnitudes[k, j] = -2.5 * np.log10(numerator / denominator)

        dset = f.create_dataset(filter_name, data=magnitudes)

        # Now that we have the magnitudes (which do not require knowing
        # the units, if any, of the transmission curve), we note that
        # the denominator is essentially the wavelength-integrated photon
        # flux that corresponds to a magnitude of zero. We will eventually
        # need this for calibration purposes. Our goal is to create mock
        # images that have the same units as the corresponding real images,
        # which depends on the instrument or survey. In order to do this,
        # below we calculate a quantity called "fluxmag0", which gives the
        # "flux" -- in image units -- that corresponds to a magnitude of zero.
        # This requires knowing the units, if any, of the transmission curve
        # (e.g. a capture cross-section in m^2 electrons/photon in the case
        # of Pan-STARRS and GALEX) or, alternatively, applying appropriate
        # zeropoints (MAG = -2.5 * log10(data) + ZP) for each filter.
        if mock_set == 'pogs':
            # In Pan-STARRS, the filter response is given as a capture
            # cross-section in m^2 electrons/photon (Tonry et al. 2012).
            # This is ideal because the "denominator" calculated above
            # divided by hc is already the number of electrons/s
            # registered by the CCD:
            fluxmag0 = float(denominator) / (h*c)
        elif mock_set.startswith('sdss'):
            # SDSS filter curves are expressed as an adimensional quantum
            # efficiency (electrons per photon), so we need to
            # multiply by the area of the 2.5 m primary mirror:
            area = np.pi * (2.5/2.0)**2  # m^2
            fluxmag0 = float(denominator) / (h*c) * area
        elif mock_set.startswith('kids'):
            # The science images from the Kilo-Degree Survey (KiDS) DR4
            # are in units of ADU/s and appear to have been calibrated
            # so that the zeropoint is zero (according to the r-band
            # image headers; I have not checked other bands). Therefore:
            fluxmag0 = 1.0
        elif mock_set == 'galex':
            # Thankfully, GALEX filter curves are already expressed as
            # an effective area (just like Pan-STARRS). We just need to
            # convert the units from cm^2 to m^2:
            fluxmag0 = float(denominator) / (h*c) * 1e-4
        elif mock_set == 'candels_acs' or mock_set == 'candels_wfc3':
            # Although HST doesn't use AB magnitudes, we use them here
            # for consistency with the rest of the code. The final images
            # have units of counts/s, so the magnitude system used here
            # is unimportant. The calculations below assume that the units
            # of the transmission curve are electrons/photon. These assumptions
            # were verified by comparing with the "PHOTFNU" and "PHOTFLAM"
            # attributes found in actual FITS headers from CANDELS.
            area = np.pi * (2.4/2.0)**2  # m^2
            fluxmag0 = float(denominator) / (h*c) * area
        elif mock_set == 'hsc':
            # Like SDSS, Hyper Suprime-Cam filter curves are also
            # expressed as a quantum efficiency (electrons per photon).
            # I was originally assuming a collecting area with a diameter = 8.2 m,
            # which yielded the following zero-points (ZP = 2.5 * log10(fluxmag0)),
            # aka the magnitude that corresponds to 1 electron/s:
            # ZP(hsc_g) = 29.147
            # ZP(hsc_r) = 29.222
            # ZP(hsc_i) = 28.807
            # ZP(hsc_z) = 27.902
            # ZP(hsc_y) = 27.525
            # However, I now use the slightly different (and hopefully
            # more accurate) values from the official HSC website
            # (https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html):
            ZP = {
                'subaru/hsc_g': 28.8,
                'subaru/hsc_r': 28.8,  # for the old r-filter, not the new r2-filter
                'subaru/hsc_i': 28.5,  # for the old i-filter, not the new i2-filter
                'subaru/hsc_z': 27.5,
                'subaru/hsc_y': 27.2,
            }
            fluxmag0 = 10.0**(ZP[filter_name] / 2.5)
        else:
            print('mock_set not understood.')
            sys.exit()

        dset.attrs['fluxmag0'] = fluxmag0  # in electrons/s
        # We also store the assumed redshift and luminosity distance:
        dset.attrs['use_z'] = use_z
        dset.attrs['d_L'] = d_L

        print('Finished for filter %s.' % (filter_name))

    f.close()
