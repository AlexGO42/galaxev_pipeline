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
from astropy.cosmology import FlatLambdaCDM

import illustris_python as il

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
                f.readline()
            # The next line has the wavelengths:
            line = f.readline()
            words = line.split()
            wavelengths = np.array(words[1:], dtype=np.float64)
            assert len(wavelengths) == int(words[0]) == num_wavelengths
            # The next 221 lines are spectra
            for j in range(num_stellar_ages):
                line = f.readline()
                cur_sed = np.array(line.split()[1:], dtype=np.float64)
                # There are some extra data points in each line, so we
                # truncate to the supposed number of wavelengths:
                datacube[k, j, :] = cur_sed[:num_wavelengths]
            # There are a few more lines after this (with 221 values per line),
            # which we don't need.
    print('Time: %g s.' % (time.time() - start))

    return datacube, metallicities, stellar_ages, wavelengths


def apply_cf00(datacube, stellar_ages, wavelengths):
    """
    Apply dust model from Charlot & Fall (2000).
    """
    t_BC = 1e7  # yr
    tau_BC = 1.0
    tau_ISM = 0.3
    n = -0.7
    lambda_0 = 5500.0  # angstroms

    tau_cube = np.zeros_like(datacube)
    tau_cube[:, stellar_ages <= t_BC, :] = tau_BC
    tau_cube[:, stellar_ages >  t_BC, :] = tau_ISM

    datacube *= np.exp(-tau_cube * (wavelengths[np.newaxis, np.newaxis, :] /
                                    lambda_0)**n)

    return datacube


def calculate_magnitudes(
        suite, use_z, use_cf00, filter_dir, filename_filters, filename_out,
        datacube, metallicities, stellar_ages, wavelengths):
    """
    Calculate AB magnitudes and store them in an HDF5 file,
    along with some attributes.
    """
    num_stellar_ages = len(stellar_ages)
    num_metallicities = len(metallicities)

    # Get luminosity distance (if redshift is too small, use 10 Mpc)
    if use_z < 2.5e-3:
        use_z = 0.0
        d_L = 10.0 * 3.086e22  # 10 Mpc in meters
        print('WARNING: Observation redshift is too small.',
              'Assuming that source is at 10 Mpc...')
    else:
        d_L = acosmo.luminosity_distance(use_z).value * 3.086e22  # meters

    # Apply Charlot & Fall (2000) model
    if use_cf00:
        datacube = apply_cf00(datacube, stellar_ages, wavelengths)

    # Shift rest-frame wavelengths to observer-frame; convert to meters
    wavelengths *= (1.0 + use_z) * 1e-10  # m

    # AB magnitude system in wavelength units
    FAB_nu = 3631.0 * 1e-26  # W/m^2/Hz
    c = 299792458.0  # m/s
    FAB_lambda = FAB_nu * c / wavelengths**2  # W/m^2/m

    # Convert rest-frame spectra from Lsun/angstrom to W/m
    datacube *= 3.826e26 * 1e10

    # Convert luminosity to flux in observer-frame.
    # Note that the (1+z) factor comes from the stretching
    # of dlambda (spreading of photons in wavelength).
    datacube *= 1.0 / ((4.0 * np.pi * d_L**2) * (1.0 + use_z))  # W/m^2/m

    # Read filter names
    with open(filename_filters, 'r') as f:
        filter_names = list(map(lambda s: s.strip('\n'), f.readlines()))

    # Open HDF5 for writing
    f = h5py.File(filename_out, 'w')
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

        # Iterate for every Z, age combination:
        for k in range(num_metallicities):
            for j in range(num_stellar_ages):
                F_lambda = datacube[k, j, :]  # observer-frame SED in W/m^2/m
                numerator = it.trapz(
                    F_lambda * wavelengths * R, x=wavelengths)
                magnitudes[k, j] = -2.5 * np.log10(numerator / denominator)

        # Create a dataset for the magnitudes and include some attributes.
        dset = f.create_dataset(filter_name, data=magnitudes)

        # Store the denominator from eq. (8) from the BC03 manual,
        # in case we need it later for calibration purposes:
        dset.attrs['denominator'] = float(denominator)

        # Also store the assumed redshift and the luminosity distance:
        dset.attrs['use_z'] = use_z
        dset.attrs['d_L'] = d_L

        print('Finished for filter %s.' % (filter_name))

    f.close()


if __name__ == '__main__':
    try:
        suite = sys.argv[1]
        basedir = sys.argv[2]
        writedir = sys.argv[3]
        bc03_model_dir = sys.argv[4]
        use_cf00 = bool(int(sys.argv[5]))
        snapnum = int(sys.argv[6])
        use_z = float(sys.argv[7])
        mock_set = sys.argv[8]  # 'hsc', etc.
    except:
        print('Arguments: suite basedir writedir bc03_model_dir',
              'use_cf00 snapnum use_z mock_set')
        sys.exit()

    # If True, use high resolution data (at 3 angstrom intervals) in the
    # wavelength range from 3200 to 9500 angstroms:
    high_resolution = False

    # Some additional directories and filenames
    suitedir = '%s/%s' % (writedir, suite)
    filter_dir = '%s/filter_curves' % (writedir,)
    filename_filters = '%s/filters.txt' % (writedir,)
    if use_cf00:
        filename_out = '%s/stellar_photometrics_cf00_%03d.hdf5' % (suitedir, snapnum)
    else:
        filename_out = '%s/stellar_photometrics_%03d.hdf5' % (suitedir, snapnum)

    # Make sure that write directories exist
    if not os.path.lexists(suitedir):
        os.makedirs(suitedir)

    # Cosmology
    if suite == 'IllustrisTNG':  # Planck 2015 XIII (Table 4, last column)
        acosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486)
    elif suite == 'Illustris':  # WMAP-7, Komatsu et al. 2011 (Table 1, v2)
        acosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
    else:
        raise Exception("Cosmology not specified.")

    # Load snapshot redshift from header
    header = il.groupcat.loadHeader(basedir, snapnum)
    z = header['Redshift']

    # If use_z is not specified, use the intrinsic snapshot redshift:
    if use_z == -1:
        use_z = z
        # However, if for some reason we are given the last snapshot,
        # set an arbitrary redshift:
        if ((suite == 'Illustris' and snapnum == 135) or (
             suite == 'IllustrisTNG' and snapnum == 99)):
            use_z = 0.0994018026302  # corresponds to snapnum_last - 8
            print('WARNING: use_z is too small. Setting use_z = %g.' % (use_z,))

    # Read BC03 model data
    datacube, metallicities, stellar_ages, wavelengths = read_bc03(
        bc03_model_dir, high_resolution)

    # Calculate magnitudes and store to HDF5 file
    calculate_magnitudes(
        suite, use_z, use_cf00, filter_dir, filename_filters, filename_out,
        datacube, metallicities, stellar_ages, wavelengths)
