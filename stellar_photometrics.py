import numpy as np
import h5py
import sys
import os
import time
import scipy.interpolate as ip
import scipy.integrate as it

from astropy.io import fits
import astropy.constants as const
import cosmology as cosmo

def read_bc03():
    """
    Read single stellar population (SSP) model data from GALAXEV.
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
            # but I'm not sure what those are.
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
        filter_dir = sys.argv[4]
        filename_filters = sys.argv[5]
        codedir = sys.argv[6]
        use_cf00 = bool(int(sys.argv[7]))
        use_z = float(sys.argv[8])
    except:
        print('Arguments: suite writedir bc03_model_dir filter_dir filename_filters codedir use_cf00 use_z')
        sys.exit()

    # If True, use high resolution data (at 3 angstrom intervals) in the
    # wavelength range from 3200 to 9500 angstroms:
    high_resolution = False

    # Make sure write directory exists
    if not os.path.lexists(writedir):
        os.makedirs(writedir)

    # Get luminosity distance (if redshift is too small, use 10 Mpc)
    if use_z < 2.5e-3:
        use_z = 0.0
        d_L = 10.0 * 3.086e22  # 10 Mpc in meters
        print('WARNING: Assuming that source is at 10 Mpc...')
    else:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_L = cosmo.luminosity_distance(use_z, params)  # meters

    # Read BC03 model data
    datacube, metallicities, stellar_ages, wavelengths = read_bc03()
    num_stellar_ages = len(stellar_ages)
    num_metallicities = len(metallicities)

    # Apply Charlot & Fall (2000) model
    if use_cf00:
        datacube = apply_cf00(datacube, stellar_ages, wavelengths)

    # Shift rest-frame wavelengths to observer-frame; convert to meters
    wavelengths *= (1.0 + use_z) * 1e-10  # m
    
    # AB magnitude system in wavelength units
    FAB_nu = 3631.0 * 1e-26  # W/m^2/Hz
    FAB_lambda = FAB_nu * const.c.value / (wavelengths)**2  # W/m^2/m

    # Convert rest-frame spectra from Lsun/angstrom to W/m
    datacube *= 3.826e26 * 1e10
    
    # Convert luminosity to flux
    datacube *= 1.0 / (4.0 * np.pi * d_L**2)  # W/m^2/m

    # Read filter names
    with open(filename_filters, 'r') as f:
        filter_names = list(map(lambda s: s.strip('\n'), f.readlines()))

    # Write output to this file
    if use_cf00:
        f = h5py.File('%s/stellar_photometrics_cf00.hdf5' % (writedir), 'w')
    else:
        f = h5py.File('%s/stellar_photometrics.hdf5' % (writedir), 'w')
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

        # Apply eq. (8) from the BC03 manual
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

        # Note that this denominator is proportional to the reference flux,
        # for each filter, which we store for convenience.
        # If the filter response is given as a capture cross-section in
        # m^2 counts/photon (e.g. Pan-STARRS), then the units are counts/s.
        # Otherwise, if the filter response is an adimensional quantum efficiency,
        # then the units are photons/s/m^2.
        zeropoint_phot = float(denominator) / (const.h.value * const.c.value)
        dset.attrs['zeropoint_phot'] = zeropoint_phot
        # We also store the assumed redshift and luminosity distance:
        dset.attrs['use_z'] = use_z
        dset.attrs['d_L'] = d_L

        print('Finished for filter %s.' % (filter_name))

    f.close()
