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
    model_filenames = [
        'bc2003_hr_stelib_m22_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m32_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m42_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m52_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m62_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m72_chab_ssp.ised_ASCII',
        'bc2003_hr_stelib_m82_chab_ssp.ised_ASCII',
    ]
    metallicities = np.array([
        0.0001,
        0.0004,
        0.004,
        0.008,
        0.02,
        0.05,
        0.1,
    ], dtype=np.float32)
    num_metallicities = len(metallicities)  # 7
    num_stellar_ages = 221
    num_wavelengths_hr = 6917  # original/high resolution
    datacube_hr = np.zeros(
        (num_metallicities, num_stellar_ages, num_wavelengths_hr), dtype=np.float64)
    
    print('Reading BC03 ASCII data...')
    start = time.time()
    for k in range(num_metallicities):
        with open('%s/%s' % (bc03_model_dir, model_filenames[k]), 'r') as f:
            # The first line has the stellar ages
            line = f.readline()
            stellar_ages = np.array(line.split()[1:], dtype=np.float32)
            assert len(stellar_ages) == num_stellar_ages
            # Next 5 lines are not needed:
            for i in range(5):
                line = f.readline()
            # The next line has the wavelengths:
            line = f.readline()
            wavelengths_hr = np.array(line.split()[1:], dtype=np.float64)
            assert len(wavelengths_hr) == num_wavelengths_hr
            # The next 221 lines are spectra
            for j in range(num_stellar_ages):
                line = f.readline()
                cur_sed = np.array(line.split()[1:], dtype=np.float32)
                # There are some extra data points in each line:
                datacube_hr[k, j, :] = cur_sed[:num_wavelengths_hr]
            # There are a few more lines after this (with 221 values per line),
            # but I'm not sure what those are.
    print('Time: %g s.' % (time.time() - start))

    if degrade_resolution:
        # Degrade resolution just to make it similar to my SKIRT runs
        print('Degrading resolution...')
        start = time.time()
        datacube = np.zeros(
            (num_metallicities, num_stellar_ages, num_wavelengths), dtype=np.float64)
        for i in range(num_wavelengths):
            i_hr = np.searchsorted(wavelengths_hr, wavelengths[i]) - 1
            if (i_hr < 0) or (i_hr + 1 > num_wavelengths_hr - 1):
                raise Exception('Out of wavelength range.')
            fraction = (wavelengths[i] - wavelengths_hr[i_hr]) / (wavelengths_hr[i_hr+1] - wavelengths_hr[i_hr])
            assert (fraction >= 0.0) & (fraction <= 1.0)
            datacube[:, :, i] = (1.0 - fraction) * datacube_hr[:, :, i_hr] + fraction * datacube_hr[:, :, i_hr+1]
        print('Time: %g s.' % (time.time() - start))

        return datacube, stellar_ages, metallicities

    else:
        return datacube_hr, metallicities, stellar_ages, wavelengths_hr

if __name__ == '__main__':
    
    try:
        suite = sys.argv[1]
        snapnum = int(sys.argv[2])
        writedir = sys.argv[3]
        bc03_model_dir = sys.argv[4]
        filter_dir = sys.argv[5]
        filename_filters = sys.argv[6]
    except:
        print('Arguments: suite snapnum writedir bc03_model_dir filter_dir filename_filters')
        sys.exit()

    # ~ suite = 'IllustrisTNG'
    # ~ snapnum = 99
    # ~ writedir = '/n/ghernquist/vrodrigu/SyntheticImages/galaxev/IllustrisTNG/stellar_photometrics'
    # ~ bc03_model_dir = '/n/home10/vrodrigu/galaxev_code/bc03/Padova1994/chabrier'
    # ~ filter_dir = '/n/home10/vrodrigu/SyntheticImages/broadband_filters'
    # ~ filename_filters = '/n/home10/vrodrigu/SyntheticImages/broadband_filters/panstarrs.txt'

    degrade_resolution = False
    
    # Get redshift for given snapshot
    redshifts_all = np.loadtxt('Redshifts%s.txt' % (suite))
    z = redshifts_all[snapnum]

    # Make sure write directory exists
    if not os.path.lexists(writedir):
        os.makedirs(writedir)

    # Get luminosity distance (if redshift too small, use 10 Mpc)
    if z < 2.5e-3:
        d_L = 10.0 * 3.086e22  # 10 Mpc in meters
        print('WARNING: assuming that source is at 10 Mpc.')
    else:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_L = cosmo.luminosity_distance(z, params)  # meters

    if degrade_resolution:
        print('WARNING: degrading BC03 spectral resolution...')
        min_wavelength = 3000.0  # angstroms
        max_wavelength = 10000.0  # angstroms
        num_wavelengths = 100
        wavelengths = np.logspace(
            np.log10(min_wavelength), np.log10(max_wavelength), num_wavelengths)

    # Read BC03 model data
    datacube, metallicities, stellar_ages, wavelengths = read_bc03()
    num_stellar_ages = len(stellar_ages)
    num_metallicities = len(metallicities)
    
    # Shift rest-frame wavelengths to observer-frame; convert to meters
    wavelengths *= (1.0 + z) * 1e-10  # m
    
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
    f = h5py.File('%s/stellar_photometrics_%03d.hdf5' % (writedir, snapnum), 'w')
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
                    F_lambda / (1.0 + z) * wavelengths * R, x=wavelengths)
                magnitudes[k, j] = -2.5 * np.log10(numerator / denominator)

        dset = f.create_dataset(filter_name, data=magnitudes)

        # Note that this denominator is proportional to the reference flux
        # in counts/s (assuming that the filter response is given in m^2 counts/photon; otherwise ),
        # which we store for convenience.
        # If the filter response is given as a capture cross-section in
        # m^2 counts/photon (e.g. Pan-STARRS), then the units are photons/s.
        # Otherwise, if the filter response is an adimensional quantum efficiency,
        # then the units are photons/s/m^2.
        zeropoint_phot = float(denominator) / (const.h.value * const.c.value)  # photons/s
        dset.attrs['zeropoint_phot'] = zeropoint_phot

        print('Finished for filter %s.' % (filter_name))

    f.close()
