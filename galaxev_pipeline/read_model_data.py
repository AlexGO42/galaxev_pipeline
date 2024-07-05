import time
import numpy as np
import pandas as pd
import astropy.units as u

from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning

from numpy.lib import recfunctions as rfn


def read_bc03(bc03_model_dir, high_resolution=False):
    """
    Read single stellar population (SSP) model data from GALAXEV.
    Most parameters are hardcoded and might only work with the
    2013 version of GALAXEV.
    """
    if high_resolution:
        num_wavelengths = 6917
        model_filenames = [
            "bc2003_hr_stelib_m22_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m32_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m42_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m52_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m62_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m72_chab_ssp.ised_ASCII",
            "bc2003_hr_stelib_m82_chab_ssp.ised_ASCII",
        ]
    else:
        num_wavelengths = 1238
        model_filenames = [
            "bc2003_lr_BaSeL_m22_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m32_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m42_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m52_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m62_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m72_chab_ssp.ised_ASCII",
            "bc2003_lr_BaSeL_m82_chab_ssp.ised_ASCII",
        ]
    metallicities = np.array(
        [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1], dtype=np.float64
    )
    num_metallicities = len(metallicities)  # 7
    num_stellar_ages = 221
    datacube = np.zeros(
        (num_metallicities, num_stellar_ages, num_wavelengths),
        dtype=np.float64,
    )

    # Read and parse ASCII file
    print("Reading BC03 ASCII data...")
    start = time.time()
    for k in range(num_metallicities):
        with open("%s/%s" % (bc03_model_dir, model_filenames[k]), "r") as f:
            # The first line has the stellar ages
            line = f.readline()
            words = line.split()
            stellar_ages = np.array(words[1:], dtype=np.float64) * u.yr
            assert len(stellar_ages) == int(words[0]) == num_stellar_ages
            # Next 5 lines are not needed:
            for i in range(5):
                f.readline()
            # The next line has the wavelengths:
            line = f.readline()
            words = line.split()
            wavelengths = np.array(words[1:], dtype=np.float64) * u.angstrom
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
    datacube = datacube * (u.solLum / u.angstrom)
    print(f"Time: {time.time() - start:g} s.")

    return datacube, metallicities, stellar_ages, wavelengths


def read_cb19(cb19_model_dir):
    """
    Read single stellar population (SSP) model data from GALAXEV.
    Most parameters are hardcoded and might only work with the
    2019 version of GALAXEV.
    """
    num_wavelengths = 16902
    model_filenames = [
        "cb2019_z0001_chab_hr_xmilesi_ssp.fits",
        "cb2019_z0002_chab_hr_xmilesi_ssp.fits",
        "cb2019_z0005_chab_hr_xmilesi_ssp.fits",
        "cb2019_z001_chab_hr_xmilesi_ssp.fits",
        "cb2019_z002_chab_hr_xmilesi_ssp.fits",
        "cb2019_z004_chab_hr_xmilesi_ssp.fits",
        "cb2019_z006_chab_hr_xmilesi_ssp.fits",
        "cb2019_z008_chab_hr_xmilesi_ssp.fits",
        "cb2019_z010_chab_hr_xmilesi_ssp.fits",
        "cb2019_z014_chab_hr_xmilesi_ssp.fits",
        "cb2019_z017_chab_hr_xmilesi_ssp.fits",
        "cb2019_z020_chab_hr_xmilesi_ssp.fits",
        "cb2019_z030_chab_hr_xmilesi_ssp.fits",
        "cb2019_z040_chab_hr_xmilesi_ssp.fits",
        "cb2019_z060_chab_hr_xmilesi_ssp.fits",
    ]
    metallicities = np.array(
        [
            0.0001,
            0.0002,
            0.0005,
            0.001,
            0.002,
            0.004,
            0.006,
            0.008,
            0.010,
            0.014,
            0.017,
            0.020,
            0.030,
            0.040,
            0.060,
        ],
        dtype=np.float64,
    )
    num_metallicities = len(metallicities)
    num_stellar_ages = 221
    datacube = np.zeros(
        (num_metallicities, num_stellar_ages, num_wavelengths),
        dtype=np.float64,
    )
    # Read and parse FITS file
    print("Reading CB19 FITS data...")
    start = time.time()
    for k in range(num_metallicities):
        with fits.open(f"{cb19_model_dir}/{model_filenames[k]}") as hdulist:
            wav_sed_table = Table(hdulist[1].data)  # 16902R x 222C
            age_table = Table(hdulist[5].data)
            # First column contains the wavelengths; the rest are SEDs
            sed_cols = wav_sed_table.colnames[1:]
            wavelengths = wav_sed_table["Wavelength"].data * u.angstrom
            # Structured array is converted to an ordinary numpy array
            sed = rfn.structured_to_unstructured(
                wav_sed_table[sed_cols].as_array()
            )
            assert len(wavelengths) == num_wavelengths
            stellar_ages = age_table["age-yr"].data * u.yr
            assert len(stellar_ages) == num_stellar_ages
            datacube[k, ...] = sed.T
    datacube = datacube * (u.solLum / u.angstrom)
    print(f"Time: {time.time() - start:g} s.")

    return datacube, metallicities, stellar_ages, wavelengths
