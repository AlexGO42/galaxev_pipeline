import numpy as np
import time


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
    print("Time: %g s." % (time.time() - start))

    return datacube, metallicities, stellar_ages, wavelengths


def read_cb19():
    """
    T.B.D.
    """
    pass
