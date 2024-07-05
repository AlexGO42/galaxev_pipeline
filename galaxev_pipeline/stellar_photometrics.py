"""
Calculate magnitudes and other photometric properties for a set of
broadband filters based on models from GALAXEV.
"""

# Author: Vicente Rodriguez-Gomez <vrodgom.astro@gmail.com>
# Licensed under a 3-Clause BSD License.
import os
import sys
import time
import h5py

import numpy as np
import astropy.units as u
import scipy.integrate as it
import scipy.interpolate as ip
import astropy.constants as const

import illustris_python as il

from astropy.cosmology import LambdaCDM

from read_model_data import read_bc03, read_cb19


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
    tau_cube[:, stellar_ages > t_BC, :] = tau_ISM

    datacube *= np.exp(
        -tau_cube * (wavelengths[np.newaxis, np.newaxis, :] / lambda_0) ** n
    )

    return datacube


def calculate_magnitudes(
    suite,
    use_z,
    use_cf00,
    filter_dir,
    filename_filters,
    filename_out,
    datacube,
    metallicities,
    stellar_ages,
    wavelengths,
):
    """
    Calculate AB magnitudes and store them in an HDF5 file,
    along with some attributes.
    """
    num_stellar_ages = len(stellar_ages)
    num_metallicities = len(metallicities)

    # Get luminosity distance (if redshift is too small, use 10 Mpc)
    if use_z < 2.5e-3:
        use_z = 0.0
        d_L = (10.0 * u.Mpc).to(u.m)  # meters
        print(
            "WARNING: Observation redshift is too small.",
            "Assuming that source is at 10 Mpc...",
        )
    else:
        d_L = acosmo.luminosity_distance(use_z).to(u.m)

    # Apply Charlot & Fall (2000) model
    if use_cf00:
        datacube = apply_cf00(datacube, stellar_ages, wavelengths)

    # Shift rest-frame wavelengths to observer-frame; convert to meters
    wavelengths = (wavelengths * (1.0 + use_z)).to(u.m)

    # AB magnitude system in wavelength units
    FAB_nu = 3631.0 * u.Jy
    FAB_lambda = FAB_nu.to(
        u.W / u.m**2 / u.m, equivalencies=u.spectral_density(wav=wavelengths)
    )

    # Convert rest-frame spectra from Lsun/angstrom to W/m
    datacube = datacube.to(u.W / u.m)

    # Convert luminosity to flux in observer-frame.
    # Note that the (1+z) factor comes from the stretching
    # of dlambda (spreading of photons in wavelength).
    datacube *= 1.0 / ((4.0 * np.pi * d_L**2) * (1.0 + use_z))  # W/m^2/m

    # Read filter names
    with open(filename_filters, "r") as f:
        filter_names = list(map(lambda s: s.strip("\n"), f.readlines()))

    # Drop units (scipy's splines do not yet support astropy units)
    datacube = datacube.value
    wavelengths = wavelengths.value
    stellar_ages = stellar_ages.value
    FAB_lambda = FAB_lambda.value

    # Open HDF5 for writing
    f = h5py.File(filename_out, "w")
    f.create_dataset("metallicities", data=metallicities)
    f.create_dataset("stellar_ages", data=stellar_ages)

    # Iterate over filters
    for filter_name in filter_names:
        # Store magnitudes here:
        magnitudes = np.zeros(
            (num_metallicities, num_stellar_ages), dtype=np.float32
        )

        # Read filter response function
        filter_data = np.loadtxt(os.path.join(filter_dir, filter_name))
        filter_lambda, filter_response = filter_data.T
        filter_lambda *= 1e-10  # from AA to meters
        filter_interp = ip.interp1d(
            filter_lambda, filter_response, bounds_error=False, fill_value=0.0
        )

        # Apply eq. (8) from the BC03 manual to calculate the apparent
        # magnitude of the integrated photon flux collected by a detector
        # with filter response R(lambda).
        R = filter_interp(wavelengths)
        denominator = it.trapz(FAB_lambda * wavelengths * R, x=wavelengths)

        # Iterate for every Z, age combination:
        for k in range(num_metallicities):
            for j in range(num_stellar_ages):
                F_lambda = datacube[k, j, :]  # observer-frame SED in W/m^2/m
                numerator = it.trapz(F_lambda * wavelengths * R, x=wavelengths)
                magnitudes[k, j] = -2.5 * np.log10(numerator / denominator)

        # Create a dataset for the magnitudes and include some attributes.
        dset = f.create_dataset(filter_name, data=magnitudes)

        # Store the denominator from eq. (8) from the BC03 manual,
        # in case we need it later for calibration purposes:
        dset.attrs["denominator"] = float(denominator)

        # Also store the assumed redshift and the luminosity distance:
        dset.attrs["use_z"] = use_z
        dset.attrs["d_L"] = d_L.value

        print(f"Finished for filter {filter_name}.")

    f.close()


if __name__ == "__main__":
    try:
        suite = sys.argv[1]
        basedir = sys.argv[2]
        writedir = sys.argv[3]
        model_version = sys.argv[4]  # 'bc03', 'cb19'
        model_dir = sys.argv[5]
        use_cf00 = bool(int(sys.argv[6]))
        snapnum = int(sys.argv[7])
        use_z = float(sys.argv[8])
        mock_set = sys.argv[9]  # 'hsc', etc.
    except:
        print(
            "Arguments: suite basedir writedir model_version model_dir",
            "use_cf00 snapnum use_z mock_set",
        )
        sys.exit()

    # If True, use high resolution data (at 3 angstrom intervals) in the
    # wavelength range from 3200 to 9500 angstroms:
    high_resolution = False

    # Some additional directories and filenames
    suitedir = os.path.join(writedir, suite)
    filter_dir = os.path.join(writedir, "filter_curves")
    filename_filters = os.path.join(writedir, "filters.txt")
    if use_cf00:
        filename_out = os.path.join(
            suitedir, f"stellar_photometrics_cf00_{snapnum:03d}.hdf5"
        )
    else:
        filename_out = os.path.join(
            suitedir, f"stellar_photometrics_{snapnum:03d}.hdf5"
        )

    # Make sure that write directories exist
    if not os.path.lexists(suitedir):
        os.makedirs(suitedir)

    # Cosmology
    if suite == "IllustrisTNG":
        acosmo = LambdaCDM(
            name=suite,
            H0=67.74,
            Om0=0.3089,
            Ob0=0.0486,
            Ode0=0.6911,
            meta={
                "n_s": 0.9667,
                "sigma8": 0.8159,
                "reference": "Planck 2015 XIII (Table 4, last column)",
            },
        )
    elif suite == "Illustris":
        acosmo = LambdaCDM(
            name=suite,
            H0=70.4,
            Om0=0.2726,
            Ob0=0.0456,
            Ode0=0.7274,
            meta={
                "n_s": 0.963,
                "sigma8": 0.809,
                "reference": "WMAP-7, Komatsu et al. 2011 (Table 1, v2)",
            },
        )
    else:
        raise Exception("Cosmology not specified.")

    # Load snapshot redshift from header
    header = il.groupcat.loadHeader(basedir, snapnum)
    z = header["Redshift"]

    # If use_z is not specified, use the intrinsic snapshot redshift:
    if use_z == -1:
        use_z = z
        # However, if for some reason we are given the last snapshot,
        # set an arbitrary redshift:
        if (suite == "Illustris" and snapnum == 135) or (
            suite == "IllustrisTNG" and snapnum == 99
        ):
            use_z = 0.0994018026302  # corresponds to snapnum_last - 8
            print(f"WARNING: use_z is too small. Setting {use_z = :g}.")

    # Read model data
    if model_version == "bc03":
        datacube, metallicities, stellar_ages, wavelengths = read_bc03(
            model_dir, high_resolution
        )
    elif model_version == "cb19":
        datacube, metallicities, stellar_ages, wavelengths = read_cb19(
            model_dir
        )
    else:
        raise NotImplementedError(model_version)

    # Calculate magnitudes and store to HDF5 file
    calculate_magnitudes(
        suite,
        use_z,
        use_cf00,
        filter_dir,
        filename_filters,
        filename_out,
        datacube,
        metallicities,
        stellar_ages,
        wavelengths,
    )
