import numpy as np
import h5py
import os
import sys
import scipy.interpolate as ip
from scipy.spatial import cKDTree
from astropy.io import fits
from multiprocessing import Process, Queue, Pool
import time
import ctypes

import illustris_python as il
import cosmology as cosmo

parttype_stars = 4

def get_fluxes(initial_masses_Msun, metallicities, stellar_ages_yr, filter_name):
    """
    Return "fluxes" in electron counts per second (as collected by the CCD),
    using the BC03 table generated by stellar_photometrics.py.
    """
    if use_cf00:
        filename = '%s/stellar_photometrics_cf00.hdf5' % (stellar_photometrics_dir)
    else:
        filename = '%s/stellar_photometrics.hdf5' % (stellar_photometrics_dir)

    with h5py.File(filename, 'r') as f:
        bc03_metallicities = f['metallicities'][:]
        bc03_stellar_ages = f['stellar_ages'][:]
        bc03_magnitudes = f[filter_name][:]
        fluxmag0 = f[filter_name].attrs['fluxmag0']

    spline = ip.RectBivariateSpline(
        bc03_metallicities, bc03_stellar_ages, bc03_magnitudes, kx=1, ky=1, s=0)

    magnitudes = spline.ev(metallicities, stellar_ages_yr) - 2.5 * np.log10(initial_masses_Msun)

    # Convert apparent magnitudes to "fluxes" in counts/s. For more details,
    # see computation of "fluxmag0" in stellar_photometrics.py.
    fluxes = 10.0**(-2.0/5.0*magnitudes) * fluxmag0

    return fluxes


def transform(x, jvec, proj_kind='xy'):
    """
    Return a projection of the particle positions. In all cases we only
    return the first two coordinates (x, y).

    Parameters
    ----------
    x : array-like
        2-dimensional array (Nx3) with the particle positions.
    jvec : array-like
        Direction of the galaxy's total stellar angular momentum.
    proj_kind : str, optional
        Specify which kind of projection. The possible values are:

        'yz' : Project onto the yz plane.
        'zx' : Project onto the zx plane.
        'xy' : Project onto the xy plane.
        'planar' : Rotate around z axis until jvec is in the yz plane. Return xy.
        'edgeon' : jvec perpendicular to z-axis.
        'faceon' : jvec parallel to z-axis.

    Returns
    -------
    x_new : array-like
        2-dimensional array (Nx2) with the projected particle positions.

    """
    assert len(jvec) == x.shape[1] == 3

    # Define unit basis vectors of new reference frame
    if proj_kind == 'yz':
        e1 = np.array([0,1,0])
        e2 = np.array([0,0,1])
    elif proj_kind == 'zx':
        e1 = np.array([0,0,1])
        e2 = np.array([1,0,0])
    elif proj_kind == 'xy':
        e1 = np.array([1,0,0])
        e2 = np.array([0,1,0])
    elif proj_kind == 'planar':
        # New y-axis is the projection of jvec onto the xy plane
        e2 = np.array([jvec[0], jvec[1], 0.0])
        e2 = e2 / np.linalg.norm(e2)  # normalize
        e1 = np.cross(e2, np.array([0,0,1]))
    elif proj_kind == 'edgeon':
        # New y-axis is aligned with jvec
        e2 = jvec[:] / np.linalg.norm(jvec)  # normalize
        e1 = np.cross(e2, np.array([0,0,1]))
    elif proj_kind == 'faceon':
        # New z-axis is aligned with jvec
        e3 = jvec[:] / np.linalg.norm(jvec)  # normalize
        # New x-axis is chosen to coincide with edge-on projection
        e1 = np.cross(e3, np.array([0,0,1]))
        e2 = np.cross(e3, e1)
    else:
        raise Exception('Projection kind not understood.')

    # Project onto new axes
    x_new = np.zeros((x.shape[0], 2), dtype=np.float64)
    x_new[:,0] = np.dot(x, e1)
    x_new[:,1] = np.dot(x, e2)

    return x_new

def get_hsml(x, y, z, num_neighbors):
    """
    Get distance to the Nth (usually 16th) nearest neighbor in 3D.

    Parameters
    ----------
    x : array-like
        x-coordinates of the particles.
    y : array-like
        y-coordinates of the particles.
    z : array-like
        z-coordinates of the particles.
    num_neighbors : int
        Specifies how many neighbors to search for.

    Returns
    -------
    hsml : array-like
        Distances to the Nth nearest neighbors.

    """
    data = np.empty((len(x), 3))
    data[:,0] = x.ravel()
    data[:,1] = y.ravel()
    data[:,2] = z.ravel()

    tree = cKDTree(data)
    res = tree.query(data, k=num_neighbors+1)
    hsml = res[0][:,-1]

    return hsml

def adaptive_smoothing(x, y, hsml, xcenters, ycenters, num_rhalfs, weights=None):
    """
    Do adaptive smoothing similar to Torrey et al. (2015).

    Parameters
    ----------
    x : array-like
        x-coordinates of the particles.
    y : array-like
        y-coordinates of the particles.
    hsml : array-like
        Smoothing lengths (same units as x and y).
    xcenters : array-like
        1-d array with the pixel centers along the x-axis
    ycenters : array-like
        1-d array with the pixel centers along the y-axis
    weights : array-like, optional
        Array of the same size as ``x`` and ``y`` with the particle
        weights, e.g., particle masses or fluxes. If ``None``, the
        particle number density is calculated.

    Returns
    -------
    H : array-like
        A 2D array with the density at each pixel center.

    """
    assert x.shape == y.shape
    if weights is None:
        weights = np.ones_like(x)

    # Make everything double
    x = np.float64(x)
    y = np.float64(y)
    hsml = np.float64(hsml)
    weights = np.float64(weights)

    # Ignore out-of-range particles
    locs_withinrange = (np.abs(x) < num_rhalfs) | (np.abs(y) < num_rhalfs)
    x = x[locs_withinrange]
    y = y[locs_withinrange]
    hsml = hsml[locs_withinrange]
    weights = weights[locs_withinrange]

    # Compile as:
    # gcc -o adaptive_smoothing.so -shared -fPIC adaptive_smoothing.c
    sphlib = np.ctypeslib.load_library('adaptive_smoothing', codedir)
    sphlib.add.restype = None
    sphlib.add.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
    ]

    start = time.time()
    print('Doing adaptive smoothing...')

    X, Y = np.meshgrid(xcenters, ycenters)
    ny, nx = X.shape
    Y_flat, X_flat = Y.ravel(), X.ravel()
    Z_flat = np.zeros_like(X_flat)
    sphlib.add(
        X_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        Y_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        Z_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(X.shape[1]),
        ctypes.c_int(X.shape[0]),
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        weights.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        hsml.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(x.size),
        ctypes.c_double(num_rhalfs))

    H = Z_flat.reshape(X.shape)

    print('Time: %g s.' % (time.time() - start))

    return H

def get_subfind_ids(snapnum, log_mstar_bin_lower, log_mstar_bin_upper, mstar):
    nsubs = len(mstar)

    # Mass bins
    mstar_bin_lower = 10.0**log_mstar_bin_lower / 1e10 * h
    mstar_bin_upper = 10.0**log_mstar_bin_upper / 1e10 * h
    num_mstar_bins = len(log_mstar_bin_lower)

    # Iterate over mass bins
    subfind_ids = []
    for mstar_bin_index in range(num_mstar_bins):
        mstar_min = 10.0**log_mstar_bin_lower[mstar_bin_index] / 1e10 * h
        mstar_max = 10.0**log_mstar_bin_upper[mstar_bin_index] / 1e10 * h

        # Only proceed if there are enough galaxies
        locs_valid = ((mstar >= mstar_min) * (mstar < mstar_max))
        if np.sum(locs_valid) == 0:
            print('Not enough galaxies. Skipping...')
            return

        # Iterate over subhalos
        for subfind_id in range(nsubs):
            # Only proceed if current subhalo is within mass range
            if locs_valid[subfind_id]:
                subfind_ids.append(subfind_id)

    return np.array(subfind_ids, dtype=np.int32)

def get_num_rhalfs_npixels(subfind_id):
    """
    Helper function to get the current values of num_rhalfs and npixels,
    considering that one and only one of the two is defined.
    """
    # Already checked this, but check again just in case:
    if num_rhalfs > 0 and npixels > 0:
        raise Exception('Only one of num_rhalfs and npixels should be defined ' +
                        '(the other should be -1).')
    if num_rhalfs > 0:
        cur_num_rhalfs = num_rhalfs
        cur_npixels = int(np.ceil(2.0*num_rhalfs*rhalf[subfind_id]/kpc_h_per_pixel))
    elif npixels > 0:
        assert rhalf[subfind_id] > 0
        cur_num_rhalfs = npixels*kpc_h_per_pixel/(2.0*rhalf[subfind_id])
        cur_npixels = npixels

    return cur_num_rhalfs, cur_npixels

def create_image_single_sub(subfind_id, pos, hsml_ckpc_h, fluxes):
    """Create image for a single subhalo."""
    cur_num_rhalfs, cur_npixels = get_num_rhalfs_npixels(subfind_id)
    print('cur_num_rhalfs = %.1f' % (cur_num_rhalfs))
    print('cur_npixels = %d' % (cur_npixels))

    # Periodic boundary conditions
    dx = pos[:] - sub_pos[subfind_id]
    dx = dx - (np.abs(dx) > 0.5*box_size) * np.copysign(box_size, dx - 0.5*box_size)

    # Normalize by rhalf
    dx = dx / rhalf[subfind_id]
    hsml = hsml_ckpc_h / rhalf[subfind_id]

    # Transform particle positions according to 'proj_kind' (2D projection)
    dx_new = transform(dx, jstar_direction[subfind_id], proj_kind=proj_kind)

    # Define 2D bins (in units of rhalf)
    xedges = np.linspace(-cur_num_rhalfs, cur_num_rhalfs, num=cur_npixels+1)
    yedges = np.linspace(-cur_num_rhalfs, cur_num_rhalfs, num=cur_npixels+1)
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])

    # Store images here
    image = np.zeros((num_filters,cur_npixels,cur_npixels), dtype=np.float32)

    # Iterate over broadband filters
    for i, filter_name in enumerate(filter_names):
        H = adaptive_smoothing(
            dx_new[:,0], dx_new[:,1], hsml, xcenters, ycenters, cur_num_rhalfs,
            weights=fluxes[i,:])
        # Store in array
        image[i,:,:] = H

    # Convert area units from rhalf^{-2} to pixel_size^{-2}
    pixel_size_rhalfs =  2.0 * cur_num_rhalfs / float(cur_npixels)  # in rhalfs
    image *= pixel_size_rhalfs**2

    # Create some header attributes
    header = fits.Header()
    header["BUNIT"] = ("counts/s/pixel", "Unit of the array values")
    header["CDELT1"] = (kpc_h_per_pixel * 1000.0 / h / (1.0 + z), "Coordinate increment along X-axis")
    header["CTYPE1"] = ("pc", "Physical units of the X-axis increment")
    header["CDELT2"] = (kpc_h_per_pixel * 1000.0 / h / (1.0 + z), "Coordinate increment along Y-axis")
    header["CTYPE2"] = ("pc", "Physical units of the Y-axis increment")
    header["PIXSCALE"] = (arcsec_per_pixel, "Pixel size in arcsec")
    header["USE_Z"] = (use_z, "Observed redshift of the source")
    for k in range(num_filters):
        header["FILTER%d" % (k)] = (filter_names[k], "Broadband filter index = %d" % (k))

    # Write to FITS file
    hdu = fits.PrimaryHDU(data=image, header=header)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto('%s/broadband_%d.fits' % (datadir, subfind_id))
    hdulist.close()

    print('Finished for subhalo %d.\n' % (subfind_id))

def create_images(object_id):
    """
    Create (adaptively smoothed) images for the chosen filters.
    Consider all relevant subhalos of a given FoF group.

    Parameters
    ----------
    object_id : int
        The ID of the object of interest. If ``use_fof`` is True,
        this corresponds to the FoF group ID. Otherwise, this is
        the subhalo ID.
    """
    if use_fof:
        # Get Subfind IDs that belong to the current FoF group.
        fof_subfind_ids = subfind_ids[fof_ids == object_id]
    else:
        # A bit of a hack. This way we only process the subhalo of
        # interest without rewriting too much code.
        fof_subfind_ids = np.array([object_id], dtype=np.int32)

    # Load stellar particle info
    start = time.time()
    print('Loading info from snapshot...')
    if use_fof:
        cat = il.snapshot.loadHalo(
            basedir, snapnum, object_id, parttype_stars,
            fields=['Coordinates', 'GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime'])
    else:
        cat = il.snapshot.loadSubhalo(
            basedir, snapnum, object_id, parttype_stars,
            fields=['Coordinates', 'GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime'])
    all_pos = cat['Coordinates']  # comoving kpc/h
    all_initial_masses = cat['GFM_InitialMass']
    all_metallicities = cat['GFM_Metallicity']
    all_formtimes = cat['GFM_StellarFormationTime']  # actually the scale factor
    print('Time: %f s.' % (time.time() - start))

    # Remove wind particles
    locs_notwind = all_formtimes >= 0
    pos = all_pos[locs_notwind]
    initial_masses = all_initial_masses[locs_notwind]
    metallicities = all_metallicities[locs_notwind]
    formtimes = all_formtimes[locs_notwind]

    # Prepare input for BC03 model
    initial_masses_Msun = initial_masses * 1e10 / h
    z_form = 1.0/formtimes - 1.0
    params = cosmo.CosmologicalParameters(suite=suite)
    stellar_ages_yr = (cosmo.t_Gyr(z, params) - cosmo.t_Gyr(z_form, params)) * 1e9

    # Get smoothing lengths in 3D (before making 2D projection)
    # once and for all, in simulation units [ckpc/h].
    # We temporarily center on any subhalo of the FoF group to account
    # for periodic boundary conditions.
    dx = pos[:] - sub_pos[fof_subfind_ids[0]]
    dx = dx - (np.abs(dx) > 0.5*box_size) * np.copysign(box_size, dx - 0.5*box_size)
    start = time.time()
    print('Doing spatial search...')
    hsml_ckpc_h = get_hsml(dx[:,0], dx[:,1], dx[:,2], num_neighbors)  # in rhalfs
    print('Time: %g s.' % (time.time() - start))

    # Get all fluxes for once and for all
    fluxes = np.empty((len(filter_names), len(initial_masses_Msun)), dtype=np.float64)
    for i, filter_name in enumerate(filter_names):
        fluxes[i,:] = get_fluxes(
            initial_masses_Msun, metallicities, stellar_ages_yr, filter_name)

    for subfind_id in fof_subfind_ids:
        create_image_single_sub(subfind_id, pos, hsml_ckpc_h, fluxes)

    print('Finished for object %d.\n' % (object_id))

def slave(jobs):
    """Slave function."""
    try:
        while True:
            object_id = jobs.get_nowait()
            create_images(object_id)
    except:
        pass  # an exception is raised when job queue is empty


if __name__ == '__main__':
    try:
        suite = sys.argv[1]
        basedir = sys.argv[2]
        amdir = sys.argv[3]
        filename_filters = sys.argv[4]
        stellar_photometrics_dir = sys.argv[5]
        writedir = sys.argv[6]
        codedir = sys.argv[7]
        snapnum = int(sys.argv[8])
        use_z = float(sys.argv[9])  # if -1, use intrinsic snapshot redshift
        proj_kind = sys.argv[10]  # 'yz', 'zx', 'xy', 'planar', 'faceon', 'edgeon'
        num_neighbors = int(sys.argv[11])  # for adaptive smoothing, usually 32
        num_rhalfs = float(sys.argv[12])  # on each side from the center, usually 7.5
        npixels = int(sys.argv[13])  # total # of pixels on each side, usually -1
        log_mstar_min = float(sys.argv[14])  # minimum log10(M*) of galaxies
        use_fof = bool(int(sys.argv[15]))  # If True, load particles from FoF group
        use_cf00 = bool(int(sys.argv[16]))  # If True, apply Charlot & Fall (2000)
        mock_type = sys.argv[17]  # 'pogs', 'sdss', etc.
        nprocesses = int(sys.argv[18])
    except:
        print('Arguments: suite basedir amdir filename_filters ' + 
              'stellar_photometrics_dir writedir codedir snapnum use_z ' +
              'proj_kind num_neighbors num_rhalfs npixels log_mstar_min ' +
              'use_cf00 use_fof mock_type nprocesses')
        sys.exit()

    # Check input
    if num_rhalfs > 0 and npixels > 0:
        raise Exception('Only one of num_rhalfs and npixels should be defined ' +
                        '(the other should be -1).')

    # Save images here
    synthdir = '%s/snapnum_%03d/galaxev/%s' % (writedir, snapnum, proj_kind)
    if use_fof:
        synthdir += '_fof'
    if use_cf00:
        synthdir += '_cf00'
    datadir = '%s/data' % (synthdir)
    if not os.path.lexists(datadir):
        os.makedirs(datadir)

    # Read filter names
    with open(filename_filters, 'r') as f:
        filter_names = list(map(lambda s: s.strip('\n'), f.readlines()))
    num_filters = len(filter_names)

    # Load some info from snapshot header
    with h5py.File('%s/snapdir_%03d/snap_%03d.0.hdf5' % (basedir, snapnum, snapnum), 'r') as f:
        header = dict(f['Header'].attrs.items())
        h = header['HubbleParam']
        z = header['Redshift']
        box_size = header['BoxSize']

    # Set redshift at which the galaxy is observed (unless specified)
    if use_z == -1:
        use_z = z

    # Get pixel scale (both angular and physical) and define "use_z",
    # which is the redshift at which the galaxy is assumed to be observed,
    # but is not necessarily equal to the intrinsic redshift "z" of the galaxy
    # (in particular, use_z = 0 for rest-frame photometry).
    if mock_type == 'generic':
        use_z = 0.0  # Rest-frame
        print('WARNING: setting use_z=0.0 (rest-frame photometry).')
        kpc_h_per_pixel = 0.25  # Fixed pixel scale in ckpc/h
        # Camera is 10 Mpc away from source in physical units:
        d_A_kpc_h = 10000.0 * h * (1.0 + z)  # ckpc/h
        rad_per_pixel = kpc_h_per_pixel / d_A_kpc_h
        arcsec_per_pixel = rad_per_pixel * (3600.0 * 180.0 / np.pi)
    elif mock_type == 'pogs':
        # If at the last snapshot, set ad hoc redshift
        if ((suite == 'Illustris' and snapnum == 135) or
            (suite == 'IllustrisTNG' and snapnum == 99)):
                use_z = 0.0485236299818  # corresponds to snapnum_last - 4
                print('WARNING: Setting use_z=%g.' % (use_z))
        arcsec_per_pixel = 0.25  # Chambers et al. (2016)
        rad_per_pixel = arcsec_per_pixel / (3600.0 * 180.0 / np.pi)
        # Note that the angular-diameter distance is expressed in comoving coordinates:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_A_kpc_h = cosmo.angular_diameter_distance_Mpc(use_z, params) * 1000.0 * h * (1.0+z)  # ckpc/h
        kpc_h_per_pixel = rad_per_pixel * d_A_kpc_h  # about 0.174 (ckpc/h)/pixel at z = 0.0485
    elif mock_type == 'sdss':
        # If at the last snapshot, set ad hoc redshift
        if ((suite == 'Illustris' and snapnum == 135) or
            (suite == 'IllustrisTNG' and snapnum == 99)):
                use_z = 0.0485236299818  # corresponds to snapnum_last - 4
                print('WARNING: Setting use_z=%g.' % (use_z))
        arcsec_per_pixel = 0.396  # https://www.sdss.org/instruments/camera/
        rad_per_pixel = arcsec_per_pixel / (3600.0 * 180.0 / np.pi)
        # Note that the angular-diameter distance is expressed in comoving coordinates:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_A_kpc_h = cosmo.angular_diameter_distance_Mpc(use_z, params) * 1000.0 * h * (1.0+z)  # ckpc/h
        kpc_h_per_pixel = rad_per_pixel * d_A_kpc_h  # about 0.174 (ckpc/h)/pixel at z = 0.0485
    elif mock_type == 'kids':
        # If at the last snapshot, set ad hoc redshift
        if ((suite == 'Illustris' and snapnum == 135) or
            (suite == 'IllustrisTNG' and snapnum == 99)):
                use_z = 0.15274876890238098  # corresponds to snapnum_last - 12
                print('WARNING: Setting use_z=%g.' % (use_z))
        arcsec_per_pixel = 0.21  # Lingyu Wang, private communication
        rad_per_pixel = arcsec_per_pixel / (3600.0 * 180.0 / np.pi)
        # Note that the angular-diameter distance is expressed in comoving coordinates:
        params = cosmo.CosmologicalParameters(suite=suite)
        d_A_kpc_h = cosmo.angular_diameter_distance_Mpc(use_z, params) * 1000.0 * h * (1.0+z)  # ckpc/h
        kpc_h_per_pixel = rad_per_pixel * d_A_kpc_h  # about 0.174 (ckpc/h)/pixel at z = 0.0485
    else:
        print('mock_type not understood.')
        sys.exit()

    print('use_z = %g' % (use_z))
    print('arcsec_per_pixel = %g' % (arcsec_per_pixel))
    print('kpc_h_per_pixel = %g' % (kpc_h_per_pixel))

    # Load subhalo info
    start = time.time()
    print('Loading subhalo info...')
    mstar = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
    rhalf = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloHalfmassRadType'])[:, parttype_stars]
    sub_pos = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloPos'])
    with h5py.File('%s/jstar_%03d.hdf5' % (amdir, snapnum), 'r') as f:
        jstar_direction = f['jstar_direction'][:]
    sub_gr_nr = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloGrNr'])
    nsubs = len(mstar)
    print('Time: %f s.' % (time.time() - start))

    # For performance checks
    start_all = time.time()

    # Define stellar mass bins
    log_mstar_bin_lower = np.array([log_mstar_min])
    log_mstar_bin_upper = np.array([13.0])  # hardcoded since we don't expect larger M*
    mstar_bin_lower = 10.0**log_mstar_bin_lower / 1e10 * h
    mstar_bin_upper = 10.0**log_mstar_bin_upper / 1e10 * h

    # Get list of relevant Subfind IDs
    subfind_ids = get_subfind_ids(snapnum, log_mstar_bin_lower, log_mstar_bin_upper, mstar)
    # Get associated FoF group IDs
    fof_ids = sub_gr_nr[subfind_ids]
    unique_fof_ids = np.unique(fof_ids)
    # Print Subfind IDs to a file
    filename = '%s/subfind_ids.txt' % (synthdir)
    with open(filename, 'w') as f:
        for sub_index in subfind_ids:
            f.write('%d\n' % (sub_index))

    # Create list of "generic" objects (halo or subhalo)
    if use_fof:
        object_ids = unique_fof_ids
    else:
        object_ids = subfind_ids

    # Create images
    if nprocesses == 1:
        for object_id in object_ids:
            create_images(object_id)
    else:
        # The FoF and subhalo versions require different parallelization
        # strategies (at least with the multiprocessing package). The Queue
        # results in better load balancing, but seems to involve too much
        # communication when dealing with small tasks (i.e. per subhalo).
        if use_fof:
            # See http://stevemorphet.weebly.com/python/python-multiprocessing
            pool = []       # instantiate pool of processes
            jobs = Queue()  # instantiate job queue
            # Instantiate N slave processes
            for proc_id in range(nprocesses):
                pool.append(Process(target=slave, args=(jobs,)))
            # Populate the job queue
            for object_id in object_ids:
                jobs.put(object_id)
            # Start the slaves
            for slave in pool:
                slave.start()
            # Wait for the slaves to finish processing
            for slave in pool:
                slave.join()
        else:
            with Pool(nprocesses) as p:
                p.map(create_images, list(object_ids))

    print('Total time: %f s.\n' % (time.time() - start_all))

