"""
The idea of these "generic" images is to sacrifice some observational
realism in order to make generic images of any galaxy property. For now,
we make "generic" images of rest-frame SDSS g,r,i,z bands (pre-calculated
in the Subfind catalogs), stellar mass, and SFR.
Based on "create_images.py" and "generic_image_test.ipynb".
"""

import numpy as np
import h5py
import os
import sys
from scipy.spatial import cKDTree
from astropy.io import fits
from mpi4py import MPI
import time

import illustris_python as il
import cosmology as cosmo

from create_images import transform, get_hsml, adaptive_smoothing, get_subfind_ids

KILL_TAG = 1
WORK_TAG = 0
parttype_gas = 0
parttype_stars = 4

def get_softening_length_ckpc_h(z, simulation):
    """
    Return the gravitational softening length in *comoving* kpc/h
    for star particles (only implemented for IllustrisTNG).
    """
    if simulation == 'L75n1820TNG':
        max_softening_length_kpc_h = 0.5  # *physical* kpc/h
    elif simulation == 'L75n910TNG':
        max_softening_length_kpc_h = 1.0  # *physical* kpc/h
    elif simulation == 'L75n455TNG':
        max_softening_length_kpc_h = 2.0  # *physical* kpc/h
    elif simulation == 'L35n2160TNG':
        max_softening_length_kpc_h = 0.195  # *physical* kpc/h
    elif simulation == 'L35n1080TNG':
        max_softening_length_kpc_h = 0.39  # *physical* kpc/h
    elif simulation == 'L35n540TNG':
        max_softening_length_kpc_h = 0.78  # *physical* kpc/h
    elif simulation == 'L35n270TNG':
        max_softening_length_kpc_h = 1.56  # *physical* kpc/h
    elif simulation == 'L205n2500TNG':
        max_softening_length_kpc_h = 1.0  # *physical* kpc/h
    elif simulation == 'L205n1250TNG':
        max_softening_length_kpc_h = 2.0  # *physical* kpc/h
    elif simulation == 'L205n625TNG':
        max_softening_length_kpc_h = 4.0  # *physical* kpc/h
    else:
        raise NotImplementedError(simulation)

    if z < 1:
        softening_length_ckpc_h = (1.0 + z) * max_softening_length_kpc_h  # comoving kpc/h
    else:
        softening_length_ckpc_h = 2.0 * max_softening_length_kpc_h  # comoving kpc/h

    return softening_length_ckpc_h

def create_image_single_sub(subfind_id, pos, hsml_ckpc_h, fluxes, xcenters, ycenters):
    """
    Create generic image for a single subhalo. The "fluxes" can also
    refer to masses, sfr, etc.
    """
    # Periodic boundary conditions (center at center of current subhalo)
    dx = pos[:] - sub_pos[subfind_id]
    dx = dx - (np.abs(dx) > 0.5*box_size) * np.copysign(box_size, dx - 0.5*box_size)

    # Normalize by rhalf
    dx = dx / rhalf[subfind_id]
    hsml = hsml_ckpc_h / rhalf[subfind_id]

    # Transform particle positions according to 'proj_kind' (2D projection)
    dx_new = transform(dx, jstar_direction[subfind_id], proj_kind=proj_kind)

    # Create image "layer"
    H = adaptive_smoothing(
        dx_new[:,0], dx_new[:,1], hsml, xcenters, ycenters, num_rhalfs,
        codedir, weights=fluxes)

    return H

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

    # --------------------------- STARS --------------------------------

    # Load stellar particle info
    start = time.time()
    print('Loading stellar info from snapshot...')
    if use_fof:
        cat_stars = il.snapshot.loadHalo(
            basedir, snapnum, object_id, parttype_stars,
            fields=['Coordinates', 'Masses', 'GFM_StellarFormationTime', 'GFM_StellarPhotometrics'])
    else:
        cat_stars = il.snapshot.loadSubhalo(
            basedir, snapnum, object_id, parttype_stars,
            fields=['Coordinates', 'Masses', 'GFM_StellarFormationTime', 'GFM_StellarPhotometrics'])
    all_pos = cat_stars['Coordinates']  # comoving kpc/h
    all_masses = cat_stars['Masses']
    all_formtimes = cat_stars['GFM_StellarFormationTime']  # actually the scale factor
    all_magnitudes = cat_stars['GFM_StellarPhotometrics']  # U, B, V, K, g, r, i, z -- use the last 4
    print('Time: %f s.' % (time.time() - start))

    # Remove wind particles
    locs_notwind = all_formtimes >= 0
    pos = all_pos[locs_notwind]
    masses = all_masses[locs_notwind]
    formtimes = all_formtimes[locs_notwind]
    magnitudes = all_magnitudes[locs_notwind]
    fluxes = 10.0**(-2.0/5.0*magnitudes)

    # Get smoothing lengths in 3D (before making 2D projection)
    # once and for all, in simulation units [ckpc/h].
    # We temporarily center on any subhalo of the FoF group to account
    # for periodic boundary conditions.
    dx = pos[:] - pos[0]
    dx = dx - (np.abs(dx) > 0.5*box_size) * np.copysign(box_size, dx - 0.5*box_size)
    start = time.time()
    print('Doing spatial search...')
    hsml_ckpc_h = get_hsml(dx, num_neighbors)
    print('Time: %g s.' % (time.time() - start))

    # ---------------------------- GAS ---------------------------------

    # Load gas cell info
    start = time.time()
    print('Loading gas info from snapshot...')
    if use_fof:
        cat_gas = il.snapshot.loadHalo(
            basedir, snapnum, object_id, parttype_gas,
            fields=['Coordinates', 'Masses', 'StarFormationRate'])
    else:
        cat_gas = il.snapshot.loadSubhalo(
            basedir, snapnum, object_id, parttype_gas,
            fields=['Coordinates', 'StarFormationRate'])
    print('Time: %f s.' % (time.time() - start))

    # Only bother for galaxies with at least 2 gas cells
    has_gas = cat_gas['count'] >= 2
    if has_gas:
        pos_gas = cat_gas['Coordinates']  # comoving kpc/h
        sfr = cat_gas['StarFormationRate']

        # Get smoothing lengths in 3D (before making 2D projection)
        # once and for all, in simulation units [ckpc/h].
        # We temporarily center on any subhalo of the FoF group to account
        # for periodic boundary conditions.
        dx_gas = pos_gas[:] - pos_gas[0]
        dx_gas = dx_gas - (np.abs(dx_gas) > 0.5*box_size) * np.copysign(box_size, dx_gas - 0.5*box_size)
        start = time.time()
        print('Doing spatial search...')
        hsml_ckpc_h_gas = get_hsml(dx_gas, num_neighbors)
        print('Time: %g s.' % (time.time() - start))

    # ------------------------- MAKE IMAGES -----------------------------

    # Iterate over subhalos in current FoF
    for subfind_id in fof_subfind_ids:
        # Define 2D bins (in units of rhalf)
        cur_npixels = int(np.ceil(2.0*num_rhalfs*rhalf[subfind_id]/ckpc_h_per_pixel))
        xedges = np.linspace(-num_rhalfs, num_rhalfs, num=cur_npixels+1)
        yedges = np.linspace(-num_rhalfs, num_rhalfs, num=cur_npixels+1)
        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])

        # Store images here
        image = np.zeros((num_filters + 2, cur_npixels, cur_npixels), dtype=np.float32)

        # Iterate over broadband filters
        for i in range(num_filters):
            image[i,:,:] = create_image_single_sub(
                subfind_id, pos, hsml_ckpc_h, fluxes[:, filter_indices[i]], xcenters, ycenters)
        # Same for stellar mass
        image[num_filters,:,:] = create_image_single_sub(
            subfind_id, pos, hsml_ckpc_h, masses, xcenters, ycenters)
        # Same for SFR; only bother if there are at least 2 gas cells
        if has_gas:  # otherwise image = 0
            image[num_filters+1,:,:] = create_image_single_sub(
                subfind_id, pos_gas, hsml_ckpc_h_gas, sfr, xcenters, ycenters)

        # Convert area units from rhalf^{-2} to pixel_size^{-2}
        pixel_size_rhalfs =  2.0 * num_rhalfs / float(cur_npixels)  # in rhalfs
        image *= pixel_size_rhalfs**2

        # Create some header attributes
        header = fits.Header()
        header["BUNIT"] = ("counts/s/pixel", "Unit of the array values")
        header["CDELT1"] = (ckpc_h_per_pixel * 1000.0 / h / (1.0 + z), "Coordinate increment along X-axis")
        header["CTYPE1"] = ("pc", "Physical units of the X-axis increment")
        header["CDELT2"] = (ckpc_h_per_pixel * 1000.0 / h / (1.0 + z), "Coordinate increment along Y-axis")
        header["CTYPE2"] = ("pc", "Physical units of the Y-axis increment")
        header["PIXSCALE"] = (arcsec_per_pixel, "Pixel size in arcsec")
        header["USE_Z"] = (use_z, "Observed redshift of the source")
        for k in range(num_filters):
            header["LAYER%d" % (k)] = (filter_names[k], "Broadband filter %s" % (filter_names[k]))
        header["LAYER%d" % (num_filters)] = ('mstar', "Stellar mass distribution")
        header["LAYER%d" % (num_filters+1)] = ('sfr', "SFR distribution")

        # Write to FITS file
        hdu = fits.PrimaryHDU(data=image, header=header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto('%s/broadband_%d.fits' % (datadir, subfind_id))
        hdulist.close()

        print('Finished for subhalo %d.\n' % (subfind_id))

    print('Finished for object %d.\n' % (object_id))

def master(comm):
    """Master process (to be run by process with rank 0)."""
    size = comm.Get_size()
    status = MPI.Status()
    
    dummy_arr = np.zeros(1, dtype=np.uint32)

    # Initialize by sending one unit of work to each slave
    cur_pos = 0
    for k in range(1, size):
        object_id = object_ids[cur_pos]
        comm.Send(np.array([object_id], dtype=np.uint32), dest=k, tag=WORK_TAG)
        cur_pos += 1

    # While there is more work...
    while cur_pos < len(object_ids):
        object_id = object_ids[cur_pos]
        # Get results from slave
        comm.Recv(dummy_arr, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        # Send another unit of work to slave
        comm.Send(np.array([object_id], dtype=np.uint32), dest=status.source, tag=WORK_TAG)
        # Next iteration
        cur_pos += 1
        print('Slave %i did object %i.' % (status.source, status.tag))

    # Get remaining results and kill slaves
    for k in range(1, size):
        comm.Recv(dummy_arr, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        comm.Send(dummy_arr, dest=status.source, tag=KILL_TAG)

def slave(comm):
    """ Slave process (to process one unit of work)."""
    status = MPI.Status()
    object_id = np.zeros(1, dtype=np.uint32)

    # Iterate until slaves receives the KILL_TAG
    while True:
        comm.Recv(object_id, source=0, tag=MPI.ANY_TAG, status=status)

        if status.tag == KILL_TAG:
            return

        # Do the work
        create_images(object_id[0])

        # Let the master know that the job is done
        comm.Send(np.zeros(1, dtype=np.uint32), dest=0, tag=object_id[0])


if __name__ == '__main__':
    try:
        suite = sys.argv[1]
        simulation = sys.argv[2]
        basedir = sys.argv[3]
        amdir = sys.argv[4]
        writedir = sys.argv[5]
        codedir = sys.argv[6]
        snapnum = int(sys.argv[7])
        proj_kind = sys.argv[8]  # 'yz', 'zx', 'xy', 'planar', 'faceon', 'edgeon'
        num_neighbors = int(sys.argv[9])  # for adaptive smoothing, usually 32
        num_rhalfs = float(sys.argv[10])  # on each side from the center, usually 7.5
        log_mstar_min = float(sys.argv[11])  # minimum log10(M*) of galaxies
        use_fof = bool(int(sys.argv[12]))  # If True, load particles from FoF group
        nprocesses = int(sys.argv[13])
    except:
        print('Arguments: suite basedir amdir writedir codedir snapnum ' +
              'proj_kind num_neighbors num_rhalfs log_mstar_min ' +
              'use_fof nprocesses')
        sys.exit()

    # MPI stuff
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # Save images here
    synthdir = '%s/snapnum_%03d/galaxev/%s' % (writedir, snapnum, proj_kind)
    if use_fof:
        synthdir += '_fof'
    datadir = '%s/data' % (synthdir)

    # Create write directory if it does not exist
    if rank == 0:
        if not os.path.lexists(datadir):
            os.makedirs(datadir)
    comm.Barrier()

    # Define filter names
    filter_names = ['sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']
    filter_indices = [4, 5, 6, 7]
    num_filters = len(filter_names)

    # Load some info from snapshot header
    with h5py.File('%s/snapdir_%03d/snap_%03d.0.hdf5' % (basedir, snapnum, snapnum), 'r') as f:
        header = dict(f['Header'].attrs.items())
        h = header['HubbleParam']
        z = header['Redshift']
        box_size = header['BoxSize']

    # Get pixel scale (both angular and physical) and define "use_z",
    # which is the redshift at which the galaxy is assumed to be observed,
    # but is not necessarily equal to the intrinsic redshift "z" of the galaxy.
    use_z = 0.0  # Rest-frame
    softening_length_ckpc_h = get_softening_length_ckpc_h(z, simulation)
    # Fixed pixel scale in ckpc/h -- should be 0.25 ckpc/h for TNG100 at z=0:
    ckpc_h_per_pixel = softening_length_ckpc_h / 2.0
    # Camera is 10 Mpc away from source in physical units:
    d_A_ckpc_h = 10000.0 * h * (1.0 + z)  # ckpc/h
    rad_per_pixel = ckpc_h_per_pixel / d_A_ckpc_h
    arcsec_per_pixel = rad_per_pixel * (3600.0 * 180.0 / np.pi)


    # ------------ MPI PARALLELIZATION STARTS HERE ------------
    
    if rank == 0:
        print('use_z = %g' % (use_z))
        print('arcsec_per_pixel = %g' % (arcsec_per_pixel))
        print('ckpc_h_per_pixel = %g' % (ckpc_h_per_pixel))

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
    else:
        rhalf = None
        jstar_direction = None

    # For simplicity, all processes will have a copy of these arrays:
    comm.Barrier()
    rhalf = comm.bcast(rhalf, root=0)
    sub_pos = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloPos'])
    jstar_direction = comm.bcast(jstar_direction, root=0)

    if rank == 0:
        # For performance checks
        start_all = time.time()

        # Define stellar mass bins
        log_mstar_bin_lower = np.array([log_mstar_min])
        log_mstar_bin_upper = np.array([13.0])  # hardcoded since we don't expect larger M*
        mstar_bin_lower = 10.0**log_mstar_bin_lower / 1e10 * h
        mstar_bin_upper = 10.0**log_mstar_bin_upper / 1e10 * h

        # Get list of relevant Subfind IDs
        subfind_ids = get_subfind_ids(snapnum, log_mstar_bin_lower, log_mstar_bin_upper, mstar, h)
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
            start_time = MPI.Wtime()
            master(comm)
            end_time = MPI.Wtime()
            print("MPI Wtime: %f s.\n" % (end_time - start_time))

        print('Total time: %f s.\n' % (time.time() - start_all))

    else:
        slave(comm)

