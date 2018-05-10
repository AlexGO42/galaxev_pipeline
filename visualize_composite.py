import numpy as np
import os
import sys
if sys.version_info[0] == 2:  # Python 2
    if not sys.modules.has_key('matplotlib'):
        import matplotlib
        matplotlib.use('agg')
elif sys.version_info[0] == 3:  # Python 3
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from astropy.visualization import LogStretch, AsinhStretch, make_lupton_rgb
import h5py
import time
import ctypes
import illustris_python as il

def normalize(image):
    locs = np.isfinite(image)
    m, M = np.min(image[locs]), np.max(image[locs])
    return (image-m) / (M-m)

def scale_image(image, kind='arcsinh'):
    """
    Scale image for visualization.

    Parameters
    ----------
    image : array-like
        The image.
    kind : str, optional
        The type of scaling. Valid options are:
        
        'linear' : Linear scaling.
        'log' : Log scaling
        'arcsinh' : Arcsinh scaling
        'lupton' : Lupton (2004) scaling
    
    Returns
    -------
    retval : array-like
        Scaled image.

    """
    if kind == 'linear':
        return normalize(image)
    elif kind == 'log':
        log_stretch = LogStretch(a=1000.0)
        return log_stretch(normalize(image))
    elif kind == 'arcsinh':
        asinh_stretch = AsinhStretch(a=0.002)
        return asinh_stretch(normalize(image))
    elif kind == 'lupton':
        assert image.shape[2] == 3
        normalized_image = normalize(image)
        image_r = normalized_image[:,:,0]
        image_g = normalized_image[:,:,1]
        image_b = normalized_image[:,:,2]
        return make_lupton_rgb(image_r, image_g, image_b, stretch=0.02, Q=8)
    else:
        raise Exception('Not implemented: %s' % (kind))


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
        for sub_index in range(nsubs):
            # Only proceed if current subhalo is within mass range
            if locs_valid[sub_index]:
                subfind_ids.append(sub_index)

    return subfind_ids

def read_file(sub_index):
    """
    Read data from an HDF5 file, apply noise and save figure.
    """
    datafilename = '%s/composite_%d_rhalfs_ngb_%d_sub_%d_%s.hdf5' % (
        datadir, int(num_rhalfs), num_neighbors, sub_index, proj_kind)

    nbins = int(np.ceil(nbins_per_softening_length * (2*num_rhalfs) *
                          rhalf[sub_index] / softening_length_kpc_h))
    bin_size_kpc_h = 2.0 * num_rhalfs * rhalf[sub_index] / float(nbins)  # in kpc/h
    npixels = nbins**2

    f_read = h5py.File(datafilename, 'r')

    # Get number of particles in image
    counts_map = f_read['CountsMap'][:]  # counts/pixel
    ncounts = np.sum(counts_map)
    
    # Sanity checks
    ny, nx = counts_map.shape
    assert nbins == nx
    assert npixels == nx*ny

    fig = plt.figure(figsize=(4.5, 4.5))
    ax = fig.add_subplot(111)

    # Iterate over g,r,i bands
    composite = np.zeros((nbins,nbins,3), dtype=np.float64)
    for i in range(len(bands)):
        H = f_read['Composite'][i,:,:]

        # Add "consistent" noise -- corresponds to 1 count/pixel
        total_flux = np.sum(H)
        mean_flux = total_flux / float(npixels)
        sigma_noise = total_flux / float(ncounts)
        print('Average S/N:', mean_flux / sigma_noise)

        H += np.random.normal(loc=0.0, scale=sigma_noise, size=H.shape)

        composite[:,:,-(i+1)] = H

    # Scale image
    composite = scale_image(composite, kind='log')
    
    # Plot 2D histogram using imshow
    ax.imshow(composite, origin="lower", interpolation="gaussian")
    ax.set_facecolor('k')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # Some text
    ax.text(0.5, 0.97, 'Stellar Light Composite', fontsize=12, fontweight='bold', color='w',
            horizontalalignment='center', verticalalignment='top',
            transform=ax.transAxes)
    # Show scale
    ax.plot([0.75*nbins, 0.95*nbins], [0.05*nbins, 0.05*nbins], 'w', lw=2)
    ax.text(0.85*nbins, 0.1*nbins, '%.0f kpc' % (0.2*nbins*bin_size_kpc_h/h),
            fontsize=9, fontweight='bold', color='w',
            horizontalalignment='center', verticalalignment='center')
    ax.set_xlim(0,nbins)
    ax.set_ylim(0,nbins)

    figname = '%s/composite_mstar_%.1f_%d_rhalfs_ngb_%d_sub_%d_%s.png' % (
        imgdir, np.log10(mstar[sub_index]*1e10/h), int(num_rhalfs),
        num_neighbors, sub_index, proj_kind)
    fig.savefig(figname)
    plt.close(fig)

    f_read.close()

    print('Finished for subhalo %d.' % (sub_index))

if __name__ == '__main__':
    try:
        basedir = sys.argv[1]
        amdir = sys.argv[2]
        writedir = sys.argv[3]
        snapnum = int(sys.argv[4])
        num_rhalfs = float(sys.argv[5])  # how many rhalfs on each side from the center
        proj_kind = sys.argv[6]  # 'none', 'planar', 'faceon', 'edgeon'
        nprocesses = int(sys.argv[7])
    except:
        print('Arguments: basedir amdir writedir snapnum num_rhalfs proj_kind nprocesses')
        sys.exit()

    parttype_stars = 4
    nbins_per_softening_length = 2.0  # roughly corresponds to Pan-STARRS r-band
    num_neighbors = 16  # for adaptive smoothing

    # SDSS g,r,i bands correspond to entries 4,5,6
    bands = [4, 5, 6]
    num_bands = len(bands)

    datadir = '%s/snap_%03d/data' % (writedir, snapnum)
    imgdir = '%s/snap_%03d/images' % (writedir, snapnum)
    # Make sure that write directories exist
    if not os.path.exists(imgdir):
        os.popen('mkdir -p %s' % (imgdir))
    if not os.path.exists(datadir):
        os.popen('mkdir -p %s' % (datadir))

    # Load some info from snapshot header
    with h5py.File('%s/snapdir_%03d/snap_%03d.0.hdf5' % (basedir, snapnum, snapnum), 'r') as f:
        header = dict(f['Header'].attrs.items())
        h = header['HubbleParam']
        box_size = header['BoxSize']

    # Define softening length at z=0 (maybe should modify for z>0)
    if 'L75n1820' in basedir:
        softening_length_kpc_h = 0.5  # kpc/h
    elif ('L75n910' in basedir) or ('L205n2500' in basedir):
        softening_length_kpc_h = 1.0  # kpc/h
    elif ('L75n455' in basedir) or ('L205n1250' in basedir):
        softening_length_kpc_h = 2.0  # kpc/h
    elif 'L205n625' in basedir:
        softening_length_kpc_h = 4.0  # kpc/h
    else:
        raise Exception('Please define softening length.')

    # For performance checks
    start_all = time.time()

    # Define stellar mass bins
    log_mstar_bin_lower = np.array([9.0])
    log_mstar_bin_upper = np.array([13.0])
    mstar_bin_lower = 10.0**log_mstar_bin_lower / 1e10 * h
    mstar_bin_upper = 10.0**log_mstar_bin_upper / 1e10 * h
    num_mstar_bins = len(log_mstar_bin_lower)

    # Load subhalo info
    start = time.time()
    print('Loading subhalo info...')
    mstar = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
    rhalf = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloHalfmassRadType'])[:, parttype_stars]
    with h5py.File('%s/jstar_%03d.hdf5' % (amdir, snapnum), 'r') as f:
        jstar_direction = f['jstar_direction'][:]
    nsubs = len(mstar)
    print('Time: %f s.' % (time.time() - start))

    # ~ # Get list of relevant Subfind IDs
    # ~ subfind_ids = get_subfind_ids(snapnum, log_mstar_bin_lower, log_mstar_bin_upper, mstar)

    subfind_ids = [6]

    if nprocesses == 1:
        for sub_index in subfind_ids:
            read_file(sub_index)
    else:
        p = Pool(nprocesses)
        p.map(read_file, subfind_ids)

    print('Total time: %f s.\n' % (time.time() - start_all))
