#!/bin/bash

# This batch script can be modified to retrieve a set of filter
# curves from the SVO Filter Profile Service via Astroquery.

# Any user-specified name to identify the current set of synthetic images:
MOCK_SET=hsc

# A folder to store all the output corresponding to this set of synthetic images
# (its subdirectories will include all the simulations, snapshots, etc.)
WRITEDIR=/path/to/imagedir/${MOCK_SET}

# Store filter curves in this directory (and create it):
FILTER_DIR=${WRITEDIR}/filter_curves
mkdir -p ${FILTER_DIR}

# Specify the filter IDs of interest here:
FILTER_IDS=(
"Subaru/HSC.g"
"Subaru/HSC.r"
"Subaru/HSC.i"
"Subaru/HSC.z"
"Subaru/HSC.Y"
)

# Filter IDs will be stored in filters.txt (should not exist yet):
FILENAME=${WRITEDIR}/filters.txt
if [ -f "${FILENAME}" ]; then
    echo "File ${FILENAME} already exists!"
    exit 1
fi
touch ${FILENAME}

# Retrieve filter curves and list them in filters.txt:
for FILTER_ID in ${FILTER_IDS[@]}; do
    echo ${FILTER_ID} >> ${FILENAME}
    # Create new directory for the facility (e.g. Subaru), if necessary:
    FACILITY=$(echo "${FILTER_ID}" | cut -d "/" -f 1)
    if [ ! -d "${FILTER_DIR}/${FACILITY}" ]; then
        mkdir ${FILTER_DIR}/${FACILITY}
    fi
    # Retrieve filter curve via Astroquery and store it in ASCII format:
    python -c "from astroquery.svo_fps import SvoFps; \
               data = SvoFps.get_transmission_data('${FILTER_ID}'); \
               assert data['Wavelength'].unit.name == 'AA'; \
               data.write('${FILTER_DIR}/${FILTER_ID}', format='ascii.no_header')"
    if [ $? -eq 0 ]; then echo "Retrieved filter ${FILTER_ID}"; fi
done
