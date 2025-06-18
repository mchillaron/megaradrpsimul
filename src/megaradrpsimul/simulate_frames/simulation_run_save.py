#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Simulate images using the generator and save images."""

from astropy.io import fits
import shutil
import numpy as np


def simulation_run_save(data_work_dir, img_path, generator, type_image): 
    """Simulate images using the generator and save images.

    Parameters
    ----------
    data_work_dir : Path instance
        Path to the data work directory.
    img_path : Path instance
        Path to the image to be simulated.
    type_image : str
        Type of image ('bias', 'object').
    generator : instance of SimulateCCDExposure
        Generator object for simulating images.
    """
    destination_dir = data_work_dir
    image_exposure = generator.run(imgtype=type_image)

    new_filename = img_path.name  
    destination_file = destination_dir / new_filename  # save the simulated image with the same name as the original

    # Copying the original image to the destination directory
    # then we only update the data with the simulated data
    shutil.copy(img_path, destination_file)
    print(f"Image {img_path} copied at: {destination_file}")
    
    with fits.open(destination_file, mode='update') as hdul:
        data_round = np.round(image_exposure.data).astype(np.uint16)  
        hdul[0].data = data_round
        hdul.flush()
        print(f"Updated data for: {destination_file}")