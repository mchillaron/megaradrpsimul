#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Smooth the images using Savitzky-Golay and median filters."""

import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter

from megaradrpsimul import CCDregions

def smooth_frames(cleaned_images):
    """Smooth the images using Savitzky-Golay and median filters.
    
    Parameters
    ----------
    cleaned_images : dict
        Dictionary with cleaned images to be smoothed.
    Returns
    -------
    smoothed_images : dict
        Dictionary with smoothed images."""
    

    regions_kernel = CCDregions.regions_kernel
    smoothed_images = {}

    for idx, (name, image) in enumerate(cleaned_images.items()):
        data = image
        naxis2, naxis1 = data.shape

        # create an image of the same dimensions full of zeros
        zero_raw_image = np.zeros((naxis2, naxis1))

        for i, (region_name, region_params) in enumerate(regions_kernel.items()): 
            slice2d = region_params['slice2d']
            num_filters_SG = region_params['num_filters_SG']
            
            region_data = data[slice2d.python]
            
            if num_filters_SG == 0:
                if 'median_size' in region_params:
                    median_size = region_params['median_size']
                    filtered_data = median_filter(region_data, size=median_size)  # vertical median filter in these regions
                    print(f"Median filter applied to {region_name}")
                else:
                    raise ValueError(f"'median_size' is missing in region {region_name}")
                
            elif num_filters_SG != 0:
                size1_SG = region_params['size1_SG']
                pol_order1_SG = region_params['pol_order1_SG']
                axis1_SG = region_params['axis1_SG']

                if size1_SG % 2 == 0:
                    size1_SG += 1    # make sure it is odd
                
                filtered_data = savgol_filter(region_data, window_length=size1_SG, polyorder=pol_order1_SG, axis=axis1_SG)   # Apply only one Savitzky-Golay filter
                print(f"Savitzky-Golay applied in {region_name}")
        
                if num_filters_SG == 2:
                    size2_SG = region_params['size2_SG']
                    pol_order2_SG = region_params['pol_order2_SG']
                    axis2_SG = region_params['axis2_SG']

                    if size2_SG % 2 == 0:
                        size2_SG += 1 

                    filtered_data = savgol_filter(filtered_data, window_length=size2_SG, 
                                                  polyorder=pol_order2_SG, axis=axis2_SG)
                    print(f"Second Savitzky-Golay filter applied in {region_name}")

            zero_raw_image[slice2d.python] = filtered_data
            
        # Save the image in the dictionary "smoothed_{idx}"
        smoothed_image_name = f"smoothed_{idx}"
        smoothed_images[smoothed_image_name] = zero_raw_image

        if np.any(zero_raw_image == 0):
            print("There are still zeros in the smoothed bias image.")
        else:
            print("No zeros in the smoothed bias image.")
            
    return smoothed_images