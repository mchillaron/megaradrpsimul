#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Clean the images from cosmic rays using the regions defined in CCDregions.py"""

from astropy.io import fits
from tqdm import tqdm
import teareduce as tea

from megaradrpsimul import CCDregions

def cosmicray_cleaning(list_bias):
    """Clean the images from cosmic rays using the regions defined in CCDregions.py.
    
    This function uses the tea.cr2images function to clean the images, and the
    regions are defined in the CCDregions.py file. The function returns a dictionary 
    with the cleaned images, whose keys are the names of the images.

    Parameters
    ----------
    list_bias : list of Path instances
        List of image paths to be cleaned.
        
    Returns
    -------
    cleaned_images : dict
        Dictionary with cleaned images."""

    print('cleaning cosmic rays')
    regions_cosmicrays = CCDregions.regions_cosmicrays
    cleaned_images = {}

    for idx, img_path in tqdm(enumerate(list_bias), total=len(list_bias), desc="Image processing", unit="img"):
        with fits.open(img_path) as hdul:
            data = hdul[0].data.astype(float)  # we convert to float to avoid errors
            naxis2, naxis1 = data.shape        # dimensions of the image

            print(f"{idx}.Processing image {img_path} with shape ({naxis1}, {naxis2})")

            cleaned_data = data.copy()         # creating a copy of the image to do changes

            for region_name, region_params in regions_cosmicrays.items():
                slice2d = region_params['slice2d']
                print(f"Processing region {slice2d}")
                median_size = region_params['median_size']
                tsigma_peak = region_params['tsigma_peak']
                tsigma_tail = region_params['tsigma_tail']

                # cleaning the region defined by slice2d using tea.cr2images
                cleaned_region = tea.cr2images(
                    data1=data, 
                    median_size=median_size,
                    tsigma_peak=tsigma_peak,
                    tsigma_tail=tsigma_tail,
                    image_region = slice2d,
                    return_masks = False,
                    debug_level = 0
                )
                cleaned_data[slice2d.python] = cleaned_region[0][slice2d.python]
            
            # We save the cleaned image in a dictionary:
            cleaned_image_name = f"cleaned_bias_{idx}"
            cleaned_images[cleaned_image_name] = cleaned_data
            print(f"Image {img_path} cleaned from cosmis rays and saved as ({cleaned_image_name})")
            print("................................................................................")
    return cleaned_images