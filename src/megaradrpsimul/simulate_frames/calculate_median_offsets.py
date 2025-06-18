#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Calculate the median image and offset values for a set of CCD images."""

from astropy.io import fits
import numpy as np

from megaradrpsimul import CCDregions


def calculate_median_offsets(list_images, bias_smoothed_median):
    """
    Calculate the median image and offset values for a set of CCD images.

    This function computes the median image from a list of input CCD images,
    and compares overscan regions (top and bottom) against a smoothed bias
    median to evaluate offset consistency. It returns the image median,
    individual image offsets, and the difference between the median image
    and the reference smoothed bias.

    Parameters
    ----------
    list_images : list of Path
        List of file paths to FITS images used for median calculation.
    bias_smoothed_median : ndarray
        A reference image representing the smoothed median bias, used to
        calculate absolute offsets.

    Returns
    -------
    median_image : ndarray
        The median image computed from all images in `list_images`.
    image_median_offsets : list of tuple of float
        List of (top_offset, bottom_offset) tuples representing how much
        each image deviates from the median image in the overscan regions.
    top_bias_median_offset : float
        Offset between the top overscan region of the median image and the bias.
    bottom_bias_median_offset : float
        Offset between the bottom overscan region of the median image and the bias.
    """
    top_bias_constant = np.median(bias_smoothed_median[CCDregions.overscan3_top.python])
    bottom_bias_constant = np.median(bias_smoothed_median[CCDregions.overscan3_bottom.python])

    image_data = []
    image_offsets = []

    for path in list_images:
        with fits.open(path) as hdul:
            data = hdul[0].data
            image_data.append(data)

            top_offset = np.median(data[CCDregions.overscan3_top.python])
            bottom_offset = np.median(data[CCDregions.overscan3_bottom.python])
            image_offsets.append((top_offset, bottom_offset))
            print(f"Individual image || Top constant: {top_offset}, Bottom constant: {bottom_offset}")

    stacked = np.stack(image_data, axis=0)
    median_image = np.median(stacked, axis=0)
    top_median_constant = np.median(median_image[CCDregions.overscan3_top.python])
    bottom_median_constant = np.median(median_image[CCDregions.overscan3_bottom.python])

    top_bias_median_offset = top_bias_constant - top_median_constant
    bottom_bias_median_offset = bottom_bias_constant - bottom_median_constant

    image_median_single_offsets = [
        (top_median_constant - top, bottom_median_constant - bottom)
        for top, bottom in image_offsets
    ]

    print(f"MasterBias || Top constant: {top_bias_constant}, Bottom constant: {bottom_bias_constant}")
    print(f"Median image || Top constant: {top_median_constant}, Bottom constant: {bottom_median_constant}")
    print("Offsets of Median image - Individual image:")
    for i, offsets in enumerate(image_median_single_offsets):
        print(f"Image {i}: Top offset = {offsets[0]}, Bottom offset = {offsets[1]}")
    print("Offset of MasterBias - Median image:")
    print(f"Top offset = {top_bias_median_offset}, Bottom offset = {bottom_bias_median_offset}")

    return median_image, image_median_single_offsets, top_bias_median_offset, bottom_bias_median_offset