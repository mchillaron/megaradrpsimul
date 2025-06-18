#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Calculate the readout noise across the CCD based on top and bottom bias regions."""

import numpy as np
import teareduce as tea

from megaradrpsimul import CCDregions

def calculate_bias_readout_noise(data, data_smoothed, naxis2, naxis1):
    """
    Calculate the readout noise across the CCD based on top and bottom bias regions.

    The function computes the pixel-wise difference between the raw and smoothed
    data, then evaluates the robust standard deviation in predefined top and bottom
    CCD bias regions to estimate readout noise. The resulting noise is propagated
    into a full-sized array, assigning different values to top and bottom sections.

    Parameters
    ----------
    data : ndarray
        The original 2D CCD image data.
    data_smoothed : ndarray
        The smoothed version of the CCD image, used as a reference for bias subtraction.
    naxis2 : int
        Number of rows (Y-axis size) of the CCD image.
    naxis1 : int
        Number of columns (X-axis size) of the CCD image.

    Returns
    -------
    noise_array : ndarray
        A 2D array of the same shape as the input, filled with the estimated
        readout noise values, using different estimates for top and bottom regions.
    """
    data_bias_diff = data - data_smoothed

    regions_kernel = CCDregions.regions_kernel
    # Top CCD region
    top_slice = regions_kernel["bias_topCCD"]["slice2d"].python
    top_data = data_bias_diff[top_slice]
    noise_top = tea.robust_std(top_data)
    # Bottom CCD region
    bottom_slice = regions_kernel["bias_bottomCCD"]["slice2d"].python
    bottom_data = data_bias_diff[bottom_slice]
    noise_bottom = tea.robust_std(bottom_data)

    #Create the readout_noise array
    noise_array = np.full((naxis2, naxis1), noise_top)
    noise_array[0:2106, :] = noise_bottom

    print(f"Robust deviation for TopCCD is: {noise_top} and for BottomCCD is: {noise_bottom}")

    return noise_array