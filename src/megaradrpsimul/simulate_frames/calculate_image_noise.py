#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Calculate the total noise in an image considering photon shot noise and readout noise."""

import numpy as np

def calculate_image_noise(data_model_traces, gain_array, readout_noise_median):
    """
    Calculate the total noise in an image considering photon shot noise and readout noise.

    This function estimates the pixel-wise noise in units of ADU (Analog-to-Digital Units),
    combining the contribution from photon statistics (shot noise) and the detector's 
    readout noise. The photon noise is computed as the square root of the signal divided 
    by the gain, and the readout noise is included in quadrature.

    Parameters
    ----------
    data_model_traces : ndarray
        2D array (or broadcastable shape) of the modelled signal in ADU.
    gain_array : ndarray or float
        2D array or scalar representing the gain (electrons per ADU) for each pixel.
        Must be broadcast-compatible with `data_model_traces`.
    readout_noise_median : float or ndarray
        Scalar or array with the median readout noise (in ADU). Can be a single value
        or an array matching the shape of the image.

    Returns
    -------
    noise_array : ndarray
        Array of the same shape as `data_model_traces` containing the estimated noise per pixel,
        in ADU.

    Notes
    -----
    The total noise is computed using:

        noise = sqrt( (signal / gain) + readout_noise**2 )

    where:
    - signal / gain → photon shot noise (in electrons)
    - readout_noise**2 → detector noise (in ADU^2), added in quadrature
    """
    ratio_data_gain = data_model_traces / gain_array
    readout_noise_square = readout_noise_median ** 2
    noise_array = np.sqrt(ratio_data_gain + readout_noise_square)
    return noise_array