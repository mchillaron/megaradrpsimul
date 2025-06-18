#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Simulate frames for a specific step (e.g., TraceMap, ArcCalibration)."""

from astropy import units as u
from tqdm import tqdm
import numpy as np
import teareduce as tea

from megaradrpsimul import CCDregions
from .calculate_image_noise import calculate_image_noise
from .calculate_median_offsets import calculate_median_offsets
from .simulation_run_save import simulation_run_save

def simulate_frames_step(list_images, 
                        naxis1, naxis2,
                        data_work_dir,
                        name_step, 
                        bias_smoothed_median, 
                        gain_array, 
                        readout_noise_median):
        
        """Simulate frames for a specific step (e.g., TraceMap, ArcCalibration).

        This function processes a list of images, calculates the median image,
        and generates corrected data models for each image. It also calculates
        the noise array based on the data model, gain, and readout noise.
        Parameters
        ----------
        list_images : list of Path instances
            List of image paths to be processed.
        naxis1 : int
            Number of columns (X-axis size) of the CCD image.
        naxis2 : int
            Number of rows (Y-axis size) of the CCD image.
        data_work_dir : Path instance
            Path to the data work directory where the simulated images will be saved.
        name_step : str
            Name of the step for which the frames are being simulated (e.g., 'TraceMap').
        bias_smoothed_median : ndarray
            A reference image representing the smoothed median bias, used to
            calculate absolute offsets.
        gain_array : ndarray
            2D array representing the gain (electrons per ADU) for each pixel.
        readout_noise_median : ndarray
            2D array with the median readout noise (in ADU) for each pixel. 

        """
        median_image, image_single_median_offsets, top_bias_median_offset, bottom_bias_median_offset = calculate_median_offsets(list_images, bias_smoothed_median)
        print(f'median image for {name_step} frames calculated')

        data_model = median_image - bias_smoothed_median
        data_model_corrected = data_model.copy()
        data_model_corrected[CCDregions.topCCD_full.python] += top_bias_median_offset 
        data_model_corrected[CCDregions.bottomCCD_full.python] += bottom_bias_median_offset
        data_model_corrected[data_model_corrected < 0] = 0
        print(f'Corrected {name_step} data model obtained for generator')
        
        noise_array = calculate_image_noise(data_model_corrected, gain_array, readout_noise_median)
        
        generators = {}
        for idx, (image_path, (top_offset, bottom_offset)) in enumerate(zip(list_images, image_single_median_offsets)):
            # create a copy of the bias smoothed median and apply the correction of the offsets
            bias_smoothed_median_corrected = np.array(bias_smoothed_median, copy=True)
            correction_top = np.uint16(top_bias_median_offset + top_offset)
            correction_bottom = np.uint16(bottom_bias_median_offset + bottom_offset)
            bias_smoothed_median_corrected[CCDregions.topCCD_full.python] -= correction_top            
            bias_smoothed_median_corrected[CCDregions.bottomCCD_full.python] -= correction_bottom 

            # now I calculate the median of the corrected bias just to print it and compare with the original values of the images.
            top_bias_corrected_constant = np.median(bias_smoothed_median_corrected[CCDregions.regions_kernel["bias_topCCD"]["slice2d"].python])
            bottom_bias_corrected_constant = np.median(bias_smoothed_median_corrected[CCDregions.regions_kernel["bias_bottomCCD"]["slice2d"].python])
            print(f'The values of the corrected MasterBias are: {top_bias_corrected_constant} and {bottom_bias_corrected_constant}')

            # create the new generator
            image_generator = tea.SimulateCCDExposure(
                naxis1=naxis1, naxis2=naxis2, bitpix=16,
                bias=bias_smoothed_median_corrected * u.adu,
                gain=gain_array * (u.electron/u.adu),
                readout_noise=noise_array * u.adu,
                flatfield=1,
                dark=0 * u.adu,
                data_model=data_model_corrected * u.adu,
            )

            # Save the generator as "generator_{idx}"
            generator_name = f"generator_{idx}"
            generators[generator_name] = image_generator

            print(f'{name_step} generator created')
    
        # run the simulation for each image
        for idx, img_path in tqdm(enumerate(list_images), total=len(list_images), desc=f"Simulating {name_step} images", unit="img"):
            generator_name = f"generator_{idx}"
            generator = generators[generator_name]
            simulation_run_save(data_work_dir, img_path, generator, type_image='object')
