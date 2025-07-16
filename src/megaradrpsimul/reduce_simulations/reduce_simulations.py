#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Reduce the simulated images using the MEGARA pipeline"""

import subprocess
from astropy.io import fits
import os
import shutil

from .conversion_into_cube import conversion_into_cube
from .get_step_name import get_step_name
from .step_reduction import step_reduction
from .healing_traces import healing_traces
from .diffuselight_determination import diffuselight_determination

def reduce_simulations(niter, config, nstart, 
                       abs_results_dir, 
                       run_modelmap=False, run_twilight=False, 
                       run_healing=False, run_LRU=False, 
                       run_diffuselight=False, 
                       pixel_size=0.4, 
                       history_line_command=None):
    """Reduce the simulated images using the TEA pipeline.

    Parameters
    ----------
    niter : int
        Number of the current simulation.
    config : dict
        Dictionary containing the configuration for the simulation.
    nstart : int
        Number from which the simulation process starts.
    abs_results_dir : Path instance
        Path to the absolute results directory.
    run_modelmap : bool, optional
        If True, run the ModelMap step (default is False).
    run_twilight : bool, optional
        If True, run the TwilightFlat step (default is False).
    run_healing : bool, optional
        If True, run the healing of the traces (default is False).
    run_LRU : bool, optional
        If True, run the special TraceMap template for LR-U (default is False).
    run_diffuselight : bool, optional
        If True, run the DiffuseLight step (default is False).
    pixel_size : float, optional
        Pixel size in arcseconds for the conversion of RSS into a cube (default is 0.4).
    history_line_command : str, optional
        Command line to be added to the HISTORY keyword of the final_rss.fits file (default is None).
    """
    
    print('we start the reduction process')
    num = niter + nstart
    vph_name = config["VPH"]

    print('........... Step 0: Bias image ...........')
    bias_filename = config["0_Bias"] + '.yaml'
    step_name = get_step_name(bias_filename)
    print('Step name:', step_name)
    print('Step name:', step_name)
    step_reduction(bias_filename, step_name, product_file="master_bias.fits", calib_folder_path="MasterBias")
    
    print('........... Step 1: TraceMap ...........')
    tracemap_filename = config["1_TraceMap"] + '.yaml'
    step_name = get_step_name(tracemap_filename)
    print('Step name:', step_name)

    if run_LRU == True:
        print('Running special TraceMap template for LR-U')
        tracesU_filename = config["master_traces_LRU_20220325_healed"] + '.json'
        command_copy_traces_LRU_list = [
            """cp""",
            f"""{tracesU_filename}""",
            f"""ca3558e3-e50d-4bbc-86bd-da50a0998a48/TraceMap/LCB/{vph_name}/"""
            ]
        print('\033[1m\033[35m ' + f"$ {' '.join(command_copy_traces_LRU_list)}" + '\033[0m\n')
        subprocess.run(command_copy_traces_LRU_list, capture_output=True, text=True)
    else:
        step_reduction(tracemap_filename, step_name, product_file="master_traces.json", calib_folder_path=f"TraceMap/LCB/{vph_name}")
        if run_healing == True:
                print("Healing of the traces required")
                healing_traces(vph_name, step_name)
    
        
    if run_modelmap == True:
        print('........... Step 2: ModelMap ...........')
        modelmap_filename = config["2_ModelMap"] + '.yaml'
        step_name = get_step_name(modelmap_filename)
        step_reduction(modelmap_filename, step_name, product_file="master_model.json", calib_folder_path=f"ModelMap/LCB/{vph_name}")

    print('........... Step 3: WavelengthCalibration ...........')
    wavecalib_filename = config["3_WaveCalib"] + '.yaml'
    step_name = get_step_name(wavecalib_filename)
    print('Step name:', step_name)
    step_reduction(wavecalib_filename, step_name, product_file="master_wlcalib.json", calib_folder_path=f"WavelengthCalibration/LCB/{vph_name}")

    print('........... Step 3: WavelengthCalibration - Check ...........')
    wavecalibcheck_filename = config["3_WaveCalib_check"] + '.yaml'
    step_name = get_step_name(wavecalibcheck_filename)
    step_reduction(wavecalibcheck_filename, step_name)

    print('........... Step 4: FiberFlat ...........')
    fiberflat_filename = config["4_FiberFlat"] + '.yaml'
    step_name = get_step_name(fiberflat_filename)
    print('Step name:', step_name)
    step_reduction(fiberflat_filename, step_name, product_file="master_fiberflat.fits", calib_folder_path=f"MasterFiberFlat/LCB/{vph_name}")

    if run_twilight == True:
        print('........... Step 5: Bias image ...........')
        twilight_filename = config["5_TwilightFlat"] + '.yaml'
        step_name = get_step_name(twilight_filename)
        print('Step name:', step_name)
        step_reduction(twilight_filename, step_name, product_file="master_twilightflat.fits", calib_folder_path=f"MasterTwilightFlat/LCB/{vph_name}")
        
    print('........... Step 6: LCBadquisition ...........')
    lcbadquisition_filename = config["6_LcbAdquisition"] + '.yaml'
    step_name = get_step_name(lcbadquisition_filename)
    print('Step name:', step_name)
    step_reduction(lcbadquisition_filename, step_name)

    print('........... Step 7: StandardStar ...........')
    standardstar_filename = config["7_StandardStar"] + '.yaml'
    step_name = get_step_name(standardstar_filename)
    print('Step name:', step_name)
    step_reduction(standardstar_filename, step_name, product_file="master_sensitivity.fits", calib_folder_path=f"MasterSensitivity/LCB/{vph_name}")
    
    print('........... Step 8: Reduce LCB ...........')
    reduce_filename = config["8_LcbImage"] + '.yaml'
    step_name_8, extraction_offset = get_step_name(reduce_filename, extraction_offset=True)
    print('Step name:', step_name_8)
    print('Extraction offset:', extraction_offset)
    step_reduction(reduce_filename, step_name_8)

    if run_diffuselight == True:
        print("Correcting for diffuse light")
        diffuselight_determination(vph_name, step_name_8, run_LRU, run_healing, extraction_offset)

        # Now we run again the 8th step of the reduction process, but using the new yaml file:
        diffuselight_filename = config["8_LcbImage_diffuse_light"] + '.yaml'
        step_name_difflight = get_step_name(diffuselight_filename)
        print('Step name:', step_name_difflight)
        step_reduction(diffuselight_filename, step_name_difflight)

        # Change the name of the final_rss.fits file to the name of the simulation
        original_file = f"obsid{step_name_difflight}_results/final_rss.fits"
        new_file = f"obsid{step_name_difflight}_results/final_rss_{num:04d}.fits"
    else:
        # Change the name of the final_rss.fits file to the name of the simulation
        original_file = f"obsid{step_name_8}_results/final_rss.fits"
        new_file = f"obsid{step_name_8}_results/final_rss_{num:04d}.fits"
    
    print('End of the reduction process')

    if os.path.exists(original_file):
        os.rename(original_file, new_file)
        with fits.open(new_file, mode='update') as hdul:
            hdul[0].header.add_history(history_line_command)
            hdul.flush() 
        print('Command line added to final_rss.fits HISTORY keyword')
    else:
        print(f"The file {original_file} does not exist.")


    dest_file = os.path.join(abs_results_dir, os.path.basename(new_file))
    # If there is already a file with the same name in the results directory, we overwrite it
    if os.path.exists(dest_file):
        os.remove(dest_file)  # Eliminate the existing file to move the new one
        print(f"Existing {dest_file} being deleted.")

    shutil.move(new_file, abs_results_dir)
    print(f"Moved {new_file} to {abs_results_dir}")

    # Conversion of the RSS into a cube with the pixel size specified
    print('Converting the RSS into a cube with a pixel size of', pixel_size, 'arcseconds')
    cube_name = f"final_cube_{num:04d}.fits"
    dest_cube_file = os.path.join(abs_results_dir, cube_name)
    print(dest_cube_file)

    conversion_into_cube(pixel_size, dest_cube_file, dest_file)


    print('the simulation and reduction processes have finished')
    print('................................................................................')