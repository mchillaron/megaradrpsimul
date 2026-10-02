#
# Copyright 2025-2026 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Reduce the simulated images using the MEGARA pipeline"""

from astropy.io import fits
from pathlib import Path

import subprocess
import logging
import os
import shutil

from .conversion_into_cube import conversion_into_cube
from .get_calibration_label import get_calibration_label
from .get_step_name import get_step_name
from .step_reduction import step_reduction
from .healing_traces import healing_traces
from .diffuselight_determination import diffuselight_determination

logger = logging.getLogger(__name__)

def reduce_simulations(niter, config, nstart, abs_results_dir,
                       run_modelmap=False, run_twilight=False, 
                       run_healing=False, run_LRU=False, 
                       run_diffuselight=False, run_crclean=False, 
                       pixel_size=0.4, 
                       history_line_command=None):
    
    """Reduce the simulated images using MEGARA DRP.

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
    run_crclean : bool, optional
        If True, run the cosmic ray cleaning step (default is False).
    pixel_size : float, optional
        Pixel size in arcseconds for the conversion of RSS into a cube (default is 0.4).
    history_line_command : str, optional
        Command line to be added to the HISTORY keyword of the final_rss.fits file (default is None).
    """
    
    logger.info('Starting the reduction process')

    num = niter + nstart

    # get VPH name and calibration label
    vph_name = config["VPH"]

    sim_data_dir = Path("data")
    calib_dir = Path("calibrations") / "MEGARA"

    logger.info("Simulated data directory: %s", sim_data_dir)
    logger.info("Calibration directory: %s", calib_dir)

    calib_label_dir = get_calibration_label(sim_data_dir, calib_dir)
    logger.info("Complete calibration directory with INSCONF label: %s", calib_label_dir)

    print("\n")
    logger.info("[bold green]Step 0: Bias image [/bold green]")
    bias_filename = config["0_Bias"] + '.yaml'
    step_name = get_step_name(bias_filename)
    logger.debug('Step name: %s', step_name)
    step_reduction(bias_filename, step_name, calib_label_dir,  
                   product_file="master_bias.fits", 
                   calib_folder_path="MasterBias")


    print("\n")
    logger.info("[bold green]Step 1: TraceMap [/bold green]")
    tracemap_filename = config["1_TraceMap"] + '.yaml'
    step_name = get_step_name(tracemap_filename)
    logger.debug('Step name: %s', step_name)

    if run_LRU:
        logger.info('Running special TraceMap template for LR-U')
        tracesU_filename = config["master_traces_LRU_20220325_healed"] + '.json'
        shutil.copy2(tracesU_filename,calib_label_dir / "TraceMap" / "LCB" / vph_name)
        logger.debug("[bold red]Copying %s to %s[/bold red]", tracesU_filename,
                    calib_label_dir / "TraceMap" / "LCB" / vph_name)
    else:
        step_reduction(tracemap_filename, step_name, calib_label_dir,
                       product_file="master_traces.json", 
                       calib_folder_path=f"TraceMap/LCB/{vph_name}")
        
        if run_healing:
                logger.info("Healing of the traces required")
                healing_traces(vph_name, step_name, calib_label_dir)
        
    if run_modelmap:
        print("\n")
        logger.info("[bold green]Step 2: ModelMap - Check [/bold green]")
        modelmap_filename = config["2_ModelMap"] + '.yaml'
        step_name = get_step_name(modelmap_filename)
        logger.debug('Step name: %s', step_name)
        step_reduction(modelmap_filename, step_name, calib_label_dir,
                       product_file="master_model.json", calib_folder_path=f"ModelMap/LCB/{vph_name}")

    print("\n")
    logger.info("[bold green]Step 3: Wavelength Calibration [/bold green]")
    wavecalib_filename = config["3_WaveCalib"] + '.yaml'
    step_name = get_step_name(wavecalib_filename)
    logger.debug("Step name: %s", step_name)
    step_reduction(wavecalib_filename, step_name, calib_label_dir,
                   product_file="master_wlcalib.json", calib_folder_path=f"WavelengthCalibration/LCB/{vph_name}")

    print("\n")
    logger.info("[bold green]Step 3: Wavelength Calibration - Check [/bold green]")
    wavecalibcheck_filename = config["3_WaveCalib_check"] + '.yaml'
    step_name = get_step_name(wavecalibcheck_filename)
    step_reduction(wavecalibcheck_filename, step_name, calib_label_dir)

    print("\n")
    logger.info("[bold green]Step 4: FiberFlat[/bold green]")
    fiberflat_filename = config["4_FiberFlat"] + '.yaml'
    step_name = get_step_name(fiberflat_filename)
    logger.debug('Step name: %s', step_name)
    step_reduction(fiberflat_filename, step_name, calib_label_dir,
                   product_file="master_fiberflat.fits", calib_folder_path=f"MasterFiberFlat/LCB/{vph_name}")

    if run_twilight == True:
        print("\n")
        logger.info("[bold green]Step 5: Twilight Flat [/bold green]")
        twilight_filename = config["5_TwilightFlat"] + '.yaml'
        step_name = get_step_name(twilight_filename)
        logger.debug('Step name: %s', step_name)
        step_reduction(twilight_filename, step_name, calib_label_dir,
                       product_file="master_twilightflat.fits", calib_folder_path=f"MasterTwilightFlat/LCB/{vph_name}")

    print("\n")
    logger.info("[bold green]Step 6: LCB Adquisition [/bold green]")
    lcbadquisition_filename = config["6_LcbAdquisition"] + '.yaml'
    step_name = get_step_name(lcbadquisition_filename)
    logger.debug('Step name: %s', step_name)
    step_reduction(lcbadquisition_filename, step_name, calib_label_dir)

    print("\n")
    logger.info("[bold green]Step 7: Standard Star [/bold green]")
    standardstar_filename = config["7_StandardStar"] + '.yaml'
    step_name = get_step_name(standardstar_filename)
    logger.debug('Step name: %s', step_name)
    step_reduction(standardstar_filename, step_name, calib_label_dir,
                   product_file="master_sensitivity.fits", calib_folder_path=f"MasterSensitivity/LCB/{vph_name}")

    print("\n")
    logger.info("[bold green]Step 8: Reduce LCB [/bold green]")
    reduce_filename = config["8_LcbImage"] + '.yaml'
    step_name_8, extraction_offset = get_step_name(reduce_filename, extraction_offset=True)
    logger.debug('Step name: %s', step_name_8)
    logger.debug('Extraction offset: %s', extraction_offset)
    step_reduction(reduce_filename, step_name_8, calib_label_dir)

    # By default, the final product comes from the standard Step 8
    final_step_name = step_name_8

    # Apply diffuse light correction if requested
    if run_diffuselight:
        print("\n")
        logger.info("[bold green]Correcting for diffuse light[/bold green]")

        diffuselight_determination(vph_name, step_name_8, calib_label_dir,
                                run_LRU, run_healing, extraction_offset)

        diffuselight_filename = (config["8_LcbImage_diffuse_light"] + ".yaml")
        step_name_difflight = get_step_name(diffuselight_filename)
        logger.debug("Diffuse-light reduction step name: %s", step_name_difflight)

        step_reduction(diffuselight_filename, step_name_difflight, calib_label_dir)
        final_step_name = step_name_difflight

    # Clean cosmic rays if requested
    if run_crclean:
        print("\n")
        logger.info("[bold green]Applying cosmic-ray cleaning[/bold green]")

        crclean_filename = (config["8_LcbImage_cleaned"] + ".yaml")
        step_name_crclean = get_step_name(crclean_filename)
        logger.debug("Cosmic-ray cleaning step name: %s", step_name_crclean)

        step_reduction(crclean_filename, step_name_crclean, calib_label_dir)
        final_step_name = step_name_crclean

    logger.info("End of the reduction process")

    # Change the name of the final_rss.fits file to the name of the simulation
    original_file = Path(f"obsid{final_step_name}_results/final_rss.fits")
    new_file = Path(f"obsid{final_step_name}_results/final_rss_{num:04d}.fits")

    if not original_file.is_file():
        raise FileNotFoundError(f"Final RSS file was not found: {original_file}")
    
    original_file.rename(new_file)
    logger.debug("Renamed %s to %s", original_file, new_file)

    # Add command line to FITS HISTORY
    with fits.open(new_file, mode="update") as hdul:
        hdul[0].header.add_history(history_line_command)
        hdul.flush()

    logger.debug("Command line added to HISTORY keyword of %s", new_file)

    # Moving final file to the results directory
    dest_file = (Path(abs_results_dir) / new_file.name)

    if dest_file.exists():
        logger.warning("Existing result will be overwritten: %s", dest_file)
        dest_file.unlink()

    shutil.move(str(new_file), str(dest_file))
    logger.info("Final RSS saved to %s", dest_file)

    # Conversion of the RSS into a cube with the pixel size specified
    logger.info('Converting the RSS into a cube with a pixel size of %s arcseconds', pixel_size)
    cube_name = f"final_cube_{num:04d}.fits"
    dest_cube_file = (Path(abs_results_dir) / cube_name)

    if dest_cube_file.exists():
        logger.warning("Existing cube will be overwritten: %s", dest_cube_file)
        dest_cube_file.unlink()

    conversion_into_cube(pixel_size, dest_cube_file, dest_file)

    logger.info("Final cube saved to %s", dest_cube_file)

    logger.info("[bold green]Simulation and reduction completed successfully.[/bold green]")
    print('................................................................................')