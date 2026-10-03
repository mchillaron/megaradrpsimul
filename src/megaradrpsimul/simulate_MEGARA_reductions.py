#
# Copyright 2025-2026 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from pathlib import Path
from rich.console import Console
from rich.logging import RichHandler

import argparse
import logging
import os
import re
import shutil
import subprocess
import yaml

from .find_ob_directories import find_ob_directories
from .reduce_simulations.reduce_simulations import reduce_simulations
from .simulate_frames.simulate_frames import simulate_frames

logger = logging.getLogger(__name__)
console = Console()

def get_num_start(results_dir):
    """
    Determine the starting index for naming new simulation output files.

    If the `results_dir` does not exist, it is created and the starting index is set to 1.
    If it exists, the function scans for existing files matching the pattern 
    'final_rss_####.fits', determines the highest index, and returns the next one.

    Parameters
    ----------
    results_dir : Path
        Path instance pointing to the simulation results directory.

    Returns
    -------
    int
        The starting index for naming new simulation output files, ensuring continuity
        and avoiding overwriting existing results.
    """
    if not results_dir.is_dir():
        results_dir.mkdir(exist_ok=True)
        logger.info('New simulation results directory created: %s', results_dir)
        return 1

    logger.info('previous simulation results directory will be used to save final_rss.fits')
    existing_files = os.listdir(results_dir)
    pattern = re.compile(r"final_rss_(\d{4})\.fits")
    indices = [
        int(match.group(1))
        for filename in existing_files
        if (match := pattern.match(filename))
    ]
    return max(indices) + 1 if indices else 1


def ask_confirmation(nstart, 
                    nsimul, 
                    results_dir):
    """
    Prompt the user to confirm whether to proceed with the simulation run.

    Parameters
    ----------
    nstart : int
        The starting index for the simulation file naming.
    nsimul : int
        The number of simulations to be executed.
    results_dir : Path
        Directory where the output files will be saved.

    Returns
    -------
    bool
        True if the user confirms ('y'), False otherwise.
    """
    end_index = nstart + nsimul - 1
    msg = (
        f"You are going to run {nsimul} simulations.\n"
        f"Files will be saved at: '{results_dir}'\n"
        f"numbered from: final_rss_{nstart:04d}.fits "
        f"to: final_rss_{end_index:04d}.fits\n"
        f"Would you like to continue? (y/n): "
    )
    response = input(msg).strip().lower()
    return response == 'y'


def simulate_MEGARA_reductions(ob, 
                               config_file,
                               nsimul=1, 
                               run_modelmap=False, 
                               run_twilight=False, 
                               pixel_size=0.4,
                               history_line_command=None):
   
    """Simulate MEGARA reductions for a given observation.

    This function simulates the MEGARA reductions for a given observation directory.
    It creates a work directory, copies the necessary files from the MEGARA directory,
    and simulates the frames using the `simulate_frames` function.
    It also handles the creation and deletion of the work directory.
    It is assumed that the MEGARA directory contains the necessary YAML files and data
    files for the simulation.
    
    Parameters
    ----------
    ob : Path instance
        Path to the observation directory.
    config_file : str
        Name of the configuration file for the simulation.
    nsimul : int, optional
        Number of simulations to perform (default is 1).
    run_modelmap : bool, optional
        If True, run the ModelMap step (default is False).
    run_twilight : bool, optional
        If True, run the Twilight step (default is False).
    pixel_size : float, optional
        Pixel size in arcseconds for the conversion of RSS into a cube (default is 0.4).
    history_line_command : str, optional
        Command line history for the simulation (default is None).

    """

    logger.info('Simulating MEGARA reductions for: %s', ob)

    # Check if MEGARA directory exists
    megara_dir = ob / 'MEGARA'
    if not megara_dir.is_dir():
        raise ValueError(f"MEGARA directory does not exist for {ob}.")

    config_file_path = ob / config_file

    if not config_file_path.is_file():
        raise ValueError(
            f"Configuration file {config_file} does not exist in {ob}."
        )

    logger.info(f"Using configuration file: %s", config_file_path)

    with open(config_file_path, "r") as f:
        config = yaml.safe_load(f)

    # Validation of the basic configuration file structure
    if not isinstance(config, dict):
        raise ValueError(
            f"The configuration file {config_file} must contain "
            "a YAML mapping of keys and values."
        )

    mandatory_keys = {
        "VPH",
        "0_Bias",
        "1_TraceMap",
        "3_WaveCalib",
        "3_WaveCalib_check",
        "4_FiberFlat",
        "6_LcbAdquisition",
        "7_StandardStar",
        "8_LcbImage",
    }

    optional_keys = {
        "2_ModelMap",
        "5_TwilightFlat",
        "8_LcbImage_diffuse_light",
        "8_LcbImage_cleaned",
        "healing",
        "master_traces_LRU_20220325_healed",
    }

    expected_keys = mandatory_keys | optional_keys    # the total number of keys that should be in the config file

    # Check for unexpected or missing keys
    config_keys = set(config)

    unexpected = config_keys - expected_keys
    missing = expected_keys - config_keys

    if unexpected or missing:
        raise ValueError(
            "Unexpected or missing keys in configuration file.\n"
            f"Unexpected: {sorted(unexpected)}\n"
            f"Missing: {sorted(missing)}"
        )

    # Check types
    for key, value in config.items():
        if value is not None and not isinstance(value, str):
            raise ValueError(
                f"The value associated with '{key}' must be "
                f"a string or empty, not {type(value).__name__}."
            )

    # Check mandatory keys for non-empty values
    for key in mandatory_keys:
        value = config.get(key)
        if value is None or not value.strip():
            raise ValueError(
                f"The key '{key}' in {config_file} "
                "must have a non-empty value."
            )

    # Check specific conditions for optional keys based on other parameters
    conditional_keys = {
        "2_ModelMap": run_modelmap,
        "5_TwilightFlat": run_twilight,
    }

    for key, enabled in conditional_keys.items():

        if enabled and not config.get(key):
            raise ValueError(
                f"The key '{key}' must have a non-empty value "
                "because the corresponding reduction step has been called in the command line."
            )

    work_dir = ob / 'work'
    work_megara_dir = work_dir / 'MEGARA'
    logger.debug("The established working directory for simulations is: %s",work_dir)
    
    # Define the directory to store the final_rss.fits
    results_dir = ob / 'simulation_results'
    nstart = get_num_start(results_dir)
    logger.info("Simulations will start from number: %d", nstart)
    abs_results_dir = os.path.abspath(results_dir)
    
    if not ask_confirmation(nstart, nsimul, results_dir):
        logger.info("Action aborted. No changes were made.")
        return

    for i in range(nsimul):
        logger.info("Simulation and reduction number: %d", i+1)

        # Check if work directory exists
        if work_dir.is_dir():
            logger.debug('work directory already exists')
            
            # Delete work directory and everything inside it
            shutil.rmtree(work_dir, ignore_errors=True)
            logger.debug('work directory %s deleted', work_dir)

        # Create the directory once it has been deleted
        work_dir.mkdir()
        work_megara_dir.mkdir()
        logger.debug('work/MEGARA/ directory created')
        
        # From MEGARA/ directory, copy all files *.yaml to work/MEGARA/
        for file in megara_dir.glob('*.yaml'):
            shutil.copy(file, work_megara_dir / file.name)
            logger.debug('copying %s to %s', file.name, work_megara_dir)

        # create the calibration directory in work/MEGARA/ using inittree
        logger.info("Initializing calibration directory in work/MEGARA/ using megaradrp-inittree...")
        subprocess.run(["megaradrp-inittree"], cwd=work_megara_dir, check=True,)
        logger.debug("Calibration directory initialized in %s", work_megara_dir)

        # Now, we copy the data/ directory from MEGARA to work and exclude MEGARA raw images
        data_dir = megara_dir / 'data'
        data_work_dir = work_megara_dir / 'data'

        ignored_data_files = []
        for file in data_dir.iterdir():
            if file.name.startswith('0') and file.name.endswith('.fits'):
                ignored_data_files.append(file.name)
        logger.debug("Files to be ignored during data copy: %s", ignored_data_files)

        shutil.copytree(data_dir, data_work_dir, ignore=shutil.ignore_patterns(*ignored_data_files))
        logger.info("Auxiliary data files copied to %s", data_work_dir)

        # activate healing the traces if the file healing.yaml exists
        if config.get("healing"):
            healing_filename = config["healing"] + '.yaml'
            healing_yaml_path = work_megara_dir / healing_filename

            run_healing = healing_yaml_path.is_file()
            logger.debug("run_healing is set to: %s",run_healing)

            if not run_healing:
                raise FileNotFoundError(f"The file needed for 'healing' step: {healing_yaml_path} was not found.")
        else:
            run_healing = False
            logger.info("No healing of traces for this simulation and reduction.")

        # activate LRU if the file master_traces_LRU_20220325_healed.json exists
        if config.get("master_traces_LRU_20220325_healed"):
            master_traces_name = config["master_traces_LRU_20220325_healed"] + ".json"
            master_traces_LRU_path = megara_dir / master_traces_name

            if not master_traces_LRU_path.is_file():
                raise FileNotFoundError(f"The configured healed master traces file was not found: {master_traces_LRU_path}")

            run_LRU = True
            destination = work_megara_dir / master_traces_name
            shutil.copy(master_traces_LRU_path, destination)

            logger.info("Healed LR-U master traces enabled.")
            logger.debug("Copying %s to %s", master_traces_LRU_path, destination)

        else:
            run_LRU = False
            logger.debug("No healed LR-U master traces provided.")

        
        # activate diffuse light if the file 8_LcbImage_diffuse_light.yaml exists
        if config.get("8_LcbImage_diffuse_light"):
            diffuse_light_filename = config["8_LcbImage_diffuse_light"] + ".yaml"
            diffuse_light_yaml_path = work_megara_dir / diffuse_light_filename

            run_diffuselight = diffuse_light_yaml_path.is_file()
            logger.debug("run_diffuselight is set to: %s", run_diffuselight)

            if not run_diffuselight:
                raise FileNotFoundError(
                    f"The file needed for the diffuse-light step "
                    f"was not found: {diffuse_light_yaml_path}"
                )

            logger.info("Diffuse-light correction enabled.")
        else:
            run_diffuselight = False
            logger.info("Diffuse-light correction disabled.")


        # activate Cosmic-ray cleaning using numina/crmasks
        if config.get("8_LcbImage_cleaned"):
            cleaned_filename = config["8_LcbImage_cleaned"] + ".yaml"
            cleaned_yaml_path = work_megara_dir / cleaned_filename
            crmasks_path = data_work_dir / "crmasks.fits"

            missing_files = []
            if not cleaned_yaml_path.is_file():
                missing_files.append(str(cleaned_yaml_path))
            if not crmasks_path.is_file():
                missing_files.append(str(crmasks_path))
            if missing_files:
                raise FileNotFoundError(
                    "Cosmic-ray cleaning was requested, but the following "
                    "required file(s) were not found:\n"
                    + "\n".join(f"  - {file}" for file in missing_files))

            run_crclean = True

            logger.info("Cosmic-ray cleaning enabled.")
            logger.debug("Cosmic-ray cleaning YAML: %s", cleaned_yaml_path)
            logger.debug("Cosmic-ray mask: %s", crmasks_path)

        else:
            run_crclean = False
            logger.info("Cosmic-ray cleaning disabled.")

        #------------------------SIMULATION OF FRAMES------------------------
        # At this point, we are ready to start simulating images
        simulate_frames(megara_dir, data_work_dir, work_megara_dir, config)

        logger.info('All the frames have been simulated')


        #------------------------REDUCTION OF SIMULATED FRAMES---------------

        # Save the current working directory
        original_dir = Path.cwd()

        # Reduction must be run from work/MEGARA
        reduction_dir = work_megara_dir.resolve()

        if not reduction_dir.is_dir():
            raise FileNotFoundError(f"MEGARA reduction directory not found: {reduction_dir}")

        os.chdir(reduction_dir)
        logger.info("Starting reduction of simulated frames.")
        logger.debug("Working directory changed to: %s", Path.cwd())

        reduce_simulations(i, config, nstart, abs_results_dir,
                          run_modelmap, run_twilight, 
                          run_healing, run_LRU, run_diffuselight, run_crclean,
                          pixel_size, history_line_command)

        # we go back to the directory where the script was executed:
        os.chdir(original_dir)   # Go back to original directory
        logger.info('The directory has been changed to: %s', os.getcwd())     



def main():
    parser = argparse.ArgumentParser(description='Simulate MEGARA reductions.')
    parser.add_argument('--obj', type=str, help='Object directories must start with "obj_". Default: obj_*', default='obj_*')
    parser.add_argument('--vph', type=str, help='VPH directories must start with "VPH_". Default: VPH_*', default='VPH_*')
    parser.add_argument('-c', '--config-file', type=str, help='Name of the YAML configuration file.', default='config_simulation.yaml')
    parser.add_argument('-n', '--num-simul', type=int, help='Number of simulations to perform.', default=1)
    parser.add_argument('--run-modelmap', action='store_true', help='Run ModelMap step.')
    parser.add_argument('--run-twilight', action='store_true', help='Run Twilight step.')
    parser.add_argument('--pixel-size', type=float, help='Pixel size in arcseconds for the conversion of RSS into a cube.', default=0.4)
    parser.add_argument("--log-level", type=str, choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO", help="Set the logging level. Default: INFO.")
    args = parser.parse_args()

    # logging configuration
    if args.log_level in ["DEBUG", "WARNING", "ERROR", "CRITICAL"]:

        format_log = "%(message)s"
        handlers = [
            RichHandler(
                console=console,
                show_time=False,
                show_level=True,
                show_path=True,
                markup=True,
                rich_tracebacks=True,
            )
        ]

    else:

        format_log = "%(message)s"
        handlers = [
            RichHandler(
                console=console,
                show_time=False,
                show_level=False,
                show_path=False,
                markup=True,
            )
        ]

    logging.basicConfig(level=args.log_level, format=format_log, handlers=handlers)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    obj = args.obj
    vph = args.vph

    config_file = args.config_file
    run_modelmap = args.run_modelmap
    run_twilight = args.run_twilight
    pixel_size = args.pixel_size

    logger.info("Running ModelMap: %s", run_modelmap)
    logger.info("Running Twilight: %s", run_twilight)

    if obj is not None:
        logger.debug("Object search pattern: %s", obj)
    
    if vph is not None:
        logger.debug("VPH pattern: %s", vph)
    
    history_line_command = (
            f"$ python simulate_MEGARA_reductions.py "
            f"--obj {obj} --vph {vph} "
            f"--config-file {config_file} "
            f"--num-simul {args.num_simul} "
            f"{'--run-modelmap' if run_modelmap else ''} "
            f"{'--run-twilight' if run_twilight else ''}"
            f" --pixel-size {pixel_size} "
            f"--log-level {args.log_level}"
        )
    logger.debug("History line command: %s", history_line_command)

    # Validate arguments
    logger.debug("Validating arguments...")
    
    if args.num_simul < 1:
        raise ValueError("Number of simulations must be at least 1.")
    else:
        logger.info("Number of simulations: %s", args.num_simul)
        nsimul = args.num_simul

    if pixel_size <= 0:
        raise ValueError("Pixel size must be a positive number.")
    else:
        logger.info("Pixel size: %s arcseconds", pixel_size)

    # locate observations
    ob_list = find_ob_directories(obj_pattern=obj, vph_pattern=vph,)

    logger.debug("Found observation directories: %s", ob_list)
    if not ob_list:
        logger.warning(
            "No observation directories found from %s "
            "using object pattern '%s' and VPH pattern '%s'.",
            Path.cwd(),
            obj,
            vph,
        )
        return

    logger.info("Found %d observation(s)",len(ob_list))

    # Run simulations
    for ob in ob_list:
        logger.debug("Processing observation: %s", ob)
        simulate_MEGARA_reductions(ob, config_file, nsimul, run_modelmap, run_twilight, pixel_size, history_line_command)
        

if __name__ == "__main__":
    main()
