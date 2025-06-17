# Simulate MEGARA's reduction final_rss.fits
# Date: 28-05-2025
# María Chillarón
# Universidad Complutense de Madrid

from pathlib import Path
from .reduce_simulations import reduce_simulations
from .simulate_frames import simulate_frames

import argparse
import os
import re
import shutil
import subprocess
 

# -------------------------------------------------------------------------
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
        print('New simulation results directory created:', results_dir)
        return 1

    print('previous simulation results directory will be used to save final_rss.fits')
    existing_files = os.listdir(results_dir)
    pattern = re.compile(r"final_rss_(\d{4})\.fits")
    indices = [
        int(match.group(1))
        for filename in existing_files
        if (match := pattern.match(filename))
    ]
    return max(indices) + 1 if indices else 1

# -------------------------------------------------------------------------
def ask_confirmation(nstart, nsimul, results_dir):
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

# -------------------------------------------------------------------------
def simulate_MEGARA_reductions(ob, nsimul=1, run_modelmap=False, run_twilight=False, history_line_command=None):
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
    nsimul : int, optional
        Number of simulations to perform (default is 1).
    run_modelmap : bool, optional
        If True, run the ModelMap step (default is False).
    run_twilight : bool, optional
        If True, run the Twilight step (default is False).
    history_line_command : str, optional
        Command line history for the simulation (default is None).

    """
    
    print('Simulating MEGARA reductions for:', ob)

    # Check if MEGARA directory exists
    megara_dir = ob / 'MEGARA'
    if not megara_dir.is_dir():
        raise ValueError(f"MEGARA directory does not exist for {ob}.")
    
    work_dir = ob / 'work'
    work_megara_dir = work_dir / 'MEGARA'
    print(work_dir)
    
    # Define the directory to store the final_rss.fits
    results_dir = ob / 'simulation_results'
    nstart = get_num_start(results_dir)
    print(f'simulations will start from number: {nstart}')
    abs_results_dir = os.path.abspath(results_dir)
    
    if not ask_confirmation(nstart, nsimul, results_dir):
        print("Action aborted. No changes were made.")
        return

    for i in range(nsimul):
        print(f"Simulation and reduction number: {i+1}")

        # Check if work directory exists
        if work_dir.is_dir():
            print('work directory already exists')
            # Delete work directory and everything inside it
            shutil.rmtree(work_dir, ignore_errors=True)
            print(f'work directory {work_dir} deleted')

        # Create the directory once it has been deleted
        work_dir.mkdir()
        work_megara_dir.mkdir()
        print('work/MEGARA/ directory created')
        
        # From MEGARA/ directory, copy all files *.yaml to work/MEGARA/
        for file in megara_dir.glob('*.yaml'):
            shutil.copy(file, work_megara_dir / file.name)
            print('copying', file.name, 'to', work_megara_dir)

        # Copy the ca3558e3-e50d-4bbc-86bd-da50a0998a48 tree but empty
        calibration_dir = megara_dir / 'ca3558e3-e50d-4bbc-86bd-da50a0998a48'
        calibration_work_dir = work_megara_dir / 'ca3558e3-e50d-4bbc-86bd-da50a0998a48'
        shutil.copytree(calibration_dir, calibration_work_dir)
        print('copying', calibration_dir, 'to', calibration_work_dir)
        
        for dirpath, dirnames, filenames in os.walk(calibration_work_dir):
            for filename in filenames:
                if filename not in ['master_bpm.fits'] and not filename.endswith('.lis'):
                    file_path = os.path.join(dirpath, filename)
                    try:
                        os.remove(file_path)
                        print('removing', file_path)
                    except Exception as e:
                        print(f"Error removing {file_path}: {e}")

        # If the file healing.yaml exists, set run_healing to True
        healing_yaml_path = work_megara_dir / 'healing.yaml'
        run_healing = healing_yaml_path.is_file()
        print('run_healing is set to:', run_healing)

        # If the file master_traces_LRU_20220325_healed.json exists, copy it to megara_work_dir and set run_LRU to True
        master_traces_LRU_path = megara_dir / 'master_traces_LRU_20220325_healed.json'
        run_LRU = master_traces_LRU_path.is_file()
        if run_LRU:
            shutil.copy(master_traces_LRU_path, work_megara_dir / 'master_traces_LRU_20220325_healed.json')
            print('copying', master_traces_LRU_path, 'to', work_megara_dir / 'master_traces_LRU_20220325_healed.json')
        else:
            print('master_traces_LRU_20220325_healed.json does not exist, run_LRU is set to False')

        # Now, we copy the data/ directory from MEGARA to work
        data_dir = megara_dir / 'data'
        data_work_dir = work_megara_dir / 'data'

        ignored_data_files = []
        for file in data_dir.iterdir():
            if file.name.startswith('0') and file.name.endswith('.fits'):
                ignored_data_files.append(file.name)

        # Copy the data directory to work, ignoring the files in ignored_data_files
        shutil.copytree(data_dir, data_work_dir, ignore=shutil.ignore_patterns(*ignored_data_files))
        print('copying', data_dir, 'to', data_work_dir)
        print('ignoring files:', ignored_data_files)

        # At this point, we are ready to start simulating images
        simulate_frames(megara_dir, data_work_dir, work_megara_dir)

        print('All the frames have been simulated')

        # Now we change the directory to the work directory:
        original_dir = os.getcwd() # Save the current working directory
        megara_reduction_dir = os.path.join(original_dir, work_dir)
        os.chdir(megara_reduction_dir) 
        print('the directory has been changed to:', os.getcwd())

        # We have to change the rootdir for control.yaml:
        print('set the correct directory for work/MEGARA/control.yaml')
        command_control_yaml = """sed -i '' 's@rootdir:\ .*@rootdir:\ '"$PWD"'@' MEGARA/control.yaml"""
        print('\033[1m\033[31m ' + f"$ {command_control_yaml}" + '\033[0m\n')
        subprocess.run(command_control_yaml, shell=True)
        print('the rootdir has been set to:', os.getcwd())

        # Start the reduction process:
        reduction_dir = os.path.join(original_dir, work_megara_dir)
        os.chdir(reduction_dir)
        print('the directory has been changed to:', os.getcwd())
        
        reduce_simulations(i, nstart, abs_results_dir, run_modelmap, run_twilight, run_healing, run_LRU, history_line_command)
        
        # we go back to the directory where the script was executed:
        os.chdir(original_dir)   # Go back to original directory
        print('the directory has been changed to:', os.getcwd())




def main():
    parser = argparse.ArgumentParser(description='Simulate MEGARA reductions.')
    parser.add_argument('--obj_name', type=str, help='Name of the object to simulate.', default='obj_*')
    parser.add_argument('--vph', type=str, help='VPH name.', default='VPH_*')
    parser.add_argument('-n', '--num_simul', type=int, help='Number of simulations to perform.', default=1)
    parser.add_argument('--run_modelmap', action='store_true', help='Run ModelMap step.')
    parser.add_argument('--run_twilight', action='store_true', help='Run Twilight step.')
    args = parser.parse_args()

    obj_vph = f'{args.obj_name}/{args.vph}'
    ob_list = list(Path('.').glob(obj_vph))
    run_modelmap = args.run_modelmap
    run_twilight = args.run_twilight
    print(f"Running ModelMap: {run_modelmap}")
    print(f"Running Twilight: {run_twilight}")
    
    history_line_command = (
            f"$ python simulate_MEGARA_reductions.py "
            f"--obj_name {args.obj_name} --vph {args.vph} "
            f"--num_simul {args.num_simul} "
            f"{'--run_modelmap' if run_modelmap else ''} "
            f"{'--run_twilight' if run_twilight else ''}"
        )
    print(history_line_command)
    
    if args.num_simul < 1:
        raise ValueError("Number of simulations must be at least 1.")
    else:
        print(f"Number of simulations: {args.num_simul}")
        nsimul = args.num_simul

    if not ob_list:
        print(f"No galaxies found with name {obj_vph}")
        return
    for ob in ob_list:
        simulate_MEGARA_reductions(ob, nsimul, run_modelmap, run_twilight, history_line_command)
        

if __name__ == "__main__":
    main()
