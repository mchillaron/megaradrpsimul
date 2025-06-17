import subprocess
import yaml
from astropy.io import fits
from pathlib import Path
from datetime import datetime
import os
import shutil

def get_step_name(yaml_file):
    """Get the step name and VPH name from the YAML file.
    Parameters
    ---------- 
    yaml_file : str
        Path to the YAML file.
    Returns
    -------
    step_name : str
        The step name extracted from the YAML file.
    vph_name : str
        The VPH name extracted from the YAML file.
    """

    with open(yaml_file, 'r') as f:
        data_yaml = yaml.safe_load(f)

    step_name = data_yaml.get('id', '')
    file_name = Path(yaml_file).name
    if not file_name.startswith('0'):
        vph_name = step_name.rsplit('_', 1)[-1] if '_' in step_name else ''
    else:
        vph_name = None

    return step_name, vph_name

def step_reduction(yaml_file, step_name, product_file=None, calib_folder_path=None):
    """Run each reduction step of the pipeline.
    Parameters
    ----------
    yaml_file : str
        Path to the YAML file.
    step_name : str
        The step name extracted from the YAML file.
    product_file : str, optional
        The name of the product file to be copied (default is None).
    calib_folder_path : str, optional
        The path to the calibration folder (default is None).
    """

    calib_dir_name = "ca3558e3-e50d-4bbc-86bd-da50a0998a48"

    command_run_list = [
        """numina""",
        """run""",
        f"""{yaml_file}""",
        """--link-files""",
        """-r""",
        """control.yaml"""
    ]
    print('\033[1m\033[31m ' + f"$ {' '.join(command_run_list)}" + '\033[0m\n')

    t_start = datetime.now()
    subprocess.run(command_run_list, capture_output=False, text=True)  # capture_output=True allows to capture stdout and stderr, also add sp = subprocess ...
    t_stop = datetime.now()
    #print(f'std_err: {sp.stderr}')
    #print(f'std_out: {sp.stdout}')
    print('\033[1m\033[32m ' + f'Elapsed time: {t_stop - t_start}' + '\033[0m\n')

    if product_file!= None and calib_folder_path != None:
        command_copy_list = [
            """cp""",
            f"""obsid{step_name}_results/{product_file}""",
            f"""{calib_dir_name}/{calib_folder_path}/"""
        ]

        print('\033[1m\033[31m ' + f"$ {' '.join(command_copy_list)}" + '\033[0m\n')
        sp = subprocess.run(command_copy_list, capture_output=True, text=True)
        #print(f'std_err: {sp.stderr}')
        #print(f'std_out: {sp.stdout}')

def healing_traces(vph_name, step_name):
    # First, we delete master_traces.json from calibration directory if it exists
    calib_dir_name = "ca3558e3-e50d-4bbc-86bd-da50a0998a48"
    calib_path = Path(calib_dir_name) / f"TraceMap/LCB/{vph_name}/master_traces.json"
    if calib_path.exists():
        print(f"Deleting existing {calib_path}")
        calib_path.unlink()
    else:
        print(f"No existing {calib_path} to delete.")
    
    # Then, we run the healing command
    command_healing_list = [
        """megaradrp-heal_traces""",
        f"""obsid{step_name}_results/reduced_image.fits""",
        f"""obsid{step_name}_results/master_traces.json""",
        """--fibids""",
        """--healing""",
        """healing.yaml""",
        """--updated_traces""",
        f"""obsid{step_name}_results/master_traces_healed.json""",
        """--pdffile""",
        f"""obsid{step_name}_results/healed_traces.pdf"""
    ]
    print('\033[1m\033[35m ' + f"$ {' '.join(command_healing_list)}" + '\033[0m\n')
    subprocess.run(command_healing_list, capture_output=True, text=True)

    command_copy_healing_list = [
        """cp""",
        f"""obsid{step_name}_results/master_traces_healed.json""",
        f"""{calib_dir_name}/TraceMap/LCB/{vph_name}/"""
        ]
    print('\033[1m\033[35m ' + f"$ {' '.join(command_copy_healing_list)}" + '\033[0m\n')
    subprocess.run(command_copy_healing_list, capture_output=True, text=True)

def reduce_simulations(niter, nstart, abs_results_dir, run_modelmap=False, run_twilight=False, run_healing=False, run_LRU=False, history_line_command=None):
    """Reduce the simulated images using the TEA pipeline.

    Parameters
    ----------
    niter : int
        Number of the current simulation.
    nstart : int
        Number from which the simulation process starts.
    abs_results_dir : Path instance
        Path to the absolute results directory.
    run_modelmap : bool, optional
        If True, run the ModelMap step (default is False).
    run_twilight : bool, optional
        If True, run the TwilightFlat step (default is False).
    
    """
    
    print('we start the reduction process')

    print('........... Step 0: Bias image ...........')
    bias_yaml_file = list(Path('.').glob('0_*.yaml'))[0].name   # .name obtains the string with the name from the Path instance
    step_name, vph_name = get_step_name(bias_yaml_file)
    print('Step name:', step_name)
    step_reduction(bias_yaml_file, step_name, product_file="master_bias.fits", calib_folder_path="MasterBias")
    
    print('........... Step 1: TraceMap ...........')
    tracemap_yaml_file = list(Path('.').glob('1_*.yaml'))[0].name
    step_name, vph_name = get_step_name(tracemap_yaml_file)
    print('Step name:', step_name)
    print('VPH name:', vph_name)

    if run_LRU == True:
        print('Running special TraceMap for LR-U')
        command_copy_traces_LRU_list = [
            """cp""",
            """master_traces_LRU_20220325_healed.json""",
            f"""ca3558e3-e50d-4bbc-86bd-da50a0998a48/TraceMap/LCB/{vph_name}/"""
            ]
        print('\033[1m\033[35m ' + f"$ {' '.join(command_copy_traces_LRU_list)}" + '\033[0m\n')
        subprocess.run(command_copy_traces_LRU_list, capture_output=True, text=True)
    else:
        step_reduction(tracemap_yaml_file, step_name, product_file="master_traces.json", calib_folder_path=f"TraceMap/LCB/{vph_name}")
        if run_healing == True:
                print("Healing of the traces required")
                healing_traces(vph_name, step_name)
    
        

    if run_modelmap == True:
        print('........... Step 2: ModelMap ...........')
        modelmap_yaml_file = list(Path('.').glob('2_*.yaml'))[0].name
        step_name = get_step_name(modelmap_yaml_file)
        step_reduction(modelmap_yaml_file, step_name, product_file="master_model.json", calib_folder_path=f"ModelMap/LCB/{vph_name}")

    print('........... Step 3: WavelengthCalibration ...........')
    wavecalib_yaml_file = list(Path('.').glob('3_*b.yaml'))[0].name
    step_name, vph_name = get_step_name(wavecalib_yaml_file)
    print('Step name:', step_name)
    step_reduction(wavecalib_yaml_file, step_name, product_file="master_wlcalib.json", calib_folder_path=f"WavelengthCalibration/LCB/{vph_name}")

    print('........... Step 3: WavelengthCalibration - Check ...........')
    wavecalibcheck_yaml_file = list(Path('.').glob('3_*check.yaml'))[0].name
    step_name, vph_name = get_step_name(wavecalibcheck_yaml_file)
    step_reduction(wavecalibcheck_yaml_file, step_name)

    print('........... Step 4: FiberFlat ...........')
    fiberflat_yaml_file = list(Path('.').glob('4_*.yaml'))[0].name
    step_name, vph_name = get_step_name(fiberflat_yaml_file)
    print('Step name:', step_name)
    step_reduction(fiberflat_yaml_file, step_name, product_file="master_fiberflat.fits", calib_folder_path=f"MasterFiberFlat/LCB/{vph_name}")


    if run_twilight == True:
        print('........... Step 5: Bias image ...........')
        twilight_yaml_file = list(Path('.').glob('5_*.yaml'))[0].name
        step_name = get_step_name(twilight_yaml_file)
        print('Step name:', step_name)
        step_reduction(twilight_yaml_file, step_name, product_file="master_twilightflat.fits", calib_folder_path=f"MasterTwilightFlat/LCB/{vph_name}")
        
    print('........... Step 6: LCBadquisition ...........')
    lcbadquisition_yaml_file = list(Path('.').glob('6_*.yaml'))[0].name
    step_name, vph_name = get_step_name(lcbadquisition_yaml_file)
    print('Step name:', step_name)
    step_reduction(lcbadquisition_yaml_file, step_name)

    print('........... Step 7: StandardStar ...........')
    standardstar_yaml_file = list(Path('.').glob('7_*.yaml'))[0].name
    step_name, vph_name = get_step_name(standardstar_yaml_file)
    print('Step name:', step_name)
    step_reduction(standardstar_yaml_file, step_name, product_file="master_sensitivity.fits", calib_folder_path=f"MasterSensitivity/LCB/{vph_name}")

    print('........... Step 8: Reduce LCB ...........')
    reduce_yaml_file = list(Path('.').glob('8_*.yaml'))[0].name
    step_name, vph_name = get_step_name(reduce_yaml_file)
    print('Step name:', step_name)
    step_reduction(reduce_yaml_file, step_name)

    print('End of the reduction process')
    
    num = niter + nstart
    # Change the name of the final_rss.fits file to the name of the simulation
    original_file = f"obsid{step_name}_results/final_rss.fits"
    new_file = f"obsid{step_name}_results/final_rss_{num:04d}.fits"

    if os.path.exists(original_file):
        os.rename(original_file, new_file)
        with fits.open(new_file, mode='update') as hdul:
            hdul[0].header.add_history(history_line_command)
            hdul.flush() 
        print('Command line added to final_rss.fits HISTORY keyword')
    else:
        print(f"El archivo {original_file} no existe.")


    dest_file = os.path.join(abs_results_dir, os.path.basename(new_file))
    # If there is already a file with the same name in the results directory, we overwrite it
    if os.path.exists(dest_file):
        os.remove(dest_file)  # Elimina el archivo para poder mover el nuevo
        print(f"Existing {dest_file} being deleted.")

    shutil.move(new_file, abs_results_dir)
    print(f"Moved {new_file} to {abs_results_dir}")


    print('the simulation and reduction process has finished')
    print('................................................................................')