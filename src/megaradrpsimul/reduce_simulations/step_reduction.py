#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Run each reduction step of the pipeline."""

import subprocess
from datetime import datetime


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