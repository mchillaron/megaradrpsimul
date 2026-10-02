#
# Copyright 2025-2026 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Run one reduction step and optionally copy its product
   into the corresponding calibration tree."""

from datetime import datetime
from pathlib import Path

import logging
import shutil
import subprocess

logger = logging.getLogger(__name__)

def step_reduction(yaml_file, step_name, calib_label_dir,
                   product_file=None, 
                   calib_folder_path=None):
    
    """
    Run one reduction step and optionally copy its product
    into the corresponding calibration tree.

    Parameters
    ----------
    yaml_file : str
        Path to the YAML file.
    step_name : str
        The step name extracted from the YAML file.
    calib_label_dir : Path instance
        Path to the calibration label directory.
    product_file : str, optional
        The name of the product file to be copied (default is None).
    calib_folder_path : str, optional
        The path to the calibration folder (default is None).
    """

    command_run_list = [
        """numina""",
        """run""",
        yaml_file,
    ]
    
    logger.info("[bold red]$ %s[/bold red]", " ".join(command_run_list))

    t_start = datetime.now()
    sp = subprocess.run(command_run_list, capture_output=True, text=True, check=True) 
    t_stop = datetime.now()

    logger.debug("stdout:\n%s", sp.stdout)
    logger.debug("stderr:\n%s", sp.stderr)

    logger.info("Elapsed time: %s", t_stop - t_start)

    # Copy the product file to the calibration tree if specified
    if product_file is not None and calib_folder_path is not None:

        source = (Path(f"obsid{step_name}_results") / product_file)
        destination_dir = (calib_label_dir / calib_folder_path)

        if not source.is_file():
            raise FileNotFoundError(f"Reduction product not found: {source}")

        if not destination_dir.is_dir():
            raise FileNotFoundError(
                f"Calibration destination directory not found: "
                f"{destination_dir}"
            )
        
        destination = destination_dir / product_file
        shutil.copy2(source, destination)
        logger.info("[bold red]Copying %s -> %s[/bold red]", source, destination)