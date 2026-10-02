#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Determine the diffuse light for the given VPH and step name."""

import logging
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)
    
def diffuselight_determination(vph_name, step_name_8, calib_label_dir,
                               run_LRU, run_healing, extraction_offset):
    
    """Determine the diffuse light for the given VPH and step 8.

    Parameters
    ----------
    vph_name : str
        The name of the VPH (Variable Phase Hologram).
    step_name_8 : str
        Step name corresponding to the standard Step 8 reduction.
    calib_label_dir : Path
        Calibration directory corresponding to the current INSCONF label.
    run_LRU : bool
        If True, run the special TraceMap template for LR-U (default is False).
    run_healing : bool
        If True, run the healing of the traces (default is False)."""
    
    # Need to know the correct traces file for this step
    traces_dir = (calib_label_dir/ "TraceMap"/ "LCB"/ vph_name)

    if run_LRU:
        traces_file_for_this_step = (traces_dir / "master_traces_LRU_20220325_healed.json")
    elif run_healing:
        traces_file_for_this_step = (traces_dir / "master_traces_healed.json")
    else:
        traces_file_for_this_step = (traces_dir / "master_traces.json")
    
    if not traces_file_for_this_step.is_file():
        raise FileNotFoundError(f"Trace file required for diffuse-light correction was not found: {traces_file_for_this_step}")
    
    logger.debug("Trace file used for diffuse-light correction: %s", traces_file_for_this_step)

    # Firstly, we run the diffuse light determination step
    # Run diffuse-light determination
    command_diffuselight_list = [
        "megaratools-diffuse_light",
        "-i",
        f"obsid{step_name_8}_work/reduced_image.fits",
        "-o",
        "data/background_2D.fits",
        "-r",
        "data/residuals_2D.fits",
        "-t",
        str(traces_file_for_this_step),
        "-s",
        str(extraction_offset),
        "-p",
        "data/plots_2D.pdf",
        "-2D",
    ]

    logger.debug("[bold magenta]$ %s[/bold magenta]", " ".join(command_diffuselight_list))

    t_start = datetime.now()
    sp = subprocess.run(command_diffuselight_list, capture_output=True, text=True, check=True)
    t_stop = datetime.now()

    logger.debug("Diffuse-light stdout:\n%s", sp.stdout)
    logger.debug("Diffuse-light stderr:\n%s", sp.stderr)

    logger.info("Diffuse-light determination completed in %s", t_stop - t_start)