#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Determine the diffuse light for the given VPH and step name."""

import subprocess

def diffuselight_determination(vph_name, step_name_8, run_LRU, run_healing):
    """Determine the diffuse light for the given VPH and step name. 
    Parameters
    ----------
    vph_name : str
        The name of the VPH (Variable Phase Hologram).
    step_name_8 : str
        The name of the step for which diffuse light is being determined (Step number 8)
    run_LRU : bool
        If True, run the special TraceMap template for LR-U (default is False).
    run_healing : bool
        If True, run the healing of the traces (default is False)."""
    
    # Need to know the correct traces file for this step
    if run_LRU == True:
        traces_file_for_this_step = f"ca3558e3-e50d-4bbc-86bd-da50a0998a48/TraceMap/LCB/{vph_name}/master_traces_LRU_20220325_healed.json"
    if run_healing == True:
        traces_file_for_this_step = f"ca3558e3-e50d-4bbc-86bd-da50a0998a48/TraceMap/LCB/{vph_name}/master_traces_healed.json"
    else:
        traces_file_for_this_step = f"ca3558e3-e50d-4bbc-86bd-da50a0998a48/TraceMap/LCB/{vph_name}/master_traces.json"

    # Firstly, we run the diffuse light determination step
    command_diffuselight_list = [
        """megaratools-diffuse_light""",
        """-i""",
        f"""obsid{step_name_8}_work/reduced_image.fits""",
        """-o""",
        """data/background_2D.fits""",
        """-r""", 
        """data/residuals_2D.fits""",
        """-t""",
        f"""{traces_file_for_this_step}""",
        """-s""",
        """-3""",
        """-p""",
        """data/plots_2D.pdf""",
        """-2D"""
    ]

    print('\033[1m\033[35m ' + f"$ {' '.join(command_diffuselight_list)}" + '\033[0m\n')
    subprocess.run(command_diffuselight_list, capture_output=True, text=True)