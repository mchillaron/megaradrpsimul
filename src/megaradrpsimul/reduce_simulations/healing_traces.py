#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Heal the traces for the given VPH and step name."""

import subprocess
from pathlib import Path

def healing_traces(vph_name, step_name):
    """Heal the traces for the given VPH and step name.
    Parameters
    ---------- 
    vph_name : str
        The name of the VPH (Variable Phase Hologram).
    step_name : str
        The name of the step for which traces are being healed.
    """
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