#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Get the step name and VPH name from the YAML file."""

import yaml
from pathlib import Path

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