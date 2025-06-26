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

def get_step_name(yaml_file):
    """Get the step name and VPH name from the YAML file.
    Parameters
    ---------- 
    yaml_file : str
        Name of the YAML file.
    Returns
    -------
    step_name : str
        The step name extracted from the YAML file.
    """

    with open(yaml_file, 'r') as f:
        data_yaml = yaml.safe_load(f)

    step_name = data_yaml.get('id', '')

    return step_name