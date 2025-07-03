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

def get_step_name(yaml_file, extraction_offset=False):
    """Get the step name and VPH name from the YAML file.
    Parameters
    ---------- 
    yaml_file : str
        Name of the YAML file.
    extraction_offset : bool, optional
        If True, also extract the extraction offset from the YAML file (default is False).
    Returns
    -------
    step_name : str
        The step name extracted from the YAML file.
    offset : int or None
        The extraction offset if it exists, otherwise None.
    """

    with open(yaml_file, 'r') as f:
        data_yaml = yaml.safe_load(f)

    step_name = data_yaml.get('id', '')

    if extraction_offset:
        try:
            offset = data_yaml['requirements']['extraction_offset'][0]
        except (KeyError, IndexError, TypeError):
            print("extraction_offset not found")
            offset = None
        return step_name, offset

    return step_name