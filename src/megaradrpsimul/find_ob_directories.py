#
# Copyright 2025-2026 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Locate MEGARA observation directories starting from the current working directory."""

from pathlib import Path

import logging

logger = logging.getLogger(__name__)

def find_ob_directories(
    obj_pattern="obj_*",
    vph_pattern="VPH_*"
):
    """
    Locate MEGARA observation directories starting from the
    current working directory.

    The current directory can be:

    1. A root directory containing obj_* directories.
    2. An obj_* directory containing VPH_* directories.
    3. A VPH_* directory itself.

    Parameters
    ----------
    obj_pattern : str
        Pattern used to identify object directories.

    vph_pattern : str
        Pattern used to identify VPH directories.

    Returns
    -------
    list of Path
        Observation/VPH directories to process.
    """

    current_dir = Path.cwd()

    logger.debug("Current working directory: %s", current_dir)

    # Case 1: already inside a VPH_* directory
    if current_dir.name.startswith("VPH_"):
        logger.debug("Current directory identified as a VPH directory.")
        return [current_dir]

    # Case 2: inside an obj_* directory
    if current_dir.name.startswith("obj_"):
        logger.debug("Current directory identified as an object directory.")

        return sorted(
            directory
            for directory in current_dir.glob(vph_pattern)
            if directory.is_dir()
        )


    # Case 3: root directory containing obj_* directories
    logger.debug("Current directory identified as a root directory.")

    ob_list = []
    for obj_dir in sorted(current_dir.glob(obj_pattern)):

        if not obj_dir.is_dir():
            continue

        logger.debug("Object directory found: %s", obj_dir)

        ob_list.extend(
            sorted(
                directory
                for directory in obj_dir.glob(vph_pattern)
                if directory.is_dir()
            )
        )

    return ob_list