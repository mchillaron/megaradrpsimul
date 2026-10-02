#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Convert the RSS file into a cube with the specified pixel size."""

import logging
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)

def conversion_into_cube(pixel_size, dest_cube_file, dest_file):

    """Convert the RSS file into a cube with the specified pixel size.
    Parameters
    ----------
    pixel_size : float
        Pixel size in arcseconds for the conversion of RSS into a cube.
    dest_cube_file : str
        Destination file path for the output cube.
    dest_file : str
        Path to the RSS file to be converted into a cube.
    """
    command_convert_rss_to_cube = [
        "megaradrp-cube",
        "-p",
        str(pixel_size),
        "-o",
        str(dest_cube_file),
        str(dest_file),
    ]

    logger.debug("[bold red]$ %s[/bold red]", " ".join(command_convert_rss_to_cube))

    t_start = datetime.now()
    sp = subprocess.run(command_convert_rss_to_cube, capture_output=True,text=True, check=True)
    t_stop = datetime.now()

    logger.debug("megaradrp-cube stdout:\n%s",sp.stdout)
    logger.debug("megaradrp-cube stderr:\n%s",sp.stderr)

    logger.info("RSS-to-cube conversion completed in %s", t_stop - t_start)