#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Convert the RSS file into a cube with the specified pixel size."""

import subprocess
from datetime import datetime

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
        """megaradrp-cube""",
        """-p""", 
        f"""{pixel_size}""",
        """-o""", 
        f"""{dest_cube_file}""",
        f"""{dest_file}"""]
    
    print('\033[1m\033[31m ' + f"$ {' '.join(command_convert_rss_to_cube)}" + '\033[0m\n')
    t_start = datetime.now()
    subprocess.run(command_convert_rss_to_cube, capture_output=False, text=True)  # capture_output=True allows to capture stdout and stderr, also add sp = subprocess ...
    t_stop = datetime.now()
    print('\033[1m\033[32m ' + f'Elapsed time: {t_stop - t_start}' + '\033[0m\n')