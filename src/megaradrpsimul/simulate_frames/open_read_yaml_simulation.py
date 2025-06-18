#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of megaradrpsimul.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Open and read a YAML file, returning the list of images."""

import yaml

def open_read_yaml_simulation(megara_dir, yaml_filename): 
    """Open and read a YAML file, returning the list of images.
    
    Parameters
    ----------
    megara_dir : Path instance
        Path to the MEGARA directory.
    yaml_filename : Path instance
        Path to the yaml file from which the name of the frames are obtained.
        
    Returns
    -------
    list : list of Path instances
        List of image paths."""
    
    with open(yaml_filename, 'r') as file:
        yaml_content = yaml.safe_load(file)
    
    list_frames_names = yaml_content['frames']
    # Now we generate the complete rute to the file:
    complete_list = [megara_dir/'data'/filename for filename in list_frames_names]
    return(complete_list)