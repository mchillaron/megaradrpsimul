import yaml

def open_read_yaml(directory_MEGARA, yaml_file_name):
    # Open the YAML file
    with open(directory_MEGARA/yaml_file_name, 'r') as file:
        # Load the content of the file into a dictionary
        yaml_content = yaml.safe_load(file)

    # Now `data` is a dictionary containing the YAML content
    # We obtain the images list:
    list_frames_names = yaml_content['frames']
    
    # Now we generate the complete rute to the file:
    complete_list = [directory_MEGARA/'data'/filename for filename in list_frames_names]
    return(complete_list)
    