from astropy import units as u
from astropy.io import fits
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter
from tqdm import tqdm

from megaradrpsimul import CCDregions
import numpy as np
import pickle
import shutil
import teareduce as tea
import yaml

# --------------------------------------------------------
def open_read_yaml(megara_dir, yaml_filename): 
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

# --------------------------------------------------------
def cosmicray_cleaning(list_bias):
    """Clean the images from cosmic rays using the regions defined in CCDregions.py.
    
    This function uses the tea.cr2images function to clean the images, and the
    regions are defined in the CCDregions.py file. The function returns a dictionary 
    with the cleaned images, whose keys are the names of the images.

    Parameters
    ----------
    list_bias : list of Path instances
        List of image paths to be cleaned.
        
    Returns
    -------
    cleaned_images : dict
        Dictionary with cleaned images."""

    print('cleaning cosmic rays')
    regions_cosmicrays = CCDregions.regions_cosmicrays
    cleaned_images = {}

    for idx, img_path in tqdm(enumerate(list_bias), total=len(list_bias), desc="Image processing", unit="img"):
        with fits.open(img_path) as hdul:
            data = hdul[0].data.astype(float)  # we convert to float to avoid errors
            naxis2, naxis1 = data.shape        # dimensions of the image

            print(f"{idx}.Processing image {img_path} with shape ({naxis1}, {naxis2})")

            cleaned_data = data.copy()         # creating a copy of the image to do changes

            for region_name, region_params in regions_cosmicrays.items():
                slice2d = region_params['slice2d']
                print(f"Processing region {slice2d}")
                median_size = region_params['median_size']
                tsigma_peak = region_params['tsigma_peak']
                tsigma_tail = region_params['tsigma_tail']

                # cleaning the region defined by slice2d using tea.cr2images
                cleaned_region = tea.cr2images(
                    data1=data, 
                    median_size=median_size,
                    tsigma_peak=tsigma_peak,
                    tsigma_tail=tsigma_tail,
                    image_region = slice2d,
                    return_masks = False,
                    debug_level = 0
                )
                cleaned_data[slice2d.python] = cleaned_region[0][slice2d.python]
            
            # We save the cleaned image in a dictionary:
            cleaned_image_name = f"cleaned_bias_{idx}"
            cleaned_images[cleaned_image_name] = cleaned_data
            print(f"Image {img_path} cleaned from cosmis rays and saved as ({cleaned_image_name})")
            print("................................................................................")
    return cleaned_images

# --------------------------------------------------------
def smooth_frames(cleaned_images):
    """Smooth the images using Savitzky-Golay and median filters.
    
    Parameters
    ----------
    cleaned_images : dict
        Dictionary with cleaned images to be smoothed.
    Returns
    -------
    smoothed_images : dict
        Dictionary with smoothed images."""
    

    regions_kernel = CCDregions.regions_kernel
    smoothed_images = {}

    for idx, (name, image) in enumerate(cleaned_images.items()):
        data = image
        naxis2, naxis1 = data.shape

        # create an image of the same dimensions full of zeros
        zero_raw_image = np.zeros((naxis2, naxis1))

        for i, (region_name, region_params) in enumerate(regions_kernel.items()): 
            slice2d = region_params['slice2d']
            num_filters_SG = region_params['num_filters_SG']
            
            region_data = data[slice2d.python]
            
            if num_filters_SG == 0:
                if 'median_size' in region_params:
                    median_size = region_params['median_size']
                    filtered_data = median_filter(region_data, size=median_size)  # vertical median filter in these regions
                    print(f"Median filter applied to {region_name}")
                else:
                    raise ValueError(f"'median_size' is missing in region {region_name}")
                
            elif num_filters_SG != 0:
                size1_SG = region_params['size1_SG']
                pol_order1_SG = region_params['pol_order1_SG']
                axis1_SG = region_params['axis1_SG']

                if size1_SG % 2 == 0:
                    size1_SG += 1    # make sure it is odd
                
                filtered_data = savgol_filter(region_data, window_length=size1_SG, polyorder=pol_order1_SG, axis=axis1_SG)   # Aplica solo un filtro de Savitzky-Golay
                print(f"Savitzky-Golay applied in {region_name}")
        
                if num_filters_SG == 2:
                    size2_SG = region_params['size2_SG']
                    pol_order2_SG = region_params['pol_order2_SG']
                    axis2_SG = region_params['axis2_SG']

                    if size2_SG % 2 == 0:
                        size2_SG += 1 

                    filtered_data = savgol_filter(filtered_data, window_length=size2_SG, 
                                                  polyorder=pol_order2_SG, axis=axis2_SG)
                    print(f"Second Savitzky-Golay filter applied in {region_name}")

            zero_raw_image[slice2d.python] = filtered_data
            
        # Save the image in the dictionary "smoothed_{idx}"
        smoothed_image_name = f"smoothed_{idx}"
        smoothed_images[smoothed_image_name] = zero_raw_image

        if np.any(zero_raw_image == 0):
            print("There are still zeros in the smoothed bias image.")
        else:
            print("No zeros in the smoothed bias image.")
            
    return smoothed_images

# --------------------------------------------------------
def simulation_run_save(data_work_dir, img_path, generator, type_image): 
    """Simulate using the generator and save images.

    Parameters
    ----------
    data_work_dir : Path instance
        Path to the data work directory.
    img_path : Path instance
        Path to the image to be simulated.
    type_image : str
        Type of image ('bias', 'object').
    generator : instance of SimulateCCDExposure
        Generator object for simulating images.
    """
    destination_dir = data_work_dir
    image_exposure = generator.run(imgtype=type_image)

    new_filename = img_path.name  
    destination_file = destination_dir / new_filename  # save the simulated image with the same name as the original

    # Copying the original image to the destination directory
    # then we only update the data with the simulated data
    shutil.copy(img_path, destination_file)
    print(f"Image {img_path} copied at: {destination_file}")
    
    with fits.open(destination_file, mode='update') as hdul:
        data_round = np.round(image_exposure.data).astype(np.uint16)  
        hdul[0].data = data_round
        hdul.flush()
        print(f"Updated data for: {destination_file}")
    
# --------------------------------------------------------
def calculate_bias_readout_noise(data, data_smoothed, naxis2, naxis1):
    """
    Calculate the readout noise across the CCD based on top and bottom bias regions.

    The function computes the pixel-wise difference between the raw and smoothed
    data, then evaluates the robust standard deviation in predefined top and bottom
    CCD bias regions to estimate readout noise. The resulting noise is propagated
    into a full-sized array, assigning different values to top and bottom sections.

    Parameters
    ----------
    data : ndarray
        The original 2D CCD image data.
    data_smoothed : ndarray
        The smoothed version of the CCD image, used as a reference for bias subtraction.
    naxis2 : int
        Number of rows (Y-axis size) of the CCD image.
    naxis1 : int
        Number of columns (X-axis size) of the CCD image.

    Returns
    -------
    noise_array : ndarray
        A 2D array of the same shape as the input, filled with the estimated
        readout noise values, using different estimates for top and bottom regions.
    """
    data_bias_diff = data - data_smoothed

    regions_kernel = CCDregions.regions_kernel
    # Top CCD region
    top_slice = regions_kernel["bias_topCCD"]["slice2d"].python
    top_data = data_bias_diff[top_slice]
    noise_top = tea.robust_std(top_data)
    # Bottom CCD region
    bottom_slice = regions_kernel["bias_bottomCCD"]["slice2d"].python
    bottom_data = data_bias_diff[bottom_slice]
    noise_bottom = tea.robust_std(bottom_data)

    #Create the readout_noise array
    noise_array = np.full((naxis2, naxis1), noise_top)
    noise_array[0:2106, :] = noise_bottom

    print(f"Robust deviation for TopCCD is: {noise_top} and for BottomCCD is: {noise_bottom}")

    return noise_array

# --------------------------------------------------------
def calculate_image_noise(data_model_traces, gain_array, readout_noise_median):
    """
    Calculate the total noise in an image considering photon shot noise and readout noise.

    This function estimates the pixel-wise noise in units of ADU (Analog-to-Digital Units),
    combining the contribution from photon statistics (shot noise) and the detector's 
    readout noise. The photon noise is computed as the square root of the signal divided 
    by the gain, and the readout noise is included in quadrature.

    Parameters
    ----------
    data_model_traces : ndarray
        2D array (or broadcastable shape) of the modelled signal in ADU.
    gain_array : ndarray or float
        2D array or scalar representing the gain (electrons per ADU) for each pixel.
        Must be broadcast-compatible with `data_model_traces`.
    readout_noise_median : float or ndarray
        Scalar or array with the median readout noise (in ADU). Can be a single value
        or an array matching the shape of the image.

    Returns
    -------
    noise_array : ndarray
        Array of the same shape as `data_model_traces` containing the estimated noise per pixel,
        in ADU.

    Notes
    -----
    The total noise is computed using:

        noise = sqrt( (signal / gain) + readout_noise**2 )

    where:
    - signal / gain → photon shot noise (in electrons)
    - readout_noise**2 → detector noise (in ADU^2), added in quadrature
    """
    ratio_data_gain = data_model_traces / gain_array
    readout_noise_square = readout_noise_median ** 2
    noise_array = np.sqrt(ratio_data_gain + readout_noise_square)
    return noise_array

# --------------------------------------------------------
def calculate_median_offsets(list_images, bias_smoothed_median):
    """
    Calculate the median image and offset values for a set of CCD images.

    This function computes the median image from a list of input CCD images,
    and compares overscan regions (top and bottom) against a smoothed bias
    median to evaluate offset consistency. It returns the image median,
    individual image offsets, and the difference between the median image
    and the reference smoothed bias.

    Parameters
    ----------
    list_images : list of Path
        List of file paths to FITS images used for median calculation.
    bias_smoothed_median : ndarray
        A reference image representing the smoothed median bias, used to
        calculate absolute offsets.

    Returns
    -------
    median_image : ndarray
        The median image computed from all images in `list_images`.
    image_median_offsets : list of tuple of float
        List of (top_offset, bottom_offset) tuples representing how much
        each image deviates from the median image in the overscan regions.
    top_bias_median_offset : float
        Offset between the top overscan region of the median image and the bias.
    bottom_bias_median_offset : float
        Offset between the bottom overscan region of the median image and the bias.
    """
    top_bias_constant = np.median(bias_smoothed_median[CCDregions.overscan3_top.python])
    bottom_bias_constant = np.median(bias_smoothed_median[CCDregions.overscan3_bottom.python])

    image_data = []
    image_offsets = []

    for path in list_images:
        with fits.open(path) as hdul:
            data = hdul[0].data
            image_data.append(data)

            top_offset = np.median(data[CCDregions.overscan3_top.python])
            bottom_offset = np.median(data[CCDregions.overscan3_bottom.python])
            image_offsets.append((top_offset, bottom_offset))
            print(f"Individual image || Top constant: {top_offset}, Bottom constant: {bottom_offset}")

    stacked = np.stack(image_data, axis=0)
    median_image = np.median(stacked, axis=0)
    top_median_constant = np.median(median_image[CCDregions.overscan3_top.python])
    bottom_median_constant = np.median(median_image[CCDregions.overscan3_bottom.python])

    top_bias_median_offset = top_bias_constant - top_median_constant
    bottom_bias_median_offset = bottom_bias_constant - bottom_median_constant

    image_median_single_offsets = [
        (top_median_constant - top, bottom_median_constant - bottom)
        for top, bottom in image_offsets
    ]

    print(f"MasterBias || Top constant: {top_bias_constant}, Bottom constant: {bottom_bias_constant}")
    print(f"Median image || Top constant: {top_median_constant}, Bottom constant: {bottom_median_constant}")
    print("Offsets of Median image - Individual image:")
    for i, offsets in enumerate(image_median_single_offsets):
        print(f"Image {i}: Top offset = {offsets[0]}, Bottom offset = {offsets[1]}")
    print("Offset of MasterBias - Median image:")
    print(f"Top offset = {top_bias_median_offset}, Bottom offset = {bottom_bias_median_offset}")

    return median_image, image_median_single_offsets, top_bias_median_offset, bottom_bias_median_offset

# --------------------------------------------------------
def simulate_frames_step(list_images, 
                        naxis1, naxis2,
                        data_work_dir,
                        name_step, 
                        bias_smoothed_median, 
                        gain_array, 
                        readout_noise_median):
        
        """Simulate frames for a specific step (e.g., TraceMap, ArcCalibration).

        This function processes a list of images, calculates the median image,
        and generates corrected data models for each image. It also calculates
        the noise array based on the data model, gain, and readout noise.
        Parameters
        ----------
        list_images : list of Path instances
            List of image paths to be processed.
        naxis1 : int
            Number of columns (X-axis size) of the CCD image.
        naxis2 : int
            Number of rows (Y-axis size) of the CCD image.
        data_work_dir : Path instance
            Path to the data work directory where the simulated images will be saved.
        name_step : str
            Name of the step for which the frames are being simulated (e.g., 'TraceMap').
        bias_smoothed_median : ndarray
            A reference image representing the smoothed median bias, used to
            calculate absolute offsets.
        gain_array : ndarray
            2D array representing the gain (electrons per ADU) for each pixel.
        readout_noise_median : ndarray
            2D array with the median readout noise (in ADU) for each pixel. 

        """
        median_image, image_single_median_offsets, top_bias_median_offset, bottom_bias_median_offset = calculate_median_offsets(list_images, bias_smoothed_median)
        print(f'median image for {name_step} frames calculated')

        data_model = median_image - bias_smoothed_median
        data_model_corrected = data_model.copy()
        data_model_corrected[CCDregions.topCCD_full.python] += top_bias_median_offset 
        data_model_corrected[CCDregions.bottomCCD_full.python] += bottom_bias_median_offset
        data_model_corrected[data_model_corrected < 0] = 0
        print(f'Corrected {name_step} data model obtained for generator')
        
        noise_array = calculate_image_noise(data_model_corrected, gain_array, readout_noise_median)
        
        generators = {}
        for idx, (image_path, (top_offset, bottom_offset)) in enumerate(zip(list_images, image_single_median_offsets)):
            # create a copy of the bias smoothed median and apply the correction of the offsets
            bias_smoothed_median_corrected = np.array(bias_smoothed_median, copy=True)
            correction_top = np.uint16(top_bias_median_offset + top_offset)
            correction_bottom = np.uint16(bottom_bias_median_offset + bottom_offset)
            bias_smoothed_median_corrected[CCDregions.topCCD_full.python] -= correction_top            
            bias_smoothed_median_corrected[CCDregions.bottomCCD_full.python] -= correction_bottom 

            # now I calculate the median of the corrected bias just to print it and compare with the original values of the images.
            top_bias_corrected_constant = np.median(bias_smoothed_median_corrected[CCDregions.regions_kernel["bias_topCCD"]["slice2d"].python])
            bottom_bias_corrected_constant = np.median(bias_smoothed_median_corrected[CCDregions.regions_kernel["bias_bottomCCD"]["slice2d"].python])
            print(f'The values of the corrected MasterBias are: {top_bias_corrected_constant} and {bottom_bias_corrected_constant}')

            # create the new generator
            image_generator = tea.SimulateCCDExposure(
                naxis1=naxis1, naxis2=naxis2, bitpix=16,
                bias=bias_smoothed_median_corrected * u.adu,
                gain=gain_array * (u.electron/u.adu),
                readout_noise=noise_array * u.adu,
                flatfield=1,
                dark=0 * u.adu,
                data_model=data_model_corrected * u.adu,
            )

            # Save the generator as "generator_{idx}"
            generator_name = f"generator_{idx}"
            generators[generator_name] = image_generator

            print(f'{name_step} generator created')
    
        # run the simulation for each image
        for idx, img_path in tqdm(enumerate(list_images), total=len(list_images), desc=f"Simulating {name_step} images", unit="img"):
            generator_name = f"generator_{idx}"
            generator = generators[generator_name]
            simulation_run_save(data_work_dir, img_path, generator, type_image='object')

# --------------------------------------------------------
def simulate_frames(megara_dir, data_work_dir, work_megara_dir):
    """Simulate frames for MEGARA data reduction.
    This function simulates the frames for the MEGARA data reduction process.
    It processes bias, trace map, arc calibration, LcbImage, and science object frames.
    
    Parameters
    ----------
    megara_dir : Path instance
        Path to the MEGARA directory.
    work_dir : Path instance
        Path to the work directory.
    work_megara_dir : Path instance
        Path to the work MEGARA directory.
   
    """
    
    print('starting simulation of frames')

    #.................................... BIAS ..........................................
    yaml_file_bias = list(work_megara_dir.glob('0_*.yaml'))[0]   # it is a path instance
    list_bias = open_read_yaml(megara_dir, yaml_file_bias)
    print('creating bias frames list')
    
    bias_cleaned_images = cosmicray_cleaning(list_bias)
    bias_smoothed_images = smooth_frames(bias_cleaned_images)

    with open(megara_dir/'smoothed_bias.pkl', 'wb') as f:
        pickle.dump(bias_smoothed_images, f)

    #with open(megara_dir/'smoothed_bias.pkl', 'rb') as f:
    #    bias_smoothed_images = pickle.load(f)

    # Median of all the smoothed bias images:
    list_smoothed_bias = list(bias_smoothed_images.values())
    stack_smoothed_images = np.stack(list_smoothed_bias, axis=0)      
    bias_smoothed_median = np.median(stack_smoothed_images, axis=0)

    # we convert to integer and round the values to avoid problems with the data type
    bias_smoothed_median_round = np.round(bias_smoothed_median).astype(np.uint16)

    # Bias generators
    bias_generators = {}
    all_noise_arrays = []

    for idx, img_path in tqdm(enumerate(list_bias), total=len(list_bias), desc="Image processing", unit="img"):
        with fits.open(img_path) as hdul:
            data = hdul[0].data  
            naxis2, naxis1 = data.shape
        
        # Calculate the noise for each of the CCDs
        noise_array = calculate_bias_readout_noise(data, bias_smoothed_images[f"smoothed_{idx}"], naxis2, naxis1)
        all_noise_arrays.append(noise_array)

        bias_image_generator = tea.SimulateCCDExposure(
            naxis1=naxis1, naxis2=naxis2, bitpix=16,
            bias=bias_smoothed_images[f"smoothed_{idx}"] * u.adu,
            readout_noise = noise_array * u.adu 
        )
        
        # Save the generator as "generator_{idx}"
        generator_name = f"generator_{idx}"
        bias_generators[generator_name] = bias_image_generator

    print('Bias generators created')

    all_noise_arrays = np.array(all_noise_arrays)
    #readout_noise_mean = np.mean(all_noise_arrays, axis=0) 
    readout_noise_median = np.median(all_noise_arrays, axis=0)

    # Simulation of bias images:
    for idx, img_path in tqdm(enumerate(list_bias), total=len(list_bias), desc="Simulating bias images", unit="img"):
        generator_name = f"generator_{idx}"
        generator = bias_generators[generator_name]
        simulation_run_save(data_work_dir, img_path, generator, type_image='bias')

    print('\033[1m\033[34m ' + "all bias images have been simulated" + '\033[0m\n')

    #.................................... Gain ..........................................
    # Gain values from header keywords
    gain_top = 1.6 * (u.electron/u.adu) 
    gain_bottom = 1.73 * (u.electron/u.adu)
    gain_array = np.full((naxis2, naxis1), gain_top.value)
    gain_array[0:2106, :] = gain_bottom.value

    #.................................... TRACEMAP ..........................................
    yaml_file_traces = list(work_megara_dir.glob('1_*.yaml'))[0]
    list_traces = open_read_yaml(megara_dir, yaml_file_traces)
    simulate_frames_step(list_images=list_traces, 
                        naxis1=naxis1,naxis2=naxis2,
                        data_work_dir=data_work_dir,
                        name_step='TraceMap',
                        bias_smoothed_median=bias_smoothed_median_round,
                        gain_array=gain_array,
                        readout_noise_median=readout_noise_median) 
    
    print('\033[1m\033[34m ' + "all TraceMap images have been simulated" + '\033[0m\n')


    #.................................... ArcCalibration ....................................
    yaml_file_arc = list(work_megara_dir.glob('3_*b.yaml'))[0]  # adding the b helps us distinguish between 3_WaveCalib.yaml and 3_WaveCalib_check.yaml
    list_arc = open_read_yaml(megara_dir, yaml_file_arc)
    
    simulate_frames_step(list_images=list_arc,
                        naxis1=naxis1,naxis2=naxis2,
                        data_work_dir=data_work_dir,
                        name_step='ArcCalibration',
                        bias_smoothed_median=bias_smoothed_median_round,
                        gain_array=gain_array,
                        readout_noise_median=readout_noise_median)
    
    print('\033[1m\033[34m ' + "all ArcCalibration images have been simulated" + '\033[0m\n')
    
    # ..................................LcbImage Star ..........................................
    yaml_file_lcb = list(work_megara_dir.glob('6_*.yaml'))[0]   # it is a path instance
    list_lcb = open_read_yaml(megara_dir, yaml_file_lcb)
    
    simulate_frames_step(list_images=list_lcb,
                        naxis1=naxis1,naxis2=naxis2,
                        data_work_dir=data_work_dir,
                        name_step='LcbImage',
                        bias_smoothed_median=bias_smoothed_median_round,
                        gain_array=gain_array,
                        readout_noise_median=readout_noise_median)
    
    print('\033[1m\033[34m ' + "all LcbImage images have been simulated" + '\033[0m\n')

    # ..................................Science Object ..........................................
    yaml_file_lcb_object = list(work_megara_dir.glob('8_*.yaml'))[0]   # it is a path instance
    list_lcb_object = open_read_yaml(megara_dir, yaml_file_lcb_object)
    simulate_frames_step(list_images=list_lcb_object,
                        naxis1=naxis1,naxis2=naxis2,
                        data_work_dir=data_work_dir,
                        name_step='Science LcbImage',
                        bias_smoothed_median=bias_smoothed_median_round,
                        gain_array=gain_array,
                        readout_noise_median=readout_noise_median)
    
    print('\033[1m\033[34m ' + "all LcbImage Object images have been simulated" + '\033[0m\n')