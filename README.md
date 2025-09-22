# `megaradrpsimul`

### Code for simulating and reducing MEGARA raw images.

To install the package in your environment, first clone the repository from GitHub by running in your terminal:

```text
$ git clone https://github.com/mchillaron/megaradrpsimul.git
```

Then, navigate into the megaradrpsimul/ folder:
```text
$ cd megaradrpsimul/
```
And install the package in editable mode:
```text
$ pip install -e .
```
The code will be now ready to use!

### Command line usage
```text
$ megaradrp-simulation -h
megaradrp-simulation [-h] [--obj_name OBJ_NAME] [--vph VPH]
                            [-c CONFIG_FILE] [-n NUM_SIMUL] [--run_modelmap]
                            [--run_twilight]

Simulate MEGARA reductions.

optional arguments:
  -h, --help            show this help message and exit
  --obj_name OBJ_NAME   Name of the object to simulate.
  --vph VPH             VPH name.
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Name of the configuration file.
  -n NUM_SIMUL, --num_simul NUM_SIMUL
                        Number of simulations to perform.
  --run_modelmap        Run ModelMap step.
  --run_twilight        Run Twilight step.
  --pixel_size PIXEL_SIZE
                        Pixel size in arcseconds for the conversion of RSS
                        into a cube.
  ```
