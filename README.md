# README #

## Basic usage
### Pre-process and filtering `MATLAB`
 - `make_fluo_data.m`: Launches the pre-processing of the fluorescent Calcium data with the parameters

 - `preview_image_stack.m`: Crop, rotate, and other basic pre-processing functions.
                       It requires a basic set of parameters.
                     
 - `process_image_stack.m`: Processes the whole image stack according to the
                            predefined parameters. It creates a file with the
                            filtered data and calls the plotting functions when
                            required.
 
### Post-processing and visualization `python`
 - `fluo/basic_image_visualization.py`: Basic spatial filtering and visualization. 
                                        Stimulation cycle averaging included.

### How do I get set up? ###
* Get this repo
    - `git clone git@github.com:xxx`
    - `cd hTUS`
    - `git checkout develop`
    - `git config pull.rebase false`

* To use the python scripts:
    - `cd hTUS`
    - `python3 -m venv venv` 
    - `source venv/bin/activate`
    - `pip -V`
    - `python3 -m pip install --upgrade pip`
    - `pip install scipy`
    - `pip install matplotlib`
    - `pip install ipython[all]`
    - `pip install h5py`
    - `pip install scikit-image`
    - Use `deactivate` to deactivate the python virtual environment
    - Create a file called fluo_los.sh on your main folder of the pc (or whatever folder where the command window will start) to automate starting python with the following content:
    ```
    cd hTUS path
    source venv/bin/activate
    ipython --pylab
    ```    
