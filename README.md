# Sentinel-2 Top-Of-Atmosphere Radiometric Uncertainty Propagator (stoarup)

This repository implements the complete processing and workflow required to reproduce the study about the impact of
radiometric uncertainty in Sentinel-2 Top-of-Atmosphere data on the retrieval of land surface metrics.

## Citation

```latex
@article{graf_propagating_2023,
	title = {Propagating Sentinel-2 Top-of-Atmosphere Radiometric Uncertainty into Land Surface Phenology Metrics Using a Monte Carlo Framework},
	issn = {2151-1535},
	doi = {10.1109/JSTARS.2023.3297713},
	pages = {1--41},
	journaltitle = {{IEEE} Journal of Selected Topics in Applied Earth Observations and Remote Sensing},
    journal = {IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing},
	author = {Graf, Lukas Valentin and Gorroño, Javier and Hueni, Andreas and Walter, Achim and Aasen, Helge},
	date = {2023},
	keywords = {Crops, Ecophysiological Parameters, Indexes, Land Surface Phenology, Measurement, Monte-Carlo Simulation, Radiative Transfer Modelling, Radiometry, Remote sensing, Sentinel-2, {TIMESAT}, Uncertainty, Vegetation mapping},
    year = {2023}
}
```

## OS and Software Requirements

While the Python scripts are OS-independent, the shell-scripts only work in a *nix environment. They were developed and tested under Fedora 34 (64-bit) machine using Python 3.9.7.

See the `requirements.txt` file for a list of Python packages necessary to run the scripts. All dependencies should be installed beforehand into a clean virtual environment.

Furthermore, the Sentinel-2 Radiometric Uncertainty Toolbox and Sen2Cor are required.

For S2 data download, a [CREODIAS account](https://creodias.eu/) is required (free account). The API access token must be provided according to the [module doc string in the download module](src/processing/01_download_data_creodias.py).

### Installing E:earth_africa:dal for Raster Data Handling

For **reading and writing raster data** we use the E:earth_africa:dal Python package (Earth Observation Data Analysis Library). See [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6651767.svg)](https://doi.org/10.5281/zenodo.6651767) for more information.

Installing E:earth_africa:dal will setup most of the required Python dependencies. It is recommended to install E:earth_africa:dal into a clean virtual environment to avoid any undesired side-effects.

All requirements besides E:earth_africa:dal dependencies are listed in the [requirements.txt](requirements.txt) and can be installed into the virtual environment using

```{bash}
pip3 install -r requirements.txt
```

### Installing the Sentinel-2 Radiometric Uncertainty Toolbox

#### Sentinel-2 Toolbox

The installation of the [Sentinel-2 Uncertainty Toolbox](https://github.com/senbox-org/snap-rut) is a bit cumbersome (requires Java and Python). The following steps worked for me:

1. Install the [SNAP toolbox](https://step.esa.int/main/download/snap-download/). Select a Python virtual environment when asked during the installation process. See also [here](https://senbox.atlassian.net/wiki/spaces/SNAP/pages/19300362/How+to+use+the+SNAP+API+from+Python) how to connect SNAP with Python.
2. Install [maven](https://maven.apache.org/) if not yet installed. It is required to build the packages
3. Clone the [Sentinel-2 Uncertainty Toolbox](https://github.com/senbox-org/snap-rut) from Github and follow the instructions in the README. It also links to a wiki page explaining how to integrate the toolbox into SNAP.
4. Most likely the toolbox will throw an error when trying to run it. This is because [jpy](https://github.com/jpy-consortium/jpy) is missing but required by SNAP to connect Java and Python. To overcome the issue clone the Github repository and build it following the instructions in the README. Make sure that the environmental variables are set as specified in the Wiki. For Fedora I added the following lines of code into my `~/.bashrc` file:

```{bash}
export JAVA_HOME=$(dirname $(dirname $(readlink $(readlink $(which javac)))))
export JDK_HOME=/=${JAVA_HOME}
export PATH=$PATH:${JAVA_HOME}/bin
export PATH=$PATH:${JAVA_HOME}/include
```
Don't forget to run `source ~/.bashrc` after making the changes to apply them. Then continue with building jpy. Copy the wheel package (jpy-0.10.0.dev1-cp39-cp39-linux_x86_64.whl) into /home/$USER/.snap/snap-python/snappy/

The actual algorithm of the toolbox is coded in Python and takes its inputs from the metadata xml (MTD_MSIL1C.xml). The algorithm works bandwise and assign pixel uncertainties values between 0 and 250 (250 = 25% uncertainty or higher). Thus, a value of 12 corresponds to a uncertainty of 1.2%.

#### Python Integration/Usage

Since the core algorithm of the toolbox is written in Python it is necessary to get the Java-Python bridge of SNAP (snappy) running. To make `snappy` working follow these steps (shown for Fedora):

- Create a new (i.e., clean) virtual python virtual environment:
```{bash}
virtualenv snap-python
source snap-python/bin/activate
```
- Change into the SNAP installation directory. In my case this was found in: /home/graflu/.snap/snap-python/snappy/ and run:

```{bash}
python setup.py install
```
to install snappy into the environment. Additionally install numpy:

```{bash}
pip install numpy
```

Additionally, you have to enable the Python-Java bridge by building the jpy package.
To do so, clone the source from Github and build it:

```{bash}
git clone https://github.com/bcdev/jpy.git
cd jpy
python setup.py build maven bdist_wheel
```

Then copy the created wheel packages and the shared object files (*.so) from the /build directory of jpy into the snappy directory (top-level).

**IMPORTANT NOTE ON UPDATE OF PYTHON VERSION**: If you have upgraded your Python version (e.g., from Python 3.8 to 3.9, make sure to alter the filenames in ~/snap/snap-python/build/jpyconfig.properties) so that the point towards the updated shared objects and wheel packages. This affects the `jpy.jpyLib`, `jdl.jpyLib`, and `jpy.pythonLib` entries in that file. Otherwise you will get Java errors because the old shared objects are either deleted or not working any more!

To check if snappy was installed correctly start Python in interactive mode:

```{python}
(snap-python) [graflu@kp140-208 snappy]$ python
Python 3.9.7 (default, Aug 30 2021, 00:00:00) 
[GCC 11.2.1 20210728 (Red Hat 11.2.1-1)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import snappy
INFO: org.esa.snap.python.gpf.PyOperatorSpi: Python operator 'S2RutOp' registered (Python module: 's2_rut', class: 'S2RutOp', root: '/mnt/ides/Lukas/software/snap/s2tbx/modules/org-esa-snap-snap-rut.jar')
INFO: org.esa.s2tbx.dataio.gdal.GDALVersion: GDAL 3.0.4 found on system. JNI driver will be used.
INFO: org.esa.s2tbx.dataio.gdal.GDALVersion: Installed GDAL 3.0.4 set to be used by SNAP.
INFO: org.esa.snap.core.gpf.operators.tooladapter.ToolAdapterIO: Initializing external tool adapters
INFO: org.esa.snap.core.util.EngineVersionCheckActivator: Please check regularly for new updates for the best SNAP experience.
```

### Installing Sen2Cor

Download Sen2Cor v2.9 from the [official access point](http://step.esa.int/main/snap-supported-plugins/sen2cor/sen2cor-v2-9/).

The `path` to the Sen2Cor executable must be provided in the shell script [05_execute_sen2cor.sh](src/processing/05_execute_sen2cor.sh#L59) or add the path to the binary directory to your $PATH.

## Filesystem Structure

The table below explains the file-system of this repository which will hold the Monte Carlo simulations and analysis results. `static` in the column `content type` means that the content is provided as is when cloning this repo, whereas `dynamic` means that the content is generated during execution of the workflow.

| sub-directory          | purpose                                                                                                         | content type |
|------------------------|-----------------------------------------------------------------------------------------------------------------|--------------|
| src                    | contains the Python source code organized by single modules (*.py)                                              | static       |
| log                    | all log information is placed here (ASCII files)                                                                | dynamic      |
| S2_MSIL1C_orig         | here, the original S2 data in L1C processing level downloaded from CREODIAS is stored                           | dynamic      |
| S2_MSIL1C_RUT-Scenario | sub-directory where the Monte-Carlo L1C scenarios are stored for each S2 scene                                  | dynamic      |
| S2_MSIL2A_Analysis     | contains the uncertainty analysis results after propagating uncertainty through Sen2Cor into EVI, NDVI and GLAI | dynamic      |
| S2_ProSAIL_LUTs        | contains lookup-tables of ProSAIL forward runs for each scene (50 000 spectra)                                  | static       |
| S2_TimeSeriesAnalysis  | sub-directory where the phenology time series Monte Carlo simulations and analysis results are stored           | dynamic      |
| shp                    | contains the shapefiles of the study area, crop type map and sampling points.                                   | static       |

## Executing the Workflow

To re-generate the results presented in the paper, execute the Python modules and shell scripts in the [processing sub-package](src/processing/) in the order the scripts are named: I.e., start with the [downloader script](src/processing/01_download_data_creodias.py) (named `01_download_data_creodias.py`), continue with [02_write_property_file.py](src/processing/02_write_property_file.py) and so on until you reach [the last script](src/processing/11_analyze_l4_uncertainty.sh).

**IMPORTANT**
* You have to provide the path to the SNAP graph processing tool (gpt) executable in the [S2RUT shell script](src/processing/03_s2_radiometric_uncertainty.sh)
* You have to provide the path to the L2A_Bashrc file of the Sen2Cor installation in the [shell script running Sen2Cor](src/processing/05_execute_sen2cor.sh) and [here](src/processing/05_execute_sen2cor_orig-data.sh)

**CAUTION**:
    Some of the scripts require **several days to terminate** and consume reasonable amounts of CPU time, RAM and disk storage (e.g., running Sen2Cor several hundred times). Therefore, we also do not provide a single shell script to execute everything at once. Instead we recommend careful testing beforehand of the scripts (e.g., by lowering the number of scenarios generated, starting with a single S2 scene, etc.)
