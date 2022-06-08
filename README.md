# Radiometric Atmospheric Correction Uncertainty Propagation Study

This repository implements the **complete processing and workflow required to reproduce our study about the impact of
radiometric uncertainty on the retrieval of key phenological stages from Sentinel-2 data**.

## OS and Software Requirements

While the Python scripts are OS-independent, the shell-scripts only work in a *nix environment. They were developed and tested under Fedora 34 (64-bit) machine using Python 3.9.7.

See the `requirements.txt` file for a list of Python packages necessary to run the scripts. All dependencies should be installed beforehand into a clean virtual environment.

Furthermore, the Sentinel-2 Radiometric Uncertainty Toolbox and Sen2Cor are required.

For **reading and writing raster data** we use the `eodal` Python package (Earth Observation Data Analysis Library). See ... for more information.

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
Don't forget to run `source ~\.bashrc` after making the changes to apply them. Then continue with building jpy. Copy the wheel package (jpy-0.10.0.dev1-cp39-cp39-linux_x86_64.whl) into /home/$USER/.snap/snap-python/snappy/

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

| directory                 | purpose                                                                                                                                                                                           |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|             |
| /shp                      | contains the shapefile defining the extent of the study area (AOI_Esch_EPSG32632.shp) and the single regions of interest (ROIs) within that area (ZH_Polygons_2019_EPSG32632_selected-crops.shp). |
| /src						| here, all the required Python and shell scripts are located.
| /S2A_MSIL1C_orig          | here, the original L1C scenes will be downloaded to (from Creodias).                                                                                                                              |
| /S2A_MSIL1C_RUT-Scenarios | here, the L1C scenarios (based on the radiometric uncertainty assessment) and the resulting L2A outputs (uncertainty propagation) will be stored.                                                 
| /S2A_MSIL2A_Analysis      | here, the results of the analysis of the L1C and L2A scenarios will be stored

## Executing the Workflow

coming soon

