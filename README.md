# HORAYZON

Package to efficiently compute terrain parameters (like **horizon**, **sky view factor**, **topographic openness**, slope angle/aspect) from high-resolution digital elevation model (DEM) data. The package also allows to compute **shadow maps** (and **correction factors for downwelling direct shortwave radiation**) for specific sun positions. Horizon computation is based on the high-performance ray-tracing library Intel&copy; Embree. Calculations are parallelised with OpenMP (Cython code) or Threading Building Blocks (C++ code).

When you use HORAYZON, please cite:

**Steger, C. R., Steger, B. and Schär, C (2022): HORAYZON v1.0: An efficient and flexible ray-tracing algorithm to compute horizon and sky view factor, Geosci. Model Dev., https://doi.org/10.5194/gmd-2022-58**

# Package dependencies

HORAYZON depends on multiple external libraries and packages. The essential ones are listed below under **Core dependencies**. Further listed dependencies are only needed to run the examples. The examples **horizon/gridded_curved_DEM_masked.py**, **horizon/gridded_planar_DEM_2m.py** and **shadow/gridded_curved_DEM_NASADEM.py** require more complex dependencies, which are listed under **Specific dependencies for examples**. It is recommended to handle dependencies via [Conda](https://docs.conda.io/en/latest/#), which covers all dependencies except **hmm**. A new Conda environment can be created according to the below examples.

**Core dependencies**
- [Intel Embree](https://www.embree.org) and [Threading Building Blocks (TBB)](https://github.com/oneapi-src/oneTBB)
- [NetCDF-4 C++ library](https://github.com/Unidata/netcdf-cxx4)
- Python packages: Cython, NumPy, SciPy, GeographicLib, tqdm, requests, xarray
```bash
conda create -n horayzon_core -c conda-forge embree3 tbb-devel netcdf-cxx4 cython numpy scipy geographiclib tqdm requests xarray
```

**Base dependencies of examples**
- Python packages: netCDF4, Matplotlib, Pillow, Skyfield, pyproj, IPython
```bash
conda create -n horayzon_base -c conda-forge embree3 tbb-devel netcdf-cxx4 cython numpy scipy geographiclib tqdm requests xarray netcdf4 matplotlib pillow skyfield pyproj ipython
```

**Specific dependencies for examples (masking and high-resolution DEM examples; GDAL dependency)**
- Python packages: Shapely, fiona, PyGEOS, scikit-image, Rasterio, Trimesh
- [heightmap meshing utility (hmm)](https://github.com/fogleman/hmm)
```bash
conda create -n horayzon_all -c conda-forge embree3 tbb-devel netcdf-cxx4 cython numpy scipy geographiclib tqdm requests xarray netcdf4 matplotlib pillow skyfield pyproj ipython shapely fiona pygeos scikit-image rasterio trimesh
```
Installation instruction for **hmm** can be found [here](https://github.com/fogleman/hmm). **hmm**'s dependency **glm** can be installed via a package manager (e.g. APT, MacPorts, Homebrew) or via manual building from [source code](https://glm.g-truc.net/0.9.9/index.html).

# Installation

HORAYZON has been tested on **Python 3.10** under **Linux** and **Mac OS X**. Installation requires the [GNU Compiler Collection (GCC)](https://gcc.gnu.org) and can be accomplished as follows:

## Linux

Create an appropriate Conda environment (see examples above) and activate this environment. The HORAYZON package can then be installed with:
```bash
git clone https://github.com/ChristianSteger/HORAYZON.git
cd HORAYZON  
python -m pip install .
```

## Mac OS X
Create an appropriate Conda environment (see examples above) and activate this environment. Download the HORAYZON package and change to the main directory:
```bash
git clone https://github.com/ChristianSteger/HORAYZON.git
cd HORAYZON  
```
Currently, the default method of installing HORAYZON under Mac OS X fails because the NetCDF-4 C++ library provided by Conda was build with the C++ Standard Library **libc++**, which causes an incompatibility when the source code is compiled with GCC. It is therefore necessary to build and install the NetCDF-4 C++ library with GCC as well. This is for instance possible via MacPorts or Homebrew. The appropriate command for MacPorts is:

```bash
sudo port install netcdf-cxx4 +gcc10 
```
Subsequently, the path to the newly installed NetCDF-4 C++ library has to be set in HORAYZON's setup file (**setup.py**) and the package can be installed with:
```bash 
python -m pip install .
```

# Usage

The usage of the packages is best illustrated by means of examples, which can either be run in a Python IDE (like PyCharm or Spyder) or in the terminal. To run the examples, the path **path_out** must be adapted to a location that provides enough disk space. All input data (DEM or auxiliary data; see section below) required for running the examples is downloaded automatically.

## Examples: Terrain parameters (slope, horizon and sky view factor)

Two terrain horizon functions are available, **horizon_gridded()** and **horizon_locations()**. The former function allows to compute horizon for gridded input while the latter allows to compute horizon for arbitrary selected locations. The second function can optionally output the distance to the horizon. The following five examples are provided:
- **horizon/gridded_curved_DEM.py**: Compute topographic parameters from SRTM (geodetic coordinates, ~90 m resolution) for a ~50x50 km example region in the European Alps. Earth's surface curvature is considered. Plot output of this script is shown below.
![Alt text](https://github.com/ChristianSteger/Images/blob/master/Topo_slope_SVF_new.png?raw=true "Output from horizon/gridded_curved_DEM.py")
- **horizon/gridded_planar_DEM.py**: Compute topographic parameters from swisstopo DHM25 (map projection, 25 m resolution) for a ~25x40 km example region in Switzerland. Earth's surface curvature is neglected.
- **horizon/locations_curved_DEM.py**: Compute topographic parameters (additionally distance to horizon) from SRTM (geodetic coordinates, ~90 m resolution) for 11 locations in Switzerland. Earth's surface curvature is considered. Plot output of this script for one location is shown below.
![Alt text](https://github.com/ChristianSteger/Images/blob/master/Wengen.png?raw=true "Output from horizon/locations_curved_DEM.py")
- **horizon/gridded_curved_DEM_masked.py**: Compute topographic parameters from SRTM (geodetic coordinates, ~90 m resolution) for South Georgia in the South Atlantic Ocean. Earth's surface curvature is considered. DEM grid cells, which are at least 20 km apart from land, are masked to speed-up the horizon computation.
- **horizon/gridded_planar_DEM_2m.py**:  Compute gridded topographic parameters from swissALTI3D (map projection, 2 m resolution) for a 3x3 km example region in Switzerland. Earth's surface curvature is neglected. The outer DEM domain is simplified and represented by a triangulated irregular network (TIN) to reduce the large memory footprint of the DEM data.


**A Note on Sky view factor and related parameters**<br/>
The term sky view factor (SVF) is defined ambiguously in literature. In Zakšek et al. (2011), it referes to the solid angle of the (celestial) hemisphere. We call this parameter *visible sky fraction* and its computation is performed with the function **topo_param.visible_sky_fraction()**. In applications related to radiation, the SVF is typically defined as the fraction of sky radiation received at a certain location in case of isotropic sky radiation (see e.g. Helbig et al., 2009). This parameter is called *sky view factor* in our application and its computation is performed with the function **topo_param.sky_view_factor()**. Additionally, the positive topographic openness (Yokoyama et al., 2002) can be computed with the function **topo_param.topographic_openness()**. 

## Examples: Shadow map and shortwave correction factor

The module **shadow** allows to compute shadow maps and correction factors for downwelling direct shortwave radiation for terrains and specific sun positions. A class **shadow.Terrain** can be created and initialised for this purpose, which can then be used for arbitrary sun positions. The output of the method **Terrain.shadow()** is encoded as follows: 0: illuminated, 1: self-shaded, 2: terrain-shaded, 3: not considered (respectively masked).
The correction factors for downwelling direct shortwave radiation is computed with the method **Terrain.sw_dir_cor()** according to Müller and Scherer (2005). This factor can for instance be applied to radiation output from a regional climate or general circulation model, in which radiation is only simulated along the vertical dimension and all grid cells are assumed to have a horizontal surface. The correction factor accounts for all terrain-induced modifications in radiation, like self-shading, terrain-shading, changes in angles between the sun and the surface's normal vector and the surface enlargement of grid cells due to sloping surfaces. According to Müller and Scherer (2005), it is computed as:

$f_{cor} = \left( \dfrac{1.0}{\vec{a} \times \vec{b}} \right) $

Test

- **shadow/gridded_curved_DEM_SRTM.py**:
- **shadow/gridded_curved_DEM_REMA.py**:
- **shadow/gridded_planar_DEM_artifical.py**:
- **shadow/gridded_curved_DEM_NASADEM.py**

# Digital elevation model and auxiliary data

Digital elevation model (DEM) data is available from various sources, e.g.:
- [SRTM](https://srtm.csi.cgiar.org)
- [DHM25](https://www.swisstopo.admin.ch/en/geodata/height/dhm25.html)
- [swissALTI3D](https://www.swisstopo.admin.ch/en/geodata/height/alti3d.html)
- [NASADEM](https://search.earthdata.nasa.gov/)
- [MERIT](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)
- [USGS 1/3rd arc-second DEM](https://www.sciencebase.gov/catalog/item/4f70aa9fe4b058caae3f8de5)

Auxiliary data, like geoid undulation data (EGM96 and GEOID12A), coastline polygons (GSHHG) or glacier outlines (GAMDAM) are for instance available here:
- [EGM96](https://earth-info.nga.mil)
- [GEOID12A](https://geodesy.noaa.gov/GEOID/GEOID12A/GEOID12A_AK.shtml)
- [GSHHG](https://www.soest.hawaii.edu/pwessel/gshhg/)
- [GAMDAM](https://doi.pangaea.de/10.1594/PANGAEA.891423)

# References
- Helbig, N., Löwe, H. and Lehning, M. (2009): Radiosity Approach for the Shortwave Surface Radiation Balance in Complex Terrain, Journal of the Atmospheric Sciences, 66(9), 2900-2912, https://doi.org/10.1175/2009JAS2940.1
- Yokoyama, R., Shirasawa, M. and Pike, R. J. (2002): Visualizing Topography by Openness: A New Application of Image Processing to Digital Elevation Models, Photogrammetric Engineering and Remote Sensing, 68, 257-265.
- Zakšek, K., Oštir, K. and Kokalj, Ž. (2011): Sky-View Factor as a Relief Visualization Technique, Remote Sensing, 3(2):398-415, https://doi.org/10.3390/rs3020398
- Müller, M. D., and Scherer, D. (2005): A Grid- and Subgrid-Scale Radiation Parameterization of Topographic Effects for Mesoscale Weather Forecast Models, Monthly Weather Review, 133(6), 1431-1442, https://journals.ametsoc.org/view/journals/mwre/133/6/mwr2927.1.xml

# Support 
In case of issues or questions, please contact Christian Steger (christian.steger@env.ethz.ch).