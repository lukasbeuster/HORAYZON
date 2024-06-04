import os
import numpy as np
from importlib import import_module
import xarray as xr
import glob
from osgeo import gdal, osr



def dsm(file_dem, engine="gdal"):
    """Load SRTM digital elevation model data.

    Load SRTM digital elevation model data from single GeoTIFF file.

    Parameters
    ----------
    file_dem : str
        Name of SRTM tile
    domain : dict
        Dictionary with domain boundaries (lon_min, lon_max, lat_min, lat_max)
        [degree]
    engine: str
        Backend for loading GeoTIFF file (either 'gdal' or 'pillow')

    Returns
    -------
    lon : ndarray
        Array (one-dimensional) with longitude [degree]
    lat : ndarray
        Array (one-dimensional) with latitude [degree]
    elevation : ndarray
        Array (two-dimensional) with elevation [metre]

    Notes
    -----
    Data source: https://srtm.csi.cgiar.org"""

    # Check arguments
    if engine not in ("gdal", "pillow"):
        raise ValueError("Input for 'engine' must be either "
                         "'gdal' or 'pillow'")

    # Load digital elevation model data
    if engine == "gdal":
        print("Read GeoTIFF with GDAL")
        ds = gdal.Open(file_dem)
        
        # Warp to WGS84 (because in 28992?)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)  # WGS84
        #ds = gdal.Warp('', ds, dstSRS=srs)
        ## Creates a vrt file in memory so the warp result does not need to be saved to disk.
        ds = gdal.Warp('', ds, format='vrt', dstSRS=srs)
        
        elevation = ds.GetRasterBand(1).ReadAsArray()  # 16-bit integer
        ## Added no data pull so we can mask nodata
        nd = ds.GetRasterBand(1).GetNoDataValue()
        raster_size_x, raster_size_y = ds.RasterXSize, ds.RasterYSize
        lon_ulc, lat_ulc = ds.GetGeoTransform()[0], ds.GetGeoTransform()[3]
        d_lon, d_lat = ds.GetGeoTransform()[1], ds.GetGeoTransform()[5]
    else:
        print("Read GeoTIFF with Pillow")
        if (os.path.getsize(file_dem) / (1024 ** 2)) > 500.0:
            print("Warning: reading of large GeoTIFF file with Pillow is slow")
        Image = import_module("PIL.Image")
        Image.MAX_IMAGE_PIXELS = 1300000000
        img = Image.open(file_dem)
        elevation = np.array(img)  # 32-bit integer
        raster_size_x, raster_size_y = img.tag[256][0], img.tag[257][0]
        lon_ulc, lat_ulc = img.tag[33922][3], img.tag[33922][4]
        d_lon, d_lat = img.tag[33550][0], -img.tag[33550][1]
        # Warning: unclear where sign of n-s pixel resolution is stored!
    lon_edge = np.linspace(lon_ulc, lon_ulc + d_lon * raster_size_x,
                           raster_size_x + 1)
    lat_edge = np.linspace(lat_ulc, lat_ulc + d_lat * raster_size_y,
                           raster_size_y + 1)
    lon = lon_edge[:-1] + np.diff(lon_edge / 2.0)
    lat = lat_edge[:-1] + np.diff(lat_edge / 2.0)

    # Domain size and computation settings
    domain = {
    "lon_min": lon.min(),
    "lon_max": lon.max(),
    "lat_min": lat.min(),
    "lat_max": lat.max(),
    }
    
    # Define a margin (for example, 1% of the total range)
    lat_margin = 0.1 * (domain["lat_max"] - domain["lat_min"])
    lon_margin = 0.1 * (domain["lon_max"] - domain["lon_min"])

    # Reduce the domain slightly
    reduced_domain = {
        "lon_min": domain["lon_min"] + lon_margin,
        "lon_max": domain["lon_max"] - lon_margin,
        "lat_min": domain["lat_min"] + lat_margin,
        "lat_max": domain["lat_max"] - lat_margin,
    }

#### Changes: Commented out cropping domain since we're using full domain for the urban tiles. Might have to be changed when I work with multiple tiles and consider the overlap LB 240604

    # Crop relevant domain
    # if any([domain["lon_min"] < lon_edge.min(),
    #         domain["lon_max"] > lon_edge.max(),
    #         domain["lat_min"] < lat_edge.min(),
    #         domain["lat_max"] > lat_edge.max()]):
    #     raise ValueError("Provided tile does not cover domain")
    # slice_lon = slice(np.where(lon_edge <= domain["lon_min"])[0][-1],
    #                   np.where(lon_edge >= domain["lon_max"])[0][0])
    # slice_lat = slice(np.where(lat_edge >= domain["lat_max"])[0][-1],
    #                   np.where(lat_edge <= domain["lat_min"])[0][0])
    # elevation = elevation[slice_lat, slice_lon].astype(np.float32)
    # lon, lat = lon[slice_lon], lat[slice_lat]

    # print_dem_info(elevation)

    return lon, lat, elevation, nd, reduced_domain