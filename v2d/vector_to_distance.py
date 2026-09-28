"""
-*- coding: utf-8 -*-
----------------------------------------------------------------------------
Created By: Franco Barrionuevo
Created Date: September 2026
----------------------------------------------------------------------------

v2d-Vector to Distance

Utilities for computing a Euclidean proximity (distance) raster, in meters,
from a vector layer.

Pipeline overview
------------------
1. `get_utm_epsg`     -> determine the UTM zone (EPSG code) matching the vector's location
2. `read_reprj_m`     -> read the vector file, reproject to that UTM zone (meters)
3. `get_transform`    -> build the raster's affine transform from the vector's bounds
4. `get_profile`      -> build the rasterio profile (driver, dtype, nodata, crs, transform)
5. `vector_to_raster` -> rasterize the vector into a target / no_target array
6. `array_distance`   -> Euclidean distance transform, in meters, to the nearest target
7. `export_raster`    -> write the distance array to disk as a GeoTIFF
8. `vector_to_distance` -> end-to-end pipeline, exports the result as a GeoTIFF
"""

# Import libraries
import geopandas as gpd
import numpy as np
import rasterio
from affine import Affine
from rasterio.features import rasterize
from scipy.ndimage import distance_transform_edt


def get_utm_epsg(vector_gdf):
    """
    Determine the UTM EPSG code (in meters) matching a vector layer's centroid location.

    Parameters
    ----------
    vector_gdf : geopandas.GeoDataFrame
        Vector layer in any CRS.

    Returns
    -------
    str
        EPSG code (as 'EPSG:XXXXX') of the UTM zone matching the vector's centroid location.

    Reference
    ---------
    https://sis.apache.org/howto/instantiate_utm_projection.html
    """
    # Convert to EPSG:4326
    vector_4326 = vector_gdf.to_crs('EPSG:4326')
    # Extract boundaries coordinates
    minx, miny, maxx, maxy = vector_4326.total_bounds
    # Compute mean longitude and latitude
    lon_mean = (minx + maxx) / 2
    lat_mean = (miny + maxy) / 2
    # Compute UTM zone number from longitude
    zone = int((lon_mean + 180) / 6) + 1
    # Select EPSG base according to hemisphere
    epsg_base = 32700 if lat_mean < 0 else 32600
    return f'EPSG:{epsg_base + zone}'


def read_reprj_m(vector_path):
    """
    Read a vector file and reproject it to its matching UTM CRS (in meters).

    Parameters
    ----------
    vector_path : str
        Path to the input vector file.

    Returns
    -------
    geopandas.GeoDataFrame
        Vector reprojected to the UTM CRS (in meters) matching its centroid location.
    """
    # Read vector file
    gdf = gpd.read_file(vector_path)
    # Get UTM CRS according to centroid location
    crs = get_utm_epsg(gdf)
    return gdf.to_crs(crs)


def get_transform(vector_gdf_reprjm, px_m):
    """
    Build the affine transform for a raster aligned to a vector layer's bounds.

    Parameters
    ----------
    vector_gdf_reprjm : geopandas.GeoDataFrame
        Vector layer reprojected to a projected (meters) CRS.
    px_m : float
        Pixel size in meters (square pixels).

    Returns
    -------
    affine.Affine
        Affine transform for the output raster, aligned to the vector's bounds.
    """
    # Define squared pixel
    pxx_m = pxy_m = px_m

    # Get reprojected vector boundaries coordinates
    minx, miny, maxx, maxy = vector_gdf_reprjm.total_bounds

    # Define transform matrix values
    a = pxx_m   # pixel width (x-resolution) — size of a pixel in the x-direction
    b = 0       # row rotation (usually 0 for north-up rasters)
    c = minx    # x-coordinate of the upper-left corner of the upper-left pixel
    d = 0       # column rotation (usually 0 for north-up rasters)
    e = -pxy_m  # pixel height (y-resolution) — usually negative, since y decreases as row index increases (raster origin is top-left, but y increases upward in most CRS)
    f = maxy    # y-coordinate of the upper-left corner of the upper-left pixel

    return Affine(a, b, c, d, e, f)


def get_profile(shape, crs, transform):
    """
    Build a rasterio profile for a single-band float32 GeoTIFF.

    Parameters
    ----------
    shape : tuple of int
        Raster shape as (rows, cols).
    crs : rasterio.crs.CRS or str
        Coordinate reference system of the raster.
    transform : affine.Affine
        Affine transform of the raster.

    Returns
    -------
    dict
        Rasterio profile ready to be passed to rasterio.open() as kwargs.
    """
    return {'driver': 'GTiff',
            'height': shape[0],
            'width': shape[1],
            'count': 1,
            'dtype': 'float32',
            'nodata': -9999,
            'crs': crs,
            'transform': transform}


def vector_to_raster(vector_gdf_reprjm, px_m, target, no_target):
    """
    Rasterize a vector layer into a binary target/no_target array.

    Parameters
    ----------
    vector_gdf_reprjm : geopandas.GeoDataFrame
        Vector layer reprojected to a projected (meters) CRS.
    px_m : float
        Pixel size in meters (square pixels).
    target : int
        Value assigned to rasterized vector features.
    no_target : int
        Value assigned to background (non-feature) pixels.

    Returns
    -------
    tuple
        (numpy.ndarray, dict) — the rasterized array and its matching rasterio profile.
    """
    # Get transform according the vector
    transform = get_transform(vector_gdf_reprjm, px_m)

    # Compute output raster shape
    minx, miny, maxx, maxy = vector_gdf_reprjm.total_bounds
    pxy_m = pxx_m = px_m
    shape = (int((maxy - miny) // pxy_m), int((maxx - minx) // pxx_m))

    # Convert vector to raster
    vector_array = rasterize(
        [(shape, label) for shape, label in zip(vector_gdf_reprjm['geometry'],
                                                  np.ones(vector_gdf_reprjm.shape[0], dtype=np.int16) * target)],
        out_shape=shape,
        transform=transform,
        fill=no_target,
        all_touched=True,
        dtype=rasterio.int16)

    return vector_array, get_profile(shape, vector_gdf_reprjm.crs, transform)


def array_distance(vector_array, px_m, no_target):
    """
    Compute the Euclidean distance, in meters, from each pixel to the nearest target pixel.

    Parameters
    ----------
    vector_array : numpy.ndarray
        Rasterized array containing target and no_target values.
    px_m : float
        Pixel size in meters (square pixels).
    no_target : int
        Value in vector_array representing background (non-target) pixels.

    Returns
    -------
    numpy.ndarray
        Array of Euclidean distances, in meters, from each pixel to the nearest target pixel.
    """
    # Get array of distances
    # distance_transform_edt computes distance to nearest 'no_target' from any 'target'. The method takes 0 as no-target by default
    dist_array_px = distance_transform_edt(vector_array == no_target).astype(np.float32)  # distance in pixels
    dist_array_m = dist_array_px * px_m

    return dist_array_m


def export_raster(dist_array_m, profile, file_output):
    """
    Write a distance array to disk as a single-band GeoTIFF.

    Parameters
    ----------
    dist_array_m : numpy.ndarray
        Distance array, in meters, to export.
    profile : dict
        Rasterio profile matching dist_array_m.
    file_output : str
        Output file path without extension.

    Returns
    -------
    None
    """
    # Export distance array as tif file with rasterio
    with rasterio.open(f'../../{file_output}.tif', 'w', **profile) as dst:
        dst.write(dist_array_m, 1)
        
    # Display in console completed file storaged
    print(f'Output file succesfully stored in ../../{file_output}.tif ✅')
    


def vector_to_distance(vector_path, px_m, file_output='output', target=1, no_target=0):
    """
    End-to-end pipeline: read a vector file, rasterize it, compute a Euclidean
    proximity raster in meters, and export it as a GeoTIFF (QGIS "Proximity"-style tool).

    Parameters
    ----------
    vector_path : str
        Path to the input vector file.
    px_m : float
        Pixel size in meters (square pixels).
    file_output : str
        Output file path without extension.
    target : int
        Value assigned to rasterized vector features.
    no_target : int
        Value assigned to background (non-feature) pixels.

    Returns
    -------
    numpy.ndarray
        Array of Euclidean distances, in meters, from each pixel to the nearest target pixel.
    """
    # Read the vector and reproject to UTM
    vector_gdf_reprjm = read_reprj_m(vector_path)

    # Get vector array and profile
    vector_array, profile = vector_to_raster(vector_gdf_reprjm, px_m, target, no_target)

    # Get distance array in meters (m)
    dist_array_m = array_distance(vector_array, px_m, no_target)

    # Export distance raster
    export_raster(dist_array_m, profile, file_output)

    return dist_array_m
