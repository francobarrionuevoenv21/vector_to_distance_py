# v2d — Vector to Distance

A lightweight Python package to compute a Euclidean **distance raster**,
in meters, from a vector layer based on built-in SciPy tools.

Given a vector layer (e.g. urban areas, roads, water bodies), `v2d` rasterizes it,
reprojects it automatically to the correct UTM zone (in meters) based on its location,
and computes the distance from every pixel to the nearest target feature.

## Description

The package works in four steps, all wrapped into a single function
(`vector_to_distance`):

1. **Read & reproject** — reads the input vector file and reprojects it to the UTM
   zone (in meters) matching its centroid location, so distances are computed in
   real-world units regardless of the vector's original CRS.
2. **Rasterize** — burns the vector geometries into a raster array, marking target
   pixels (e.g. `1`) and background pixels (e.g. `0`) at a user-defined pixel size.
3. **Distance transform** — computes the Euclidean distance, in meters, from every
   pixel to the nearest target pixel.
4. **Export** — writes the resulting distance raster to disk as a single-band,
   float32 GeoTIFF.

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/francobarrionuevoenv21/vector_to_distance_py 
cd vector_to_distance_py
pip install -r requirements.txt
```

## Usage

```python
import v2d

dist_array = v2d.vector_to_distance(
    vector_path='vector.geojson',   # Input vector file (All GeoPandas supported formats are accepted) 
    px_m=30,                        # Pixel size in meters
    file_output='output',           # Output file name, without extension (Default: 'output')
    target=1,                       # Value assigned to vector features (Default: 1)
    no_target=0,                    # Value assigned to background pixels (Default: 0)
)
```

Running the example above will:

- Read `vector.geojson` and reproject it to the correct UTM zone in meters.
- Rasterize it at 30 m resolution.
- Compute the distance (in meters) from every pixel to the nearest urban area.
- Save the result as `output.tif` in the (upper) current working directory.
- Return the distance array as a NumPy array (`dist_array`) for further use in
  Python (plotting, further analysis, etc.), in addition to the exported file.

### Using individual steps

Each step is also available on its own, in case you want more control over the
pipeline (e.g. reusing an already-reprojected vector, or exporting your own profile):

```python
import v2d

vector_gdf = v2d.read_reprj_m('vector.shp')
vector_array, profile = v2d.vector_to_raster(vector_gdf, px_m=30, target=1, no_target=0)
dist_array = v2d.array_distance(vector_array, px_m=30, no_target=0)
v2d.export_raster(dist_array, profile, 'output')
```

## Output

- **Return value**: a NumPy array (`float32`) of pixel-wise Euclidean distances, in
  meters, to the nearest target feature.
- **File output**: a single-band GeoTIFF (`<file_output>.tif`) with:
  - `dtype`: `float32`
  - `nodata`: `-9999`
  - CRS matching the automatically selected UTM zone (in meters)
  - Pixel size matching the `px_m` parameter used

## Requirements

- Python 3.9+
- geopandas
- numpy
- rasterio
- affine
- scipy

(see `requirements.txt`)
