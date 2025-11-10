import os
import numpy as np
from shapely import line_interpolate_point, LineString, MultiLineString
from shapely.geometry import Point,box
import math
import geopandas as gpd
import pandas as pd
import openrouteservice
from rasterstats import zonal_stats
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject
from rasterio.enums import Resampling
import folium
from tqdm import tqdm
import os
import zipfile
import requests
import urllib
from tqdm import tqdm

def download_and_extract_roads(country_roads_shp, data_path, path_roads):
    """
    Download and extract a country's roads shapefile if it does not already exist locally.
    
    This function checks whether a local directory (`path_roads`) exists.
    If not, it downloads a zipped shapefile from the given URL (`country_roads_shp`),
    displays a progress bar during download, extracts the contents into a "Roads" folder
    inside `data_path`, and removes any files or directories that do not contain 'roads'
    in their name.
    
    Parameters
    ----------
    country_roads_shp : str
        The URL pointing to the zipped shapefile of the country's roads.
        
    data_path : str
        The base directory where the data will be downloaded and extracted.
        
    path_roads : str
        The expected path of the extracted roads data. If this path exists,
        the function will skip the download and extraction.
    
    Returns
    -------
    None
        This function performs file system operations and does not return a value.
    
    Notes
    -----
    - Uses the `requests` library for streaming download.
    - Requires the `tqdm` library for displaying a progress bar.
    - Automatically deletes temporary files and non-road files after extraction.
    
    Examples
    --------
    >>> download_and_extract_roads(
    ...     country_roads_shp="https://example.com/roads.zip",
    ...     data_path="/data/geodata",
    ...     path_roads="/data/geodata/Roads/roads.shp"
    ... )
    Downloading Roads: 100%|██████████| 25.3M/25.3M [00:04<00:00, 6.2MB/s]
    ✅ Roads shapefile downloaded and extracted successfully.
    """
    
    if not os.path.exists(path_roads):
        # Stream download with a progress bar
        response = requests.get(country_roads_shp, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024  # 1 KB per iteration
        temp_zip_path = os.path.join(data_path, 'temp_roads.zip')
        
        # Write the file in chunks while updating the progress bar
        with open(temp_zip_path, 'wb') as file, tqdm(
            total=total_size, unit='B', unit_scale=True, desc='Downloading Roads'
        ) as progress_bar:
            for chunk in response.iter_content(block_size):
                progress_bar.update(len(chunk))
                file.write(chunk)
        
        # Extract ZIP contents
        with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
            extract_dir = os.path.join(data_path, 'Roads')
            zip_ref.extractall(extract_dir)
        
        # Remove non-road files or directories
        for dir_name in os.listdir(extract_dir):
            if 'roads' not in dir_name.lower():
                try:
                    os.remove(os.path.join(extract_dir, dir_name))
                except IsADirectoryError:
                    pass  # Skip directories safely
        
        # Clean up temporary ZIP
        os.remove(temp_zip_path)
        
        print("✅ Roads shapefile downloaded and extracted successfully.")


def compute_isochrones(gdf_points, distance=2000,base_url='http://localhost:8080/ors'):
    """
    Compute walking isochrones around a set of point geometries using an OpenRouteService (ORS) instance.
    
    Each point in the input GeoDataFrame is used as the origin of an isochrone query. The function
    retrieves the resulting polygons (areas reachable within a specified distance) from the ORS API,
    merges them into a single dissolved polygon, and returns the resulting reachable area.

    Parameters
    ----------
    base_url : string
        url of local server
    gdf_points : geopandas.GeoDataFrame
        GeoDataFrame containing point geometries. Must be in a projected CRS (e.g., UTM) so that
        distances are in meters.
    distance : int or float, optional
        Maximum travel distance (in meters) for the isochrone. Default is 2000 (≈ 20 min walk at ~5 km/h).

    Returns
    -------
    geopandas.GeoDataFrame
        A GeoDataFrame containing a single (dissolved) polygon geometry that represents
        the union of all isochrone areas reachable within the specified distance.
    
    Notes
    -----
    - This function assumes that an OpenRouteService backend is running locally
      and accessible at `http://localhost:8080/ors`.
    - The ORS profile used is `"foot-walking"`. You can modify this to `"driving-car"`,
      `"cycling-regular"`, etc.
    - The function prints progress every 1000 points and logs any failed requests.
    - Large inputs (e.g., >10k points) may take significant time to process;
      consider batching or reducing input points before calling this function.
    
    Example
    -------
    >>> import geopandas as gpd
    >>> points = gpd.read_file("points.gpkg")
    >>> reachable_area = compute_isochrones(points, distance=2000)
    >>> reachable_area.plot()
    """
    
    # Initialize ORS client (local instance, no API key needed)
    client = openrouteservice.Client(base_url='http://localhost:8080/ors', key=None)

    isochrones = []
    
    for i, row in gdf_points.iterrows():
        coord = (row.geometry.x, row.geometry.y)
        try:
            # Request isochrone polygon for this point
            res = client.isochrones(
                locations=[coord],
                range_type='distance',
                profile="foot-walking",
                range=[distance],
            )
            poly = gpd.GeoDataFrame.from_features(res["features"])
            isochrones.append(poly)
        except Exception as e:
            print(f"Failed point {i}: {e}")
        
        # Print progress periodically
        if i % 1000 == 0 and i > 0:
            print(f"Processed {i} points...")

    if not isochrones:
        print("No isochrones generated — check ORS connection or input points.")
        return gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs=gdf_points.crs)

    # Combine all polygons into one GeoDataFrame
    reachable = gpd.GeoDataFrame(
        pd.concat(isochrones, ignore_index=True),
        geometry='geometry',
        crs=gdf_points.crs
    )

    # Dissolve to create a single merged area
    reachable_areas = reachable.dissolve()

    return reachable_areas

def map_roads(load_country):

    """ 
    To create a new column with an aggregated list of road types. 
    
    Args:
        *load_country* : A geodataframe containing all the roads of a country.
        
    Returns:
        *load_country* : The same geodataframe but with an additional 'roads' column containing the aggregated road types.
    """

    dict_map = {
"busway": "secondary",
"disused" : "other",
"dummy" : "other",
"planned" : "other",
"platform" : "other",
"unsurfaced" : "track",
"traffic_island" : "other",
"razed" : "other",
"abandoned" : "other",
"services" : "track",
"proposed" : "other",
"corridor" : "track",
"bus_guideway" : "other",
"bus_stop" : "other",
"rest_area" : "other",
"yes" : "other",
"trail" : "other",
"escape" : "track",
"raceway" : "other",
"emergency_access_point" : "track",
"emergency_bay" : "track",
"construction" : "track",
"bridleway" : "track",
"cycleway" : "other",
"footway" : "other",
"living_street" : "tertiary",
"path" : "track",
"pedestrian" : "other",
"primary" : "primary",
"primary_link" : "primary",
"residential" : "tertiary",
"road" : "secondary",
"secondary" : "secondary",
"secondary_link" : "secondary",
"service" : "tertiary",
"steps" : "other",
"tertiary" : "tertiary",
"tertiary_link" : "tertiary",
"track" : "track",
"unclassified" : "tertiary",
"trunk" : "primary",
"motorway" : "primary",
"trunk_link" : "primary",
"motorway_link" : "primary",
"via_ferrata": "other",
"elevator": "other",
"crossing": "other",
"seasonal": "other",
"traffic_signals":"other",
"piste":"other",
"dismantled": "other",
"winter_road":"other",
"access":"other",
"ohm:military:Trench":"other",
"no":"other",
"byway":"other",
"unmarked_route":"other",
"track_grade1":"track",
"track_grade2":"track",
"track_grade3":"track",
"track_grade4":"track",
"track_grade5":"track",
"unknown":"other"
}
    
    load_country['roads'] = load_country['fclass'].map(lambda x: (dict_map[x])) 
    
    return load_country

def create_midpoints(gdf_roads, spacing=2000):
    """
    Create one or more midpoint(s) along roads longer than the given spacing.
    
    Parameters
    ----------
    gdf_roads : geopandas.GeoDataFrame
        GeoDataFrame containing LineString or MultiLineString geometries (projected CRS).
    spacing : float, optional
        Minimum road length (in CRS units, e.g. meters) above which midpoints
        are generated. Points are spaced approximately every `spacing` meters.
        Default is 2000.
    
    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame containing midpoint geometries. Each row corresponds to
        a generated point, with attributes inherited from the source road.
    
    Example
    -------
    >>> roads = gpd.read_file("simplified_roads.gpkg")
    >>> midpoints = create_midpoints(roads, spacing=2000)
    >>> midpoints.plot()
    """
    points = []
    attrs = []

    for idx, row in gdf_roads.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        # Handle both LineString and MultiLineString
        lines = geom.geoms if geom.geom_type == "MultiLineString" else [geom]

        for line in lines:
            length = line.length
            if length <= spacing:
                # Skip short roads
                continue

            # Compute distances along the line where points will be placed
            n_segments = int(length // spacing)
            distances = np.linspace(spacing / 2, length - spacing / 2, n_segments)

            # Interpolate midpoints
            for d in distances:
                pt = line.interpolate(d)
                points.append(pt)
                attrs.append(row)

    midpoints_gdf = gpd.GeoDataFrame(attrs, geometry=points, crs=gdf_roads.crs)
    return midpoints_gdf

def grid_sample(points_gdf, cell_size=200):
    """
    Spatially downsample a point GeoDataFrame by retaining only one point per grid cell.

    This function overlays a regular grid of square cells (in the same CRS as the input)
    over the extent of the input points and keeps only one point per grid cell.
    It is useful for spatial thinning or for reducing dataset density while
    maintaining broad geographic coverage.

    Parameters
    ----------
    points_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing point geometries to be sampled.
    cell_size : float, optional
        The size of each grid cell, in the same units as the GeoDataFrame’s CRS
        (typically meters for projected coordinate systems). Default is 200.

    Returns
    -------
    geopandas.GeoDataFrame
        A subset of the input GeoDataFrame containing one representative point
        per grid cell.

    Notes
    -----
    - The function uses a spatial join to assign each point to a grid cell.
    - The first point encountered in each cell is kept; others are dropped.
    - The returned GeoDataFrame preserves the original CRS.

    Example
    -------
    >>> sampled = grid_sample(points_gdf, cell_size=500)
    >>> sampled.plot()
    """
    # Extract bounding box of all points
    minx, miny, maxx, maxy = points_gdf.total_bounds

    # Determine the number of grid cells along each axis
    nx = math.ceil((maxx - minx) / cell_size)
    ny = math.ceil((maxy - miny) / cell_size)

    # Build grid polygons
    cells = []
    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * cell_size
            y0 = miny + j * cell_size
            cells.append(box(x0, y0, x0 + cell_size, y0 + cell_size))

    grid = gpd.GeoDataFrame({"geometry": cells}, crs=points_gdf.crs)

    # Spatially join points to grid cells
    joined = gpd.sjoin(
        points_gdf,
        grid.reset_index().rename(columns={"index": "cell"}),
        how="left",
        predicate="within"
    )

    # Keep the first point per grid cell
    kept = joined.drop_duplicates(subset=["cell"]).drop(columns=["index_right", "cell"])

    return kept


def extract_vertices(roads_gdf):
    """
    Extract all vertex coordinates from LineString and MultiLineString geometries.

    This function takes a GeoDataFrame of road (or linear) geometries and returns
    a new GeoDataFrame containing all vertex coordinates as individual Point
    geometries.

    Parameters
    ----------
    roads_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing LineString or MultiLineString geometries.

    Returns
    -------
    geopandas.GeoDataFrame
        A GeoDataFrame of Point geometries representing all vertices from the
        input line features. The CRS of the input is preserved.

    Notes
    -----
    - Empty or null geometries are ignored.
    - Useful for network analysis, snapping operations, or spatial sampling
      along road networks.

    Example
    -------
    >>> vertices = extract_vertices(roads_gdf)
    >>> vertices.head()
    """
    all_points = []

    for geom in roads_gdf.geometry:
        if geom is None or geom.is_empty:
            continue

        # Handle single LineString
        if geom.geom_type == "LineString":
            all_points.extend([Point(xy) for xy in geom.coords])

        # Handle MultiLineString
        elif geom.geom_type == "MultiLineString":
            for part in geom.geoms:
                all_points.extend([Point(xy) for xy in part.coords])

    return gpd.GeoDataFrame(geometry=all_points, crs=roads_gdf.crs)


def sample_points_within_roads(df_roads, simplify_tolerance=2000, midpoint_spacing=1000, grid_cell_size=1000):
    """
    Generate representative sample points from a set of road geometries.

    This function simplifies road geometries, extracts their vertices,
    creates midpoints for longer roads, and performs spatial thinning using
    a regular grid. The final result is a GeoDataFrame of sampled points that
    efficiently represent the road network for spatial analyses (e.g. isochrones,
    accessibility modeling, or network coverage).

    Parameters
    ----------
    df_roads : geopandas.GeoDataFrame
        Input GeoDataFrame containing LineString or MultiLineString road geometries.
        Must be in a projected coordinate system (e.g., meters).
    simplify_tolerance : float, optional
        Distance tolerance (in CRS units) used to simplify geometries and remove
        small detail. Larger values produce simpler geometries and fewer vertices.
        Default is 2000 (≈ 2 km).
    midpoint_spacing : float, optional
        Minimum road length threshold (in CRS units) above which midpoints
        are generated along each road (via `create_midpoints`). Default is 1000 (≈ 1 km).
    grid_cell_size : float, optional
        Cell size (in CRS units) used to thin dense vertex points by keeping
        one point per grid cell (via `grid_sample`). Default is 1000 (≈ 1 km).

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame of sampled points in EPSG:4326 (WGS84), combining:
        - Simplified road vertices (thinned via grid)
        - Midpoints along long roads

    Notes
    -----
    - Requires helper functions:
        * `extract_vertices(df)` — extracts all vertex points from line geometries.
        * `create_midpoints(df, spacing)` — generates one or more midpoints for
          roads longer than a given spacing.
        * `grid_sample(df, cell_size)` — selects one representative point per grid cell.
    - Input geometries should use a projected CRS in meters for consistent distance units.
    - The resulting points are reprojected to EPSG:4326 for interoperability (e.g. with APIs).

    Example
    -------
    >>> roads = gpd.read_file("simplified_roads.gpkg")
    >>> sampled_points = sample_points_within_roads(roads)
    >>> sampled_points.to_file("road_sample_points.gpkg", driver="GPKG")
    """

    # 1️⃣ Simplify geometries to reduce vertex density
    df_roads_simp = df_roads.copy()
    df_roads_simp["geometry"] = df_roads_simp.simplify(
        tolerance=simplify_tolerance, preserve_topology=True
    )

    # 2️⃣ Extract all vertices from simplified geometries
    points_gdf = extract_vertices(df_roads_simp)

    # 3️⃣ Generate midpoints for roads longer than `midpoint_spacing`
    midpoints = create_midpoints(df_roads_simp, spacing=midpoint_spacing)

    # 4️⃣ Perform grid-based thinning on vertex points
    points_gdf_sample = grid_sample(points_gdf, cell_size=grid_cell_size)

    # 5️⃣ Combine thinned vertices + midpoints
    sample_points_ps = pd.concat(
        [points_gdf_sample, midpoints[["geometry"]]],
        ignore_index=True
    )

    # 6️⃣ Reproject to WGS84 for API compatibility or visualization
    sample_points_ps = gpd.GeoDataFrame(sample_points_ps, geometry="geometry", crs=df_roads.crs)
    sample_points_ps = sample_points_ps.to_crs(epsg=4326)

    return sample_points_ps

def downsample_raster(scaling_factor, input_path, output_path):
    """
    Downsample a raster to a coarser spatial resolution by aggregating source pixels.

    This function reduces the spatial resolution of a raster dataset by a given 
    scaling factor. It computes a new transform and grid size, then reprojects 
    the input raster into the coarser grid using pixel aggregation. By default, 
    it uses the `sum` resampling method, meaning that the values of all source 
    pixels contributing to a new (larger) pixel are summed.

    Parameters
    ----------
    scaling_factor : float
        Factor by which to reduce the raster resolution.
        For example, `10` converts a 100 m raster into 1 km resolution.
    input_path : str
        Path to the input raster file (e.g., GeoTIFF).
    output_path : str
        Path to save the downsampled output raster (GeoTIFF).

    Notes
    -----
    - The coordinate reference system (CRS) of the output raster is the same 
      as the input raster.
    - The output raster’s resolution is increased by `scaling_factor`, and the 
      number of pixels is correspondingly reduced.
    - The function uses `Resampling.sum` to aggregate pixel values; for mean 
      values instead, replace with `Resampling.average`.

    Example
    -------
    >>> downsample_raster(
    ...     scaling_factor=10,
    ...     input_path='Data/population_100m.tif',
    ...     output_path='Data/population_1km.tif'
    ... )
    This example aggregates the population raster from 100 m to 1 km cells.

    Returns
    -------
    None
        The function writes the downsampled raster to `output_path`.

    """
    with rasterio.open(input_path) as src:

        # Compute new transform and shape
        transform, width, height = calculate_default_transform(
            src.crs, src.crs,
            src.width, src.height,
            *src.bounds,
            resolution=(src.res[0] * scaling_factor, src.res[1] * scaling_factor)
        )

        profile = src.profile
        profile.update({
            'transform': transform,
            'width': width,
            'height': height,
            'dtype': 'float32'
        })

        data_downsample = np.empty((height, width), dtype=np.float32)

        reproject(
            source=rasterio.band(src, 1),
            destination=data_downsample,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=src.crs,
            resampling=Resampling.sum  # aggregate by sum
        )

        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data_downsample, 1)

def calculate_population(morocco_shp, missing_areas, data_dir="Data"):
    """
    Calculate total and missing population for Morocco from WorldPop 2023 data.

    Downloads population rasters automatically if they do not exist locally.

    Parameters
    ----------
    morocco_shp : str or geopandas.GeoDataFrame
        Path to Morocco shapefile or a GeoDataFrame.
    missing_areas : str or geopandas.GeoDataFrame
        Path to shapefile or GeoDataFrame of missing areas.
    data_dir : str, optional
        Directory where population data is stored or downloaded to. Default is 'Data'.

    Returns
    -------
    dict
        Dictionary containing total and missing population counts.
    """

    os.makedirs(data_dir, exist_ok=True)

    files = {
        "pop_north": {
            "url": "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/2023/MAR/v1/100m/constrained/mar_pop_2023_CN_100m_R2025A_v1.tif",
            "path": os.path.join(data_dir, "pop_north.tif")
        },
        "pop_south": {
            "url": "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/2023/ESH/v1/100m/constrained/esh_pop_2023_CN_100m_R2025A_v1.tif",
            "path": os.path.join(data_dir, "pop_south.tif")
        }
    }

    # Ensure files exist (download if missing)
    for name, info in files.items():
        if not os.path.exists(info["path"]):
            print(f"Downloading {name} data...")
            urllib.request.urlretrieve(info["url"], info["path"])
            print(f"{name} data saved to {info['path']}")

    # Compute population stats
    pop_north_total = sum(item['sum'] for item in zonal_stats(morocco_shp, files["pop_north"]["path"], stats="sum") if item['sum'] is not None)
    pop_south_total = sum(item['sum'] for item in zonal_stats(morocco_shp, files["pop_south"]["path"], stats="sum") if item['sum'] is not None)
    pop_total = pop_north_total + pop_south_total

    pop_north_missing = sum(item['sum'] for item in zonal_stats(missing_areas, files["pop_north"]["path"], stats="sum", all_touched=True) if item['sum'] is not None)
    pop_south_missing = sum(item['sum'] for item in zonal_stats(missing_areas, files["pop_south"]["path"], stats="sum", all_touched=True) if item['sum'] is not None)
    pop_missing = pop_north_missing + pop_south_missing

    return {
        "pop_total": pop_total,
        "pop_excluded": pop_missing
    }
    
def create_map(gdf, threshold, column_name,line_opacity=0.3):
    """
    Create an interactive Folium choropleth map showing spatial variation
    of the "Share_enclavée" indicator, filtered by a threshold.

    This function filters the input GeoDataFrame to include only features
    where the "Share_enclavée" value exceeds a specified threshold, then
    displays the result as a color-coded choropleth map using a 
    yellow-to-red ("YlOrRd") color scale. Each feature is interactive, 
    with tooltips displaying the feature name and its value.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the geometries and associated attributes.
        Must include a numeric column named "Share_enclavée" and a unique 
        identifier column specified by `column_name`.
    threshold : float
        Minimum value of "Share_enclavée" required for a feature to be 
        included on the map.
    column_name : str
        Name of the column in `gdf` to use as the feature identifier 
        (e.g., "Commune", "Province", "Region").
    line_opacity : float
        Parameter to handle line opacity of geometries

    Returns
    -------
    folium.Map
        A Folium map object displaying the filtered choropleth layer, 
        suitable for rendering in a Jupyter notebook or exporting to HTML.

    Notes
    -----
    - The map is centered over Morocco (latitude 29.17, longitude −8.79) 
      with a zoom level of 6 by default.
    - The legend for "Taux d'enclavement (%)" is positioned in the 
      upper-left corner.
    - Features are highlighted on hover, and tooltips show both the 
      identifier and "Share_enclavée" value.
    - The color scheme follows the "YlOrRd" (Yellow-Orange-Red) palette.

    Example
    -------
    >>> m = create_map(gdf_communes, threshold=10, column_name="Commune")
    >>> m.save("map_enclavement.html")

    """
    gdf_ = gdf.copy()
    gdf_ = gdf_[gdf_.Share_enclavée > threshold]

    # Center map (you can adapt coordinates)
    m = folium.Map(location=[29.16981, -8.79054], zoom_start=6, tiles='cartodbpositron')
    gdf_['Share_enclavée_%'] = gdf_['Share_enclavée'] * 100
    # Add choropleth layer
    choropleth = folium.Choropleth(
        geo_data=gdf_,
        data=gdf_,
        columns=[column_name, "Share_enclavée_%"],
        key_on=f"feature.properties.{column_name}",
        fill_color="YlOrRd",
        fill_opacity=0.9,
        line_opacity=line_opacity,
        highlight=True,
        legend_name="Taux d'enclavement (%)"
    ).add_to(m)
    legend_style = """
    <style>
        .legend {
            font-size: 16px !important;   /* 🔹 increase font size */
            font-weight: bold;            /* optional: make it bold */
        }
    </style>
    """
    m.get_root().html.add_child(folium.Element(legend_style))
        # Add interactive tooltips (clickable features)
    if "Multi_Pove" in gdf.columns.unique():
        folium.GeoJson(
            gdf_,
            style_function=lambda feature: {
                'fillColor': 'transparent',
                'color': 'transparent',
                'weight': 0,
                'fillOpacity': 1
            },
            tooltip=folium.GeoJsonTooltip(
                fields=[column_name, "Share_enclavée","Multi_Pove"],
                aliases=[f"{column_name} :", "Taux d'enclavement :","Pauvreté Multidimensionnelle :"],
                localize=True,
                sticky=False
            )
        ).add_to(m)
    else:
        folium.GeoJson(
            gdf_,
            style_function=lambda feature: {
                'fillColor': 'transparent',
                'color': 'transparent',
                'weight': 0,
                'fillOpacity': 1
            },
            tooltip=folium.GeoJsonTooltip(
                fields=[column_name, "Share_enclavée",],
                aliases=[f"{column_name} :", "Taux d'enclavement :"],
                localize=True,
                sticky=False
            )
        ).add_to(m)

    # Move legend to upper-left using custom CSS
    legend_js = """
    <script>
        var legend = document.querySelector('.legend');
        if (legend) {
            legend.style.left = '10px';
            legend.style.right = 'auto';
            legend.style.top = '10px';
            legend.style.bottom = 'auto';
        }
    </script>
    """
    m.get_root().html.add_child(folium.Element(legend_js))

    title_html = f'''
    <div style="
        position: fixed;
        top: 10px;
        left: 50%;
        transform: translateX(-50%);
        z-index: 9999;
        background: rgba(255,255,255,0.85);
        padding: 6px 12px;
        border-radius: 6px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.3);
        font-size: 20px;
        font-weight: 700;
    ">
    Taux d'enclavement (seuil &gt; {int(100*threshold)}%) par entité
    </div>
    '''
    m.get_root().html.add_child(folium.Element(title_html))
    return m


def compute_ratio_rasters(raster_a_path, raster_b_path, output_path):
    """
    Compute the element-wise ratio of two raster datasets and save the result.

    This function reads two raster files of the same dimensions and coordinate 
    reference system, computes the ratio of their pixel values (raster_a / raster_b), 
    handles division by zero by assigning NaN, and writes the resulting raster 
    to an output file.

    Parameters
    ----------
    raster_a_path : str
        File path to the numerator raster (raster A) in GeoTIFF or compatible format.
    raster_b_path : str
        File path to the denominator raster (raster B) in GeoTIFF or compatible format.
        Must have the same dimensions and CRS as raster A.
    output_path : str
        File path where the resulting ratio raster will be saved.

    Notes
    -----
    - Division by zero in the denominator raster is handled by replacing zero values with NaN.
    - The output raster retains the metadata (CRS, transform, etc.) of raster A.
    - The output raster is saved in float32 format with a single band.

    Example
    -------
    >>> compute_ratio_rasters(
    ...     raster_a_path='Data/pop_enclavee_10km.tif',
    ...     raster_b_path='Data/pop_morocco_10km.tif',
    ...     output_path='Data/ratio_enclavement.tif'
    ... )

    Returns
    -------
    None
        The function writes the ratio raster to the specified output path.
    """
    # Read rasters
    with rasterio.open(raster_a_path) as src_a, rasterio.open(raster_b_path) as src_b:
        a = src_a.read(1).astype('float32')
        b = src_b.read(1).astype('float32')
        profile = src_a.profile  # keep metadata

    # Avoid division by zero
    b[b == 0] = np.nan

    # Compute ratio
    ratio = a / b

    profile.update(dtype=rasterio.float32, count=1)

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(ratio.astype(np.float32), 1)


def create_population_map(out_path, data_dir="Data"):
    """
    Create a unified population raster for Morocco from WorldPop (because South of Morocco is not included in WorldPop).

    This function:
      1. Downloads population raster data for northern Morocco and Western Sahara
         (if not already present).
      2. Merges the two rasters into a single continuous mosaic.
      3. Saves the result as a GeoTIFF file to the specified output path.

    Parameters
    ----------
    out_path : str
        Path where the merged population GeoTIFF will be saved.
    data_dir : str, optional
        Directory to store and read population raster files (default is `"Data"`).

    Returns
    -------
    None
        The merged raster is saved to disk; nothing is returned.

    Notes
    -----
    - Source data is downloaded from WorldPop (2023, R2025A release).
    - Missing values (-99999) are converted to NaN.
    - Output GeoTIFF is compressed, tiled, and uses `float32` precision.

    Example
    -------
    >>> create_population_map("Data/morocco_population_2023.tif")
    """
    # Define URLs and local paths for both rasters
    files = {
        "pop_north": {
            "url": "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/2023/MAR/v1/100m/constrained/mar_pop_2023_CN_100m_R2025A_v1.tif",
            "path": os.path.join(data_dir, "pop_north.tif")
        },
        "pop_south": {
            "url": "https://data.worldpop.org/GIS/Population/Global_2015_2030/R2025A/2023/ESH/v1/100m/constrained/esh_pop_2023_CN_100m_R2025A_v1.tif",
            "path": os.path.join(data_dir, "pop_south.tif")
        }
    }

    # Ensure both population rasters are available locally (download if missing)
    for name, info in files.items():
        if not os.path.exists(info["path"]):
            print(f"Downloading {name} data...")
            urllib.request.urlretrieve(info["url"], info["path"])
            print(f"{name} data saved to {info['path']}")
    
    # Open both population rasters
    raster1 = rasterio.open(files["pop_north"]["path"])
    raster2 = rasterio.open(files["pop_south"]["path"])
    
    # Merge them into one continuous raster mosaic
    mosaic, out_transform = merge([raster1, raster2])

    # Replace placeholder values with NaN
    mosaic[mosaic == -99999] = np.nan
    
    # Copy metadata and update for output
    out_meta = raster1.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform,
        "compress": "lzw",
        "tiled": True,
        "BIGTIFF": "IF_SAFER",
        "dtype": "float32",
        "nodata": np.nan
    })
    
    # Write the merged raster to disk
    with rasterio.open(out_path, "w", **out_meta) as dest:
        dest.write(mosaic)