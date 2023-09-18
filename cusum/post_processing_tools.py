from osgeo import gdal, osr, ogr
from tqdm import tqdm
import geopandas as gpd
import numpy as np
import os
import re
import warnings
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import multiprocessing as mp
import math
import pandas as pd

from cusum.utils import intersects_gpd

warnings.filterwarnings("ignore", message="pandas.Float64Index")
warnings.filterwarnings("ignore", message="pandas.Int64Index")
warnings.filterwarnings("ignore", message="The frame.append method is deprecated")


def make_paths(path):
    directories = ["Post_raster", "Post_shp", "Final_results"]

    for directory in directories:
        dir_path = os.path.join(path, directory)
        try:
            os.makedirs(dir_path)
        except FileExistsError:
            pass


def save_as_tif(raster_reference, input_array, output_filepath):
    ds = raster_reference
    nb_arrays = 1 if len(input_array.shape) == 2 else input_array.shape[2]
    nb_row, nb_col = input_array.shape[:2]

    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_filepath, nb_col, nb_row, nb_arrays, gdal.GDT_Float32)

    for band in range(nb_arrays):
        cur_array = input_array[:, :, band] if nb_arrays > 1 else input_array
        dst_ds.GetRasterBand(band + 1).WriteArray(cur_array)
        dst_ds.GetRasterBand(band + 1).SetNoDataValue(0)

    dst_ds.SetGeoTransform(ds.GetGeoTransform())
    dst_ds.SetProjection(ds.GetProjection())

    ds = None
    dst_ds = None


def cross_nb(path, ras1, ras2, mnf_shp):
    threshold = re.search(r"_\d{2}_|_\d{3}_", ras1)[0]
    period = re.findall(r"\d{8}", ras1)[0] + '_' + re.findall(r"\d{8}", ras1)[1]
    method = re.search(r"multi|single_tcs|single|", ras1)[0]
    output_filename = method + '_' + period + threshold

    output_path = os.path.join(path, "Post_raster")
    output_file_nb = os.path.join(output_path, output_filename + "nb_crossP.tif")
    output_file_rmv_nb = os.path.join(output_path, output_filename + "rmv_nb_crossP.tif")

    if not os.path.exists(output_file_rmv_nb):
        ds1, ds2 = gdal.Open(os.path.join(path, ras1)), gdal.Open(os.path.join(path, ras2))
        nb_bands = min(ds2.RasterCount, ds1.RasterCount)
        arr_sum = np.zeros((ds2.RasterYSize, ds2.RasterXSize))
        for band in tqdm(range(nb_bands), position=0, desc='Intersecting bands of: ' + ras1 + '\t' + ras2):
            ds1_band = ds1.GetRasterBand(band + 1)
            ds2_band = ds2.GetRasterBand(band + 1)
            arr1, arr2 = ds1_band.ReadAsArray(), ds2_band.ReadAsArray()
            arr_sum += np.where((arr1 != 0) & (arr2 != 0), 1, 0)

        if not os.path.exists(output_path):
            os.makedirs(output_path)
        save_as_tif(ds1, arr_sum, os.path.join(output_path, output_file_nb))

        if mnf_shp is not None:
            tqdm.write("\n Masking forest on : " + output_file_nb)
            gdal_mask(os.path.join(path, mnf_shp), output_file_nb, output_file_rmv_nb)
    else:
        tqdm.write("Process already done, skipping rasters : " + ras1 + '\t' + ras2)

    return output_file_rmv_nb if mnf_shp is not None else output_file_nb


def nb_change(path, ras, period):
    if not os.path.exists(os.path.join(path, "Post_raster", period + "_nb_change_95_VH.tif")):
        ds = gdal.Open(os.path.join(path, ras))
        nb_bands = ds.RasterCount
        arr_nb = np.zeros((ds.RasterYSize, ds.RasterXSize))
        for band in range(nb_bands):
            ds_band = ds.GetRasterBand(band + 1)
            arr = ds_band.ReadAsArray()
            arr_nb += np.where(arr != 0, 1, 0)
        save_as_tif(ds, arr_nb, os.path.join(path, "Post_raster", period + "_nb_change_95_VH.tif"))
    return os.path.join(path, "Post_raster", period + "_nb_change_95_VH.tif")


def remove_high_nb_change(ras, ras95_VH, nb_max):

    ds95_vh, ds2 = gdal.Open(ras95_VH), gdal.Open(ras, gdal.GA_Update)
    band_95_vh = ds95_vh.GetRasterBand(1)
    arr95_VH = band_95_vh.ReadAsArray()

    band = ds2.GetRasterBand(1)
    arr = band.ReadAsArray()
    arr = np.where(arr95_VH > nb_max, 0, arr)
    arr = np.where(arr > 0, 1, 0)
    ds2.GetRasterBand(1).WriteArray(arr)


def polygonize_raster(path, filename):

    tqdm.write("\n Polygonize raster : " + filename)

    src_ds = gdal.Open(os.path.join(path, "Post_raster", filename))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())
    src_band = src_ds.GetRasterBand(1)

    dst_layer_name = "POLYGONIZED_STUFF"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(os.path.join(path, "Post_shp", filename.replace(".tif", ".shp")))
    dst_layer = dst_ds.CreateLayer(dst_layer_name, srs=srs)
    dn_field = ogr.FieldDefn('DN', ogr.OFTReal)

    dst_layer.CreateField(dn_field)
    gdal.Polygonize(src_band, src_band, dst_layer, 0, [], callback=None)
    dst_ds = None


def repair_shapefile(input_filepath):
    """
    Code to solve the self-intersection polygons of the polygonize_raster output

    Parameters
    ----------
    input_filepath : str
        Path to the .shp to repair the self intersections of.

    Returns
    -------
    None.

    """

    gk500 = ogr.Open(input_filepath, 1)
    gk_lyr = gk500.GetLayer()
    for feature in gk_lyr:

        geom = feature.GetGeometryRef()
        if not geom.IsValid():
            feature.SetGeometry(geom.Buffer(0))
            gk_lyr.SetFeature(feature)
            assert feature.GetGeometryRef().IsValid()

    gk_lyr.ResetReading()
    assert all(feature.GetGeometryRef().IsValid() for feature in gk_lyr)


def area_filter(geodataframe, area_threshold):
    """
    Code to filter a geodataframe according to the polygon area in square meters.

    Parameters
    ----------
    geodataframe : geopandas.geodataframe
        Geopandas geodataframe to remove of the polygons showing a lower area
        than the threshold.
    area_threshold : int
        Area threshold in square meters.

    Returns
    -------
    geodataframe
        geopandas geodataframe composed by polygons showing a higher area than the threshold.

    """
    temp = geodataframe.to_crs('EPSG:3857')
    mask = temp.area > area_threshold
    return geodataframe.loc[mask]


def method_cross(high_tc, low_tc, area_threshold_h, area_threshold_l, chunk_number=None):
    """
    Cross-Tc computation : Working on the 'intersection' rule. If a low Tc polygon shows
    AT LEAST one intersection with a high tc polygon, it is not removed.

    Parameters
    ----------
    chunk_number : str
        Identification of chunk_number to assert multiprocess correctly
    high_tc : geopandas file
        High Tc (100) geopandas opened file.
    low_tc : geopandas file
        Low Tc (75) geopandas opened file.
    area_threshold_h : int
        Threshold based on the area (square meters) to remove the polygons showing
        a lower area.
    area_threshold_l : int
        Threshold based on the area (square meters) to remove the polygons showing
        a lower area.

    Returns
    -------
    output : geopandas geodataframe
        DESCRIPTION.

    """

    crs_tc = high_tc.to_crs('EPSG:3857')
    crs_low = low_tc.to_crs('EPSG:3857')
    filtered_high_crs = area_filter(crs_tc, area_threshold_h).to_crs('EPSG:4326')

    filtered_low_crs = area_filter(crs_low, area_threshold_l).to_crs('EPSG:4326')
    if chunk_number is not None:
        temp_high, temp_low = 'temp_high' + chunk_number + '.shp', 'temp_low' + chunk_number + '.shp'
    else:
        temp_high, temp_low = 'temp_high.shp', 'temp_low.shp'

    for filtered_crs, temp_file in [(filtered_high_crs, temp_high), (filtered_low_crs, temp_low)]:

        if filtered_crs.empty:

            tiny_polygon = Polygon([(0, 0), (0.000001, 0), (0.000001, 0.000001), (0, 0.000001), (0, 0)])

            empty_gdf = gpd.GeoDataFrame({'geometry': [tiny_polygon]}, crs="EPSG:4326")

            empty_gdf.to_file(driver='ESRI Shapefile', filename=os.path.abspath(os.path.join(os.getcwd(), temp_file)))
        else:
            filtered_crs.to_file(driver='ESRI Shapefile',
                                 filename=os.path.abspath(os.path.join(os.getcwd(), temp_file)))
    geo_high = gpd.GeoDataFrame.from_file(os.path.abspath(os.path.join(os.getcwd(), temp_high)))
    geo_low = gpd.GeoDataFrame.from_file(os.path.abspath(os.path.join(os.getcwd(), temp_low)))
    bool_output = intersects_gpd(geo_low, geo_high)

    return geo_low[bool_output]


def merge_shapefiles(path, chunks_list_shapefiles, output_path):
    dfs = []

    for shapefile_path in chunks_list_shapefiles:
        df = gpd.read_file(os.path.join(path, shapefile_path))
        dfs.append(df)
    base_filename = os.path.splitext(chunks_list_shapefiles[0])[0]
    chunk_name = re.search(r'_chunk_(\d+)_', base_filename)[0]
    base_filename = base_filename.replace(chunk_name, "_")

    merged_gdf = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=True))

    merged_gdf.to_file(os.path.join(output_path, base_filename + ".shp"), driver='ESRI Shapefile')


def gdal_mask(kml_path, input_filepath, output_filepath):
    """

    Function to set to '0' all pixels not belonging to the polygons of the vector file..
    Parameters
    ----------
    kml_path : str
        path to your kml or .shp vector file.
    input_filepath : str
        path to your .tif image to mask pixels of.
    output_filepath :str

    Returns
    -------
    str
        path to the output image.

    """
    options = {'cropToCutline': False, 'dstNodata': 0,
               'cutlineDSName': kml_path,
               'callback_data': '.'}

    gdal.Warp(output_filepath, input_filepath, options=gdal.WarpOptions(**options))


def cross_tc(path, high_tc_file, low_tc_file, thresh_low, thresh_high, area_threshold_h, area_threshold_l):

    tqdm.write("\n Cross Tc on : " + high_tc_file + '\t' + low_tc_file)

    """
    Cross-Tc spatial recombination.

    Parameters
    ----------
    polygon_in_dir : str
        Path to the polygons (.shp).
    high_tc_file : str
        Path to the high Tc VECTOR file (.shp).
    low_tc_file : str
        Path to the low Tc VECTOR file (.shp).
    area_threshold_h : int
        Minimum Mapping Unit (square meters). Removes polygons < to the threshold
    area_threshold_l : int
        Minimum Mapping Unit (square meters). Removes polygons < to the threshold
    Returns
    -------
    None.

    """

    match_high = re.search(r"_\d_", high_tc_file)
    match_low = re.search(r"_\d_", low_tc_file)
    chunk_number = match_high[0] if match_high else match_low[0] if match_low else None

    try:
        high_tc_pd = gpd.read_file(os.path.join(path, high_tc_file))
        low_tc_pd = gpd.read_file(os.path.join(path, low_tc_file))
        if chunk_number:
            output_cross = method_cross(high_tc_pd, low_tc_pd, area_threshold_h, area_threshold_l, chunk_number)
        else:
            output_cross = method_cross(high_tc_pd, low_tc_pd, area_threshold_h, area_threshold_l)

        if not os.path.exists(os.path.join(path, "crossTc")):
            os.makedirs(os.path.join(path, "crossTc"))

        output_cross_path = os.path.join(path, "crossTc")

        output_cross.reset_index(drop=True, inplace=True)
        if not output_cross.empty:
            output_cross.to_file(os.path.join(
                output_cross_path, low_tc_file.replace(str(thresh_low), str(thresh_high) + '_' + str(thresh_low))),
                driver='ESRI Shapefile')

            repair_shapefile(os.path.join(output_cross_path, low_tc_file.replace(str(thresh_low),
                                                                                 str(thresh_high) + '_' + str(
                                                                                     thresh_low))))
        if output_cross.empty:
            print("Cross Tc between : " + high_tc_file + " and " + low_tc_file + " is empty so it has been skipped")

    except ValueError:
        pass


def parallel_cross_tc(paths, high_tc_files, low_tc_files, thresh_low, thresh_high, area_threshold_h, area_threshold_l,
                      cpu_cores):
    pool = mp.Pool(processes=cpu_cores)
    pool.starmap(cross_tc, [(paths, high_tc_file, low_tc_file, thresh_low, thresh_high, area_threshold_h,
                             area_threshold_l)
                            for path, high_tc_file, low_tc_file in zip(paths, high_tc_files, low_tc_files)])

    pool.close()
    pool.join()


def run_parallel_cross_tc(path, high_tc_files, low_tc_files, thresh_low, thresh_high, area_threshold_h,
                          area_threshold_l,
                          cpu_cores):

    parallel_cross_tc(path, high_tc_files, low_tc_files, thresh_low, thresh_high, area_threshold_h, area_threshold_l,
                      cpu_cores)


def erosion_dilation(path, cuts_algo, in_res, proj, date, buffer, th, tl, method):
    cuts_algo = gpd.read_file(os.path.join(in_res, cuts_algo))
    s = gpd.GeoSeries(cuts_algo.geometry)
    s = s.to_crs(proj)
    g = s.buffer(buffer, join_style=2)
    g = g.buffer(-buffer, join_style=2)
    g.to_crs(proj)

    g.to_file(driver='ESRI Shapefile',
              filename=os.path.abspath(os.path.join(path, "Final_results",
                                                    '_'.join([method, date, str(th), str(tl), '.shp']))))


def create_custom_grid(input_shapefile, output_folder, total_cells):
    gdf = gpd.read_file(input_shapefile)

    bounds = gdf.total_bounds
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    cell_size_x = width / math.floor(width * math.sqrt(total_cells / (width * height)))
    cell_size_y = height / math.floor(height * math.sqrt(total_cells / (width * height)))

    num_cols = math.floor(width / cell_size_x)
    num_rows = math.floor(height / cell_size_y)

    grid_polygons = []

    for row in range(num_rows):
        for col in range(num_cols):
            xmin = bounds[0] + col * cell_size_x
            xmax = xmin + cell_size_x
            ymin = bounds[1] + row * cell_size_y
            ymax = ymin + cell_size_y
            polygon = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])

            if any(polygon.intersects(geom) for geom in gdf.geometry):
                grid_polygons.append(polygon)

    grid_gdf = gpd.GeoDataFrame({'geometry': grid_polygons})
    grid_gdf.crs = gdf.crs

    grid_shapefile = os.path.join(output_folder, "grid.shp")
    grid_gdf.to_file(grid_shapefile, driver='ESRI Shapefile')

    return grid_shapefile


def spatially_split_shapefile(input_shapefile, grid_shapefile, output_folder, threshold):
    gdf = gpd.read_file(input_shapefile)
    gdf.crs = "EPSG:4326"

    grid_gdf = gpd.read_file(grid_shapefile)
    grid_gdf.crs = "EPSG:4326"

    input_filename = os.path.splitext(os.path.basename(input_shapefile))[0]

    for i, grid_polygon in enumerate(grid_gdf.geometry):
        selected_geometries = []

        for j, polygon in gdf.iterrows():
            if polygon.geometry.centroid.intersects(grid_polygon):
                selected_geometries.append(polygon.geometry)

        if selected_geometries:
            if len(selected_geometries) > 1:
                merged_geometry = cascaded_union(selected_geometries)
            else:
                merged_geometry = selected_geometries[0]

            selection = gpd.GeoDataFrame({'geometry': [merged_geometry]}, crs=f"EPSG:{4326}")
            selection = selection.explode()

            if not selection.empty:
                output_shapefile = os.path.join(output_folder, f"{input_filename}_chunk_{i}_{threshold}.shp")
                selection.to_file(output_shapefile, driver='ESRI Shapefile')
        else:
            empty_shapefile = os.path.join(output_folder, f"{input_filename}_chunk_{i}_{threshold}_empty.shp")
            open(empty_shapefile, 'w').close()
