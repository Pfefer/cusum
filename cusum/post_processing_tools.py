from osgeo import gdal, osr, ogr
from tqdm import tqdm
import geopandas as gpd
import numpy as np
import os
import re
import warnings

from cusum.utils import intersects_gpd

warnings.filterwarnings("ignore", message="pandas.Float64Index")
warnings.filterwarnings("ignore", message="pandas.Int64Index")
warnings.filterwarnings("ignore", message="The frame.append method is deprecated")


def make_paths(path):

    if not os.path.exists(os.path.join(path, "Post_raster")):
        os.makedirs(os.path.join(path, "Post_raster"))

    if not os.path.exists(os.path.join(path, "Post_shp")):
        os.makedirs(os.path.join(path, "Post_shp"))

    if not os.path.exists(os.path.join(path, "Final_results")):
        os.makedirs(os.path.join(path, "Final_results"))


def save_as_dualtif(raster_reference, vv_array, vh_array, output_filepath):
    ds = raster_reference
    geo_transform = ds.GetGeoTransform()
    wkt = ds.GetProjection()
    nb_row, nb_col = vv_array.shape

    driver = gdal.GetDriverByName("GTiff")  # Create GTiff image
    dst_ds = driver.Create(output_filepath, nb_col, nb_row, 2, gdal.GDT_Float32)

    array1 = np.array(vv_array)
    array2 = np.array(vh_array)
    # writing output raster
    dst_ds.GetRasterBand(1).WriteArray(array1)
    dst_ds.GetRasterBand(2).WriteArray(array2)
    # setting nodata value
    dst_ds.GetRasterBand(1).SetNoDataValue(0)
    dst_ds.GetRasterBand(2).SetNoDataValue(0)
    # setting extension of output raster
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform(geo_transform)
    # setting spatial reference of output raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dst_ds.SetProjection(srs.ExportToWkt())
    # Close output raster dataset

    ds = None
    dst_ds = None


def save_as_tif(raster_reference, input_array, output_filepath):
    """
    Parameters
    ----------
    raster_reference : Gdal image
        Gdal raster of reference, obtained using gdal.Open(file).
    input_array : list or 3d array
        List of arrays you want to save in the file (or 3d array).
    output_filepath : string
        Path of the file to save the arrays in.

    Returns
    -------
    None.

    """
    ds = raster_reference  # Get information from reference raster dataset
    geo_transform = ds.GetGeoTransform()
    wkt = ds.GetProjection()
    if isinstance(input_array, type(list(input_array))):  # Get the numbers of column, rows and arrays in the data
        nb_arrays = len(input_array)
        nb_row, nb_col = input_array[0].shape
    else:
        nb_dim = len(input_array.shape)
        if nb_dim > 2:
            nb_arrays = input_array.shape[2]
            nb_row, nb_col = input_array[:, :, 0].shape
        else:
            nb_arrays = 1
            nb_row, nb_col = input_array.shape

    driver = gdal.GetDriverByName("GTiff")  # Create GTiff image
    dst_ds = driver.Create(output_filepath, nb_col, nb_row, nb_arrays, gdal.GDT_Float32)
    if nb_arrays > 1:
        if isinstance(input_array, type(list(input_array))):
            for cur_band in range(len(input_array)):
                cur_array = input_array[cur_band]
                new_array = np.array(cur_array)
                dst_ds.GetRasterBand(cur_band + 1).WriteArray(new_array)  # Writing output raster
                dst_ds.GetRasterBand(cur_band + 1).SetNoDataValue(0)  # Setting nodata value
                dst_ds.SetGeoTransform(geo_transform)  # Setting extension of output raster
                srs = osr.SpatialReference()  # Setting spatial reference of output raster
                srs.ImportFromWkt(wkt)
                dst_ds.SetProjection(wkt)

        elif isinstance(input_array, type(np.array(input_array))):
            for logi in range(input_array.shape[2]):
                cur_array = input_array[:, :, logi]
                new_array = np.array(cur_array)

                # writing output raster
                dst_ds.GetRasterBand(logi + 1).WriteArray(new_array)
                # setting nodata value
                dst_ds.GetRasterBand(logi + 1).SetNoDataValue(0)
                # setting extension of output raster
                # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
                dst_ds.SetGeoTransform(geo_transform)
                # setting spatial reference of output raster
                srs = osr.SpatialReference()
                srs.ImportFromWkt(wkt)
                dst_ds.SetProjection(wkt)
                # Close output raster dataset
    else:
        new_array = np.array(input_array)

        # writing output raster
        dst_ds.GetRasterBand(1).WriteArray(new_array)
        # setting nodata value
        dst_ds.GetRasterBand(1).SetNoDataValue(0)
        # setting extension of output raster
        # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform(geo_transform)
        # setting spatial reference of output raster
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
        dst_ds.SetProjection(wkt)
    # Close output raster dataset
    ds = None
    dst_ds = None


def cross_nb(path, ras1, ras2, mnf_shp):
    threshold = re.search(r"_\d{2}_|_\d{3}_", ras1)[0]
    period = re.findall(r"\d{8}", ras1)[0] + '_' + re.findall(r"\d{8}", ras1)[1]
    method = re.search(r"multi|single", ras1)[0]

    if not os.path.exists(os.path.join(path, "Post_raster", method + '_' + period + threshold + "rmv_nb_crossP.tif")):

        ds1, ds2 = gdal.Open(os.path.join(path, ras1)), gdal.Open(os.path.join(path, ras2))
        nb_bands = min(ds2.RasterCount, ds1.RasterCount)
        arr_sum = np.zeros((ds2.RasterYSize, ds2.RasterXSize))
        for band in tqdm(range(nb_bands), position=0, desc='Intersecting bands of: ' + ras1 + '\t' + ras2):
            ds1_band = ds1.GetRasterBand(band + 1)
            ds2_band = ds2.GetRasterBand(band + 1)
            arr1, arr2 = ds1_band.ReadAsArray(), ds2_band.ReadAsArray()
            arr_sum += np.where((arr1 != 0) & (arr2 != 0), 1, 0)

        if not os.path.exists(os.path.join(path, "Post_raster")):
            os.makedirs(os.path.join(path, "Post_raster"))
        save_as_tif(ds1, arr_sum, os.path.join(path, "Post_raster", method + '_' + period + threshold
                                               + "nb_crossP.tif"))

        if mnf_shp is not None:
            tqdm.write("\n Masking forest on : " + method + '_' + period + threshold + "nb_crossP.tif")
            gdal_mask(os.path.join(path, mnf_shp),
                      os.path.join(path, "Post_raster", method + '_' + period + threshold + "nb_crossP.tif"),
                      os.path.join(path, "Post_raster", method + '_' + period + threshold + "rmv_nb_crossP.tif"))
    else:
        tqdm.write("Process already done, skipping rasters : " + ras1 + '\t' + ras2)

    if mnf_shp is not None:
        return os.path.join(path, "Post_raster", method + '_' + period + threshold + "rmv_nb_crossP.tif")
    else:
        return os.path.join(path, "Post_raster", method + '_' + period + threshold + "nb_crossP.tif")


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


def polygonize_raster_2(in_dir, filename, zone, output_filepath):

    src_ds = gdal.Open(os.path.join(in_dir, filename))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())
    src_band = src_ds.GetRasterBand(1)

    #  create output datasource
    #
    dst_layer_name = "POLYGONIZED_STUFF"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(output_filepath)
    dst_layer = dst_ds.CreateLayer(dst_layer_name, srs=srs)
    dn_field = ogr.FieldDefn('DN', ogr.OFTReal)

    dst_layer.CreateField(dn_field)
    gdal.Polygonize(src_band, src_band, dst_layer, 0, [], callback=None)
    dst_ds = None


def polygonize_raster_multi(in_dir, filename):

    output_filepath = os.path.join(in_dir, "polygons", filename.replace('tif', 'shp'))

    src_ds = gdal.Open(os.path.join(in_dir, filename))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjection())
    src_band = src_ds.GetRasterBand(1)

    #  create output datasource
    #
    dst_layer_name = "POLYGONIZED_STUFF"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(output_filepath)
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
    # print('Repairing ', os.path.basename(input_filepath))

    gk500 = ogr.Open(input_filepath, 1)  # <====== 1 = update mode
    gk_lyr = gk500.GetLayer()
    # print(len(gk_lyr), ' features to check')
    for feature in gk_lyr:
        # logi += 1
        # if len(gk_lyr)!=0:
        # if logi%int(len(gk_lyr)/10) == 0:
        # print( str(logi/len(gk_lyr)*100).split('.')[0], ' % PROGRESS')
        geom = feature.GetGeometryRef()
        if not geom.IsValid():
            feature.SetGeometry(geom.Buffer(0))  # <====== SetGeometry
            gk_lyr.SetFeature(feature)  # <====== SetFeature
            assert feature.GetGeometryRef().IsValid()  # Doesn't fail

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


def method_5(high_tc, low_tc, area_threshold_h, area_threshold_l):
    """
    Cross-Tc computation : Working on the 'intersection' rule. If a low Tc polygon shows
    AT LEAST one intersection with a high tc polygon, it is not removed.

    Parameters
    ----------
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
    filtered_high_tc = area_filter(crs_tc, area_threshold_h)
    filtered_low_tc = area_filter(crs_low, area_threshold_l)
    filtered_high_crs = filtered_high_tc.to_crs('EPSG:4326')
    filtered_low_crs = filtered_low_tc.to_crs('EPSG:4326')
    filtered_high_crs.to_file(driver='ESRI Shapefile',
                              filename=os.path.abspath(os.path.join(os.getcwd(), 'temp_high.shp')))
    filtered_low_crs.to_file(driver='ESRI Shapefile',
                             filename=os.path.abspath(os.path.join(os.getcwd(), 'temp_low.shp')))

    # geo_high = PolygonLayer(os.path.abspath(os.path.join(os.getcwd(), 'temp_high.shp')))
    # geo_low = PolygonLayer(os.path.abspath(os.path.join(os.getcwd(), 'temp_low.shp')))
    geo_high = gpd.GeoDataFrame.from_file(os.path.abspath(os.path.join(os.getcwd(), 'temp_high.shp')))
    geo_low = gpd.GeoDataFrame.from_file(os.path.abspath(os.path.join(os.getcwd(), 'temp_low.shp')))
    # bool_output = geo_low.intersects(geo_high)
    bool_output = intersects_gpd(geo_low, geo_high)
    output = geo_low[bool_output]
    os.remove(os.path.abspath(os.path.join(os.getcwd(), 'temp_high.shp')))
    os.remove(os.path.abspath(os.path.join(os.getcwd(), 'temp_low.shp')))
    return output


def progress_cb(complete, cb_data):
    """
    Emit progress report in numbers for 10% intervals and dots for 3%
    """
    if int(complete * 100) % 10 == 0:
        print(f'{complete * 100:.0f}', end='', flush=True)
    elif int(complete * 100) % 3 == 0:
        print(f'{cb_data}', end='', flush=True)


def gdal_mask(kml_path, input_filepath, output_filepath):
    """

    Function to set to '0' all pixels not belonging to the polygons of the vector file..
    Parameters
    ----------
    kml_path : str
        path to your kml or .shp vector file.
    input_filepath : str
        path to your .tif image to mask pixels of.F

    Returns
    -------
    str
        path to the output image.

    """
    options = {'cropToCutline': False, 'dstNodata': 0,
               'cutlineDSName': kml_path,
               'callback_data': '.'}

    gdal.Warp(output_filepath, input_filepath, options=gdal.WarpOptions(**options))
    # exit_code = gdal.Warp(output_filepath, input_filepath,cutlineDSName = kml_path, \
    # dstNodata=0, cropToCutline = False,options =gdal.WarpOptions(**options));


def gdal_mask2(kml_path, input_filepath, output_filepath):
    """

    Function to set to '0' all pixels not belonging to the polygons of the vector file..
    Parameters
    ----------
    kml_path : str
        path to your kml or .shp vector file.
    input_filepath : str
        path to your .tif image to mask pixels of.
    output_filepath : str
        path to the output file.

    Returns
    -------
    str
        path to the output image.

    """
    options = {'cropToCutline': False, 'dstNodata': 0,
               'cutlineDSName': kml_path, 'callback': progress_cb,
               'callback_data': '.'}

    gdal.Warp(output_filepath, input_filepath, options=gdal.WarpOptions(**options))

    # exit_code = gdal.Warp(output_filepath, input_filepath,cutlineDSName = kml_path, \
    # dstNodata=0, cropToCutline = False,options =gdal.WarpOptions(**options));


def cleanup(filepath, band):
    """
    Remove the file after getting the array.

    Parameters
    ----------
    filepath : TYPE
        DESCRIPTION.
    band : TYPE
        DESCRIPTION.

    Returns
    -------
    img_arr : TYPE
        DESCRIPTION.

    """
    img = gdal.Open(filepath)
    img_band = img.GetRasterBand(band)
    img_arr = img_band.ReadAsArray()
    os.remove(filepath)
    return img_arr


def conversion(filelist_low, cross_list, in_dir, temp_folder, thresh):
    """
    Conversion of a .shp file to a .tif file using gdal_mask to remove pixels of
    the .tif file outside the polygons of the .shp file.

    Parameters
    ----------
    filelist_low : list
        List of the names of the .tif files.
    cross_list : list
        List of the cross-tc shapefiles. VECTOR FILES.
    in_dir : str
        path to the shapefile folder
    temp_folder : TYPE
        DESCRIPTION.
    thresh : int
        Tc threshold.

    Returns
    -------
    None.

    """
    for tif_low in tqdm(filelist_low):
        arr_list = list()
        img_ref = gdal.Open(os.path.join(in_dir, tif_low))
        nb_bands = img_ref.RasterCount
        for cross_shp, band in zip(cross_list, range(nb_bands)):
            file_temp = os.path.join(temp_folder, os.path.basename(cross_shp).replace('.shp', '.tif'))
            gdal_mask(cross_shp, os.path.join(in_dir, tif_low), file_temp)
            arr_list.append(cleanup(file_temp, band + 1))
        save_as_tif(img_ref, arr_list,
                    os.path.join(in_dir, tif_low.replace(str(thresh), '100_' + str(thresh))))


def cross_tc(path, high_tc_file, low_tc_file, thresh, thresh_high, area_threshold_h, area_threshold_l):

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

    high_tc_pd = gpd.read_file(os.path.join(path, "Post_shp", high_tc_file))
    low_tc_pd = gpd.read_file(os.path.join(path, "Post_shp", low_tc_file))

    output_cross = method_5(high_tc_pd, low_tc_pd, area_threshold_h, area_threshold_l)
    if not os.path.exists(os.path.join(path, "Post_shp", "crossTc")):
        os.makedirs(os.path.join(path, "Post_shp", "crossTc"))

    output_cross_path = os.path.join(path, "Post_shp", "crossTc")

    output_cross.to_file(os.path.join(
        output_cross_path, low_tc_file.replace(str(thresh), str(thresh_high) + '_' + str(thresh))),
        driver='ESRI Shapefile')

    repair_shapefile(os.path.join(output_cross_path, low_tc_file.replace(str(thresh), str(thresh_high) + '_' + str(thresh))))


def erosion_dilation(path, cuts_algo, in_res, proj, date, buffer, th, tl, method):

    cuts_algo = gpd.read_file(os.path.join(in_res, cuts_algo))
    s = gpd.GeoSeries(cuts_algo.geometry)
    s = s.to_crs("epsg:31979")
    g = s.buffer(buffer, join_style=2)
    g = g.buffer(-buffer, join_style=2)
    g.to_crs(proj)

    g.to_file(driver='ESRI Shapefile',
              filename=os.path.abspath(os.path.join(path, "Final_results", '_'.join([method, date, str(th), str(tl), '.shp']))))


# def intersection_with_plots(path, period, th, tl, shp_plots, zone):
#     if not os.path.exists(os.path.join(path, "results")):
#         os.makedirs(os.path.join(path, "results"))
#     in_inter = os.path.join(path, "results")
#
#     list = [f for f in os.listdir(os.path.join(path, "Post_shp", "crossTc")) if
#             (period in f) & ('.shp' in f) & (th + '_' + tl in f)]
#     cuts_algo = list[0]
#     geo_cuts_algo = gpd.read_file(os.path.join(path, "Post_shp", "crossTc", cuts_algo))
#     geo_plots = gpd.read_file(os.path.join(path, shp_plots))
#
#     s = geo_cuts_algo.to_crs('EPSG:3857')
#     t = geo_plots.to_crs('EPSG:3857')
#
#     s.to_file(driver='ESRI Shapefile',
#               filename=os.path.join(in_inter, "temp_cuts.shp"))
#     t.to_file(driver='ESRI Shapefile',
#               filename=os.path.join(in_inter, "temp_plots.shp"))
#
#     geo_cuts_algo = PolygonLayer(os.path.join(in_inter, "temp_cuts.shp"))
#     geo_plots = PolygonLayer(os.path.join(in_inter, "temp_plots.shp"))
#
#     geo_output = geo_cuts_algo.overlay(geo_plots, how='intersection')
#
#     geo_output.to_file(driver='ESRI Shapefile',
#                        file_path=os.path.join(in_inter, 'cuts_algo' + period + '_' + zone + '_' + tl + '.shp'))
#
#     list_to_be_removed = [f for f in os.listdir(in_inter) if ('temp' in f)]
#
#     for file in list_to_be_removed:
#         os.remove(os.path.join(in_inter, file))
