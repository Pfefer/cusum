from rtree import index
import multiprocessing as mp
from functools import partial

from osgeo import gdal
import numpy as np
from pyrasta.raster import Raster
from tqdm import tqdm

PG_DESCRIPTION = dict(single_change_detection_tcs="single change")


def cdtec_run_tcs(fhandle, sources, dates, c_level, max_samples=500,
                  nb_processes=mp.cpu_count(), window_size=100,
                  chunksize=1, data_type="float32", no_data=0,
                  progress_bar=True):
    """ Run detection change function

    Parameters
    ----------
    fhandle: function
        Change detection function to be applied
    sources: list[pyrasta.raster.RasterBase]
        list of rasters
    dates: list of numpy.ndarray
        list of dates corresponding to rasters
    c_level: list of two thresholds (confidence level)
        Confidence level (between 0 and 1)
    max_samples: int
        Max number of samples for bootstrapping
    nb_processes: int
        Number of processes for multiprocessing
    window_size: int or tuple(int, int)
        Size of window for raster calculation
    chunksize: int
        Chunk size for multiprocessing imap function
    data_type: str
        Data type for output raster
    no_data: int
        No data value
    progress_bar: bool
        if True, display progress bar

    Returns
    -------

    """
    dates = np.asarray(dates)
    nb_samples = min(max_samples, np.math.factorial(len(sources)))

    if progress_bar:
        desc = f"Compute {PG_DESCRIPTION[fhandle.__name__]} " \
               f"(CL=%s)" % ("%.2f" % c_level[0] + ' & ' + "%.2f" % c_level[1]).lstrip('0')
    else:
        desc = None

    change = Raster.raster_calculation(sources, partial(fhandle, dates=dates,
                                                            threshold_high=c_level[0],
                                                            threshold_low=c_level[1],
                                                            nb_samples=nb_samples,
                                                            no_data=no_data),
                                       output_type=gdal.GetDataTypeByName(data_type),
                                       window_size=window_size,
                                       nb_processes=nb_processes,
                                       chunksize=chunksize,
                                       description=desc)

    return change


def intersects_gpd(gdf1, gdf2):
    """ Does layer intersect with other layer ? (Not element wise)

    Parameters
    ----------
    gdf1: geopandas.GeoDataFrame
    gdf2: geopandas.GeoDataFrame

    Returns
    -------
    numpy.ndarray
        Array of boolean
    """

    is_intersecting = []
    rtree_idx = r_tree_idx(gdf2.geometry)
    for geom in tqdm(gdf1.geometry):
        is_intersecting.append(any(intersects(geom, gdf2.geometry, rtree_idx)))

    return np.asarray(is_intersecting)


def intersects(geometry, geometry_collection, r_tree=None):
    """ Return if geometry intersects with geometries of collection

    Use this function with large geometry collections
    :param geometry:
    :param geometry_collection:
    :param r_tree:
    :return: list of boolean of length = length(geometry_collection)
    """
    # Use Rtree to speed up !
    if r_tree is None:
        r_tree = r_tree_idx(geometry_collection)

    list_of_intersecting_features = list(r_tree.intersection(geometry.bounds))

    return [False if f not in list_of_intersecting_features else geometry.intersects(geometry_collection[f]) for f in
            range(len(geometry_collection))]


def is_iterable(iterable):
    """ Check if input is iterable

    :param iterable:
    :return:
    """

    try:
        iter(iterable)
        return True
    except TypeError:
        return False


def r_tree_idx(geometry_collection):
    """ Return Rtree spatial index of geometry collection

    :param geometry_collection: geometry collection (list, series, etc.)
    :return:
    """

    idx = index.Index()
    if not is_iterable(geometry_collection):
        geometry_collection = [geometry_collection]

    for i, geom in enumerate(geometry_collection):
        idx.insert(i, geom.bounds)

    return idx
