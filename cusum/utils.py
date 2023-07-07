import numpy as np

from rtree import index
from tqdm import tqdm


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
