from rtree import index
import multiprocessing as mp
from functools import partial

from osgeo import gdal
import numpy as np
from pyrasta.raster import Raster
from tqdm import tqdm
import os
import re
import shutil
from pyroSAR import identify, identify_many, ID
from pyroSAR.snap.auxil import parse_recipe, parse_node, gpt, groupbyWorkers, writer, \
    windows_fileprefix, orb_parametrize, geo_parametrize, \
    mli_parametrize, dem_parametrize
from spatialist import Vector, intersect
from spatialist.ancillary import dissolve
import logging

log = logging.getLogger(__name__)


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

    change = Raster.raster_calculation(sources, partial(fhandle, dates=dates, threshold_high=c_level[0],
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


def sub_parametrize(scene, geometry=None, offset=None, buffer=0.01, copyMetadata=True, **kwargs):
    """
    convenience function for parametrizing an `Subset` node.

    Parameters
    ----------
    scene: pyroSAR.drivers.ID
        The SAR scene to be processed
    geometry: dict or spatialist.vector.Vector or str or None
        A vector geometry for geographic subsetting (node parameter geoRegion):

         - :class:`~spatialist.vector.Vector`: a vector object in arbitrary CRS
         - :class:`str`: a name of a file that can be read with :class:`~spatialist.vector.Vector` in arbitrary CRS
         - :class:`dict`: a dictionary with keys `xmin`, `xmax`, `ymin`, `ymax` in LatLon coordinates
    offset: tuple or None
        a tuple with pixel coordinates as (left, right, top, bottom)
    buffer: int or float
        an additional buffer in degrees to add around the `geometry`
    copyMetadata: bool
        copy the metadata of the source product?
    kwargs
        further keyword arguments for node parametrization. Known options:

         - fullSwath
         - referenceBand
         - sourceBands
         - subSamplingX
         - subSamplingY
         - tiePointGrids

    Returns
    -------
    Node
        the Subset node object
    """
    subset = parse_node('Subset')
    if geometry:

        if isinstance(geometry, Vector):
            shp = geometry.clone()
        elif isinstance(geometry, str):
            shp = Vector(geometry)
        else:
            raise TypeError("argument 'geometry' must be either a dictionary, a Vector object or a filename.")
        # reproject the geometry to WGS 84 latlon
        shp.reproject(4326)
        # shp.close()

        with shp as bounds:
            inter = intersect(scene.bbox(), bounds)
            if not inter:
                raise RuntimeError('no bounding box intersection between shapefile and scene')
            inter.close()
            wkt = bounds.convert2wkt()[0]
        shp.close()
        subset.parameters['region'] = [0, 0, scene.samples, scene.lines]
        subset.parameters['geoRegion'] = wkt
    #######################
    # (optionally) configure Subset node for pixel offsets
    elif offset and not geometry:
        # left, right, top and bottom offset in pixels
        l, r, t, b = offset

        subset_values = [l, t, scene.samples - l - r, scene.lines - t - b]
        subset.parameters['region'] = subset_values
        subset.parameters['geoRegion'] = ''
    else:
        raise RuntimeError("one of 'geometry' and 'offset' must be set")

    subset.parameters['copyMetadata'] = copyMetadata
    for key, val in kwargs.items():
        subset.parameters[key] = val
    return subset


def geocode(infile, outdir, t_srs=4326, spacing=20, polarizations='all', shapefile=None, scaling='dB',
            geocoding_type='Range-Doppler', removeS1BorderNoise=True, removeS1BorderNoiseMethod='pyroSAR',
            removeS1ThermalNoise=True, offset=None, allow_RES_OSV=False, demName='SRTM 1Sec HGT',
            externalDEMFile=None, externalDEMNoDataValue=None, externalDEMApplyEGM=True, terrainFlattening=True,
            basename_extensions=None, test=False, export_extra=None, groupsize=1, cleanup=True, tmpdir=None,
            gpt_exceptions=None, gpt_args=None, returnWF=False, nodataValueAtSea=True,
            demResamplingMethod='BILINEAR_INTERPOLATION', imgResamplingMethod='BILINEAR_INTERPOLATION',
            alignToStandardGrid=False, standardGridOriginX=0, standardGridOriginY=0,
            speckleFilter=False, refarea='gamma0', clean_edges=False, clean_edges_npixels=1,
            rlks=None, azlks=None, dem_oversampling_multiple=2):
    """
    general function for geocoding of SAR backscatter images with SNAP.

    This function performs the following steps:

    - (if necessary) identify the SAR scene(s) passed via argument `infile` (:func:`pyroSAR.drivers.identify`)
    - (if necessary) create the directories defined via `outdir` and `tmpdir`
    - (if necessary) download Sentinel-1 OSV files
    - parse a SNAP workflow (:class:`pyroSAR.snap.auxil.Workflow`)
    - write the workflow to an XML file in `outdir`
    - execute the workflow (:func:`pyroSAR.snap.auxil.gpt`)

    Note
    ----
    The function may create workflows with multiple `Write` nodes. All nodes are parametrized to write data in ENVI
    format,
    in which case the node parameter `file` is going to be a directory. All nodes will use the same temporary directory,
    which will be created in `tmpdir`.
    Its name is created from the basename of the `infile` (:meth:`pyroSAR.drivers.ID.outname_base`)
    and a suffix identifying each processing node of the workflow (:meth:`pyroSAR.snap.auxil.Workflow.suffix`).

    For example: `S1A__IW___A_20180101T170648_NR_Orb_Cal_ML_TF_TC`.

    Parameters
    ----------
    infile: str or ~pyroSAR.drivers.ID or list
        The SAR scene(s) to be processed; multiple scenes are treated as consecutive acquisitions, which will be
        mosaicked with SNAP's SliceAssembly operator.
    outdir: str
        The directory to write the final files to.
    t_srs: int or str or osgeo.osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default: `4326 <https://spatialreference.org/ref/epsg/4326/>`_.
    spacing: int or float, optional
        The target pixel spacing in meters. Default is 20
    polarizations: list[str] or str
        The polarizations to be processed; can be a string for a single polarization, e.g. 'VV', or a list of several
        polarizations, e.g. ['VV', 'VH']. With the special value 'all' (default) all available polarizations are
        processed.
    shapefile: str or :py:class:`~spatialist.vector.Vector` or dict, optional
        A vector geometry for subsetting the SAR scene to a test site. Default is None.
    scaling: {'dB', 'db', 'linear'}, optional
        Should the output be in linear or decibel scaling? Default is 'dB'.
    geocoding_type: {'Range-Doppler', 'SAR simulation cross correlation'}, optional
        The type of geocoding applied; can be either 'Range-Doppler' (default) or 'SAR simulation cross correlation'
    removeS1BorderNoise: bool, optional
        Enables removal of S1 GRD border noise (default). Will be ignored if SLC scenes are processed.
    removeS1BorderNoiseMethod: str, optional
        The border noise removal method to be applied if `removeS1BorderNoise` is True.
        See :func:`pyroSAR.S1.removeGRDBorderNoise` for details. One of the following:

         - 'ESA': the pure implementation as described by ESA
         - 'pyroSAR': the ESA method plus the custom pyroSAR refinement (default)
    removeS1ThermalNoise: bool, optional
        Enables removal of S1 thermal noise (default).
    offset: tuple, optional
        A tuple defining offsets for left, right, top and bottom in pixels, e.g. (100, 100, 0, 0); this variable is
        overridden if a shapefile is defined. Default is None.
    allow_RES_OSV: bool
        (only applies to Sentinel-1) Also allow the less accurate RES orbit files to be used?
        The function first tries to download a POE file for the scene.
        If this fails and RES files are allowed, it will download the RES file.
        The selected OSV type is written to the workflow XML file.
        Processing is aborted if the correction fails (Apply-Orbit-File parameter continueOnFail set to false).
    demName: str
        The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when `externalDEMFile` is not None.
        Supported options:

         - ACE2_5Min
         - ACE30
         - ASTER 1sec GDEM
         - CDEM
         - Copernicus 30m Global DEM
         - Copernicus 90m Global DEM
         - GETASSE30
         - SRTM 1Sec Grid
         - SRTM 1Sec HGT
         - SRTM 3Sec
    externalDEMFile: str or None, optional
        The absolute path to an external DEM file. Default is None. Overrides `demName`.
    externalDEMNoDataValue: int, float or None, optional
        The no data value of the external DEM. If not specified (default) the function will try to read it from the
        specified external DEM.
    externalDEMApplyEGM: bool, optional
        Apply Earth Gravitational Model to external DEM? Default is True.
    terrainFlattening: bool
        Apply topographic normalization on the data?
    basename_extensions: list of str or None
        Names of additional parameters to append to the basename, e.g. ['orbitNumber_rel'].
    test: bool, optional
        If set to True the workflow xml file is only written and not executed. Default is False.
    export_extra: list or None
        A list of image file IDs to be exported to outdir. The following IDs are currently supported:

         - incidenceAngleFromEllipsoid
         - localIncidenceAngle
         - projectedLocalIncidenceAngle
         - DEM
         - layoverShadowMask
         - scatteringArea (requires ``terrainFlattening=True``)
         - gammaSigmaRatio (requires ``terrainFlattening=True`` and ``refarea=['sigma0', 'gamma0']``)
    groupsize: int
        The number of workers executed together in one gpt call.
    cleanup: bool
        Should all files written to the temporary directory during function execution be deleted after processing?
        Default is True.
    tmpdir: str or None
        Path of custom temporary directory, useful to separate output folder and temp folder. If `None`, the `outdir`
        location will be used. The created subdirectory will be deleted after processing if ``cleanup=True``.
    gpt_exceptions: dict or None
        A dictionary to override the configured GPT executable for certain operators;
        each (sub-)workflow containing this operator will be executed with the define executable;

         - e.g. ``{'Terrain-Flattening': '/home/user/snap/bin/gpt'}``
    gpt_args: list or None
        A list of additional arguments to be passed to the gpt call.

        - e.g. ``['-x', '-c', '2048M']`` for increased tile cache size and intermediate clearing
    returnWF: bool
        Return the full name of the written workflow XML file?
    nodataValueAtSea: bool
        Mask pixels acquired over sea? The sea mask depends on the selected DEM.
    demResamplingMethod: str
        One of the following:

         - 'NEAREST_NEIGHBOUR'
         - 'BILINEAR_INTERPOLATION'
         - 'CUBIC_CONVOLUTION'
         - 'BISINC_5_POINT_INTERPOLATION'
         - 'BISINC_11_POINT_INTERPOLATION'
         - 'BISINC_21_POINT_INTERPOLATION'
         - 'BICUBIC_INTERPOLATION'
    imgResamplingMethod: str
        The resampling method for geocoding the SAR image; the options are identical to demResamplingMethod.
    alignToStandardGrid: bool
        Align all processed images to a common grid?
    standardGridOriginX: int or float
        The x origin value for grid alignment
    standardGridOriginY: int or float
        The y origin value for grid alignment
    speckleFilter: str
        One of the following:

         - 'Boxcar'
         - 'Median'
         - 'Frost'
         - 'Gamma Map'
         - 'Refined Lee'
         - 'Lee'
         - 'Lee Sigma'
    refarea: str or list
        'sigma0', 'gamma0' or a list of both
    clean_edges: bool
        erode noisy image edges? See :func:`pyroSAR.snap.auxil.erode_edges`.
        Does not apply to layover-shadow mask.
    clean_edges_npixels: int
        the number of pixels to erode.
    rlks: int or None
        the number of range looks. If not None, overrides the computation done by function
        :func:`pyroSAR.ancillary.multilook_factors` based on the image pixel spacing and the target spacing.
    azlks: int or None
        the number of azimuth looks. Like `rlks`.
    dem_oversampling_multiple: int
        a factor to multiply the DEM oversampling factor computed by SNAP.
        Used only for terrain flattening.
        The SNAP default of 1 has been found to be insufficient with stripe
        artifacts remaining in the image.

    Returns
    -------
    str or None
        Either the name of the workflow file if ``returnWF == True`` or None otherwise


    .. figure:: figures/snap_geocode.svg
        :align: center

        Function geocode workflow diagram for processing Sentinel-1 scenes.
        Dashed lines depict optional steps. The output is sigma or gamma nought
        backscatter with ellipsoid or radiometric terrain correction (suffix elp/rtc)
        as well as several optional ancillary datasets (controlled via argument `export_extra`).

    Examples
    --------
    geocode a Sentinel-1 scene and export the local incidence angle map with it

    See Also
    --------
    :class:`pyroSAR.drivers.ID`,
    :class:`spatialist.vector.Vector`,
    :func:`spatialist.auxil.crsConvert()`
    """
    if clean_edges:
        try:
            import scipy
        except ImportError:
            raise RuntimeError('please install scipy to clean edges')

    if isinstance(infile, ID):
        id = infile
        ids = [id]
    elif isinstance(infile, str):
        id = identify(infile)
        ids = [id]
    elif isinstance(infile, list):
        ids = identify_many(infile, sortkey='start')
        id = ids[0]
    else:
        raise TypeError("'infile' must be of type str, list or pyroSAR.ID")

    if id.is_processed(outdir):
        log.info('scene {} already processed'.format(id.outname_base()))
        return

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    ############################################
    # general setup
    process_S1_SLC = False

    if id.sensor in ['ASAR', 'ERS1', 'ERS2']:
        formatName = 'ENVISAT'
    elif id.sensor in ['S1A', 'S1B']:
        if id.product == 'SLC':
            removeS1BorderNoise = False
            process_S1_SLC = True
        formatName = 'SENTINEL-1'
    else:
        raise RuntimeError('sensor not supported (yet)')

    # several options like resampling are modified globally for the whole workflow at the end of this function
    # this list gathers IDs of nodes for which this should not be done because they are configured individually
    resampling_exceptions = []
    ######################
    if isinstance(polarizations, str):
        if polarizations == 'all':
            polarizations = id.polarizations
        else:
            if polarizations in id.polarizations:
                polarizations = [polarizations]
            else:
                raise RuntimeError('polarization {} does not exists in the source product'.format(polarizations))
    elif isinstance(polarizations, list):
        polarizations = [x for x in polarizations if x in id.polarizations]
    else:
        raise RuntimeError('polarizations must be of type str or list')

    swaths = None
    if process_S1_SLC:
        if id.acquisition_mode == 'IW':
            swaths = ['IW1', 'IW2', 'IW3']
        elif id.acquisition_mode == 'EW':
            swaths = ['EW1', 'EW2', 'EW3', 'EW4', 'EW5']
        elif re.search('S[1-6]', id.acquisition_mode):
            pass
        else:
            raise RuntimeError('acquisition mode {} not supported'.format(id.acquisition_mode))

    bandnames = dict()
    bandnames['beta0'] = ['Beta0_' + x for x in polarizations]
    bandnames['gamma0'] = ['Gamma0_' + x for x in polarizations]
    bandnames['sigma0'] = ['Sigma0_' + x for x in polarizations]
    bandnames['int'] = ['Intensity_' + x for x in polarizations]
    ############################################
    ############################################
    # parse base workflow
    workflow = parse_recipe('blank')
    ############################################
    if not isinstance(infile, list):
        infile = [infile]

    last = None
    collect = []
    for i in range(0, len(infile)):
        ############################################
        # Read node configuration
        read = parse_node('Read')
        workflow.insert_node(read)
        read.parameters['file'] = ids[i].scene
        read.parameters['formatName'] = formatName
        last = read
        ############################################
        # Remove-GRD-Border-Noise node configuration
        if id.sensor in ['S1A', 'S1B'] and id.product == 'GRD' and removeS1BorderNoise:
            bn = parse_node('Remove-GRD-Border-Noise')
            workflow.insert_node(bn, before=last.id)
            bn.parameters['selectedPolarisations'] = polarizations
            last = bn
        ############################################
        # Calibration node configuration
        cal = parse_node('Calibration')
        workflow.insert_node(cal, before=last.id)
        cal.parameters['selectedPolarisations'] = polarizations
        if isinstance(refarea, str):
            refarea = [refarea]
        for item in refarea:
            if item not in ['sigma0', 'gamma0']:
                raise ValueError('unsupported value for refarea: {}'.format(item))
        if terrainFlattening:
            cal.parameters['outputBetaBand'] = True
            cal.parameters['outputSigmaBand'] = False
        else:
            for opt in refarea:
                cal.parameters['output{}Band'.format(opt[:-1].capitalize())] = True
        if id.sensor in ['ERS1', 'ERS2', 'ASAR']:
            cal.parameters['createBetaBand'] = True
        last = cal
        ############################################
        # ThermalNoiseRemoval node configuration
        if id.sensor in ['S1A', 'S1B'] and removeS1ThermalNoise:
            tn = parse_node('ThermalNoiseRemoval')
            workflow.insert_node(tn, before=last.id)
            tn.parameters['selectedPolarisations'] = polarizations
            last = tn
        collect.append(last.id)
    ############################################
    # SliceAssembly node configuration
    if len(collect) > 1:
        sliceAssembly = parse_node('SliceAssembly')
        sliceAssembly.parameters['selectedPolarisations'] = polarizations
        workflow.insert_node(sliceAssembly, before=collect)
        last = sliceAssembly
    ############################################
    # TOPSAR-Deburst node configuration
    if process_S1_SLC and swaths is not None:
        deb = parse_node('TOPSAR-Deburst')
        workflow.insert_node(deb, before=last.id)
        deb.parameters['selectedPolarisations'] = polarizations
        last = deb
    ############################################
    # Apply-Orbit-File node configuration
    orb = orb_parametrize(scene=id, formatName=formatName, allow_RES_OSV=allow_RES_OSV)
    workflow.insert_node(orb, before=last.id)
    last = orb
    ############################################
    # Subset node configuration
    if shapefile is not None or offset is not None:
        sub = sub_parametrize(scene=id, geometry=shapefile, offset=offset, buffer=0.01)
        workflow.insert_node(sub, before=last.id)
        last = sub
    ############################################
    # Multilook node configuration
    if id.sensor in ['ERS1', 'ERS2', 'ASAR']:
        bands = bandnames['beta0'] + bandnames['sigma0']
    else:
        bands = None
    ml = mli_parametrize(scene=id, spacing=spacing, rlks=rlks, azlks=azlks, sourceBands=bands)
    if ml is not None:
        workflow.insert_node(ml, before=last.id)
        last = ml
    ############################################
    # Terrain-Flattening node configuration
    tf = None
    if terrainFlattening:
        tf = parse_node('Terrain-Flattening')
        workflow.insert_node(tf, before=last.id)
        tf.parameters['sourceBands'] = bandnames['beta0']
        tf.parameters['oversamplingMultiple'] = dem_oversampling_multiple
        if 'reGridMethod' in tf.parameters.keys():
            if externalDEMFile is None:
                tf.parameters['reGridMethod'] = True
            else:
                tf.parameters['reGridMethod'] = False
        if 'sigma0' in refarea:
            try:
                tf.parameters['outputSigma0'] = True
            except KeyError:
                raise RuntimeError("The Terrain-Flattening node does not accept "
                                   "parameter 'outputSigma0'. Please update S1TBX.")
        last = tf
    ############################################
    # merge bands to pass them to Terrain-Correction
    bm_tc = None
    bands = dissolve([bandnames[opt] for opt in refarea])
    if len(refarea) > 1 and terrainFlattening and 'scatteringArea' in export_extra:
        bm_tc = parse_node('BandMerge')
        workflow.insert_node(bm_tc, before=[last.source, last.id])
        sources = bm_tc.source
        gamma_index = sources.index('Terrain-Flattening')
        sigma_index = abs(gamma_index - 1)
        s1_id = os.path.basename(os.path.splitext(id.scene)[0])
        bands_long = []
        for band in bands:
            comp = [band + '::']
            if shapefile is not None:
                comp.append('Subset_')
            comp.append(s1_id)
            if band.startswith('Gamma'):
                comp.append('_' + workflow.suffix(stop=sources[gamma_index]))
            else:
                comp.append('_' + workflow.suffix(stop=sources[sigma_index]))
            bands_long.append(''.join(comp))
        bm_tc.parameters['sourceBands'] = bands_long
        last = bm_tc
    ############################################
    # Speckle-Filter node configuration
    speckleFilter_options = ['Boxcar',
                             'Median',
                             'Frost',
                             'Gamma Map',
                             'Refined Lee',
                             'Lee',
                             'Lee Sigma']

    if speckleFilter:
        message = '{0} must be one of the following:\n- {1}'
        if speckleFilter not in speckleFilter_options:
            raise ValueError(message.format('speckleFilter', '\n- '.join(speckleFilter_options)))
        sf = parse_node('Speckle-Filter')
        workflow.insert_node(sf, before=last.id)
        sf.parameters['sourceBands'] = None
        sf.parameters['filter'] = speckleFilter
        last = sf
    ############################################
    # configuration of node sequence for specific geocoding approaches
    tc = geo_parametrize(spacing=spacing, t_srs=t_srs,
                         tc_method=geocoding_type, sourceBands=bands,
                         alignToStandardGrid=alignToStandardGrid,
                         standardGridOriginX=standardGridOriginX,
                         standardGridOriginY=standardGridOriginY)
    workflow.insert_node(tc, before=last.id)
    if isinstance(tc, list):
        last = tc = tc[-1]
    else:
        last = tc
    ############################################
    # (optionally) add node for conversion from linear to db scaling
    if scaling not in ['dB', 'db', 'linear']:
        raise RuntimeError('scaling must be  a string of either "dB", "db" or "linear"')

    if scaling in ['dB', 'db']:
        lin2db = parse_node('LinearToFromdB')
        workflow.insert_node(lin2db, before=last.id)
        lin2db.parameters['sourceBands'] = bands
        last = lin2db
    ############################################
    # parametrize write node
    # create a suffix for the output file to identify processing steps performed in the workflow
    suffix = workflow.suffix()
    if tmpdir is None:
        tmpdir = outdir
    basename = os.path.join(tmpdir, id.outname_base(basename_extensions))

    outname = basename + '_' + suffix

    write = parse_node('Write')
    workflow.insert_node(write, before=last.id)
    write.parameters['file'] = outname
    write.parameters['formatName'] = 'ENVI'
    ############################################
    ############################################
    if export_extra is not None:
        tc_options = ['incidenceAngleFromEllipsoid',
                      'localIncidenceAngle',
                      'projectedLocalIncidenceAngle',
                      'DEM',
                      'layoverShadowMask']
        tc_selection = []
        for item in export_extra:
            if item in tc_options:
                key = 'save{}{}'.format(item[0].upper(), item[1:])
                tc.parameters[key] = True
                if item == 'DEM':
                    tc_selection.append('elevation')
                else:
                    tc_selection.append(item)
            elif item == 'scatteringArea':
                if not terrainFlattening:
                    raise RuntimeError('scatteringArea can only be created if terrain flattening is performed')
                area_select = parse_node('BandSelect')
                workflow.insert_node(area_select, before=tf.source, resetSuccessorSource=False)
                area_select.parameters['sourceBands'] = bandnames['beta0']

                area_merge1 = parse_node('BandMerge')
                workflow.insert_node(area_merge1, before=[tf.id, area_select.id], resetSuccessorSource=False)

                math = parse_node('BandMaths')
                workflow.insert_node(math, before=area_merge1.id, resetSuccessorSource=False)

                pol = polarizations[0]  # the result will be the same for each polarization
                area = 'scatteringArea_{0}'.format(pol)
                expression = 'Beta0_{0} / Gamma0_{0}'.format(pol)

                math.parameters.clear_variables()
                exp = math.parameters['targetBands'][0]
                exp['name'] = area
                exp['type'] = 'float32'
                exp['expression'] = expression
                exp['noDataValue'] = 0.0

                if len(refarea) > 1:
                    bm_tc.source = bm_tc.source + [math.id]
                else:
                    bm_tc = parse_node('BandMerge')
                    workflow.insert_node(bm_tc, before=[tf.id, math.id], resetSuccessorSource=False)
                    tc.source = bm_tc.id

                # modify Terrain-Correction source bands
                tc_bands = tc.parameters['sourceBands'] + ',' + area
                tc.parameters['sourceBands'] = tc_bands

                # add scattering Area to list of band directly written from Terrain-Correction
                tc_selection.append(area)
            elif item == 'gammaSigmaRatio':
                if not terrainFlattening:
                    raise RuntimeError('gammaSigmaRatio can only be created if terrain flattening is performed')
                if sorted(refarea) != ['gamma0', 'sigma0']:
                    raise ValueError("For export_extra layer 'gammaSigmaRatio' 'refarea' "
                                     "must contain both sigma0 and gamma0")
                math = parse_node('BandMaths')
                workflow.insert_node(math, before=tf.id, resetSuccessorSource=False)

                pol = polarizations[0]  # the result will be the same for each polarization
                ratio = 'gammaSigmaRatio_{0}'.format(pol)
                expression = 'Sigma0_{0} / Gamma0_{0}'.format(pol)

                math.parameters.clear_variables()
                exp = math.parameters['targetBands'][0]
                exp['name'] = ratio
                exp['type'] = 'float32'
                exp['expression'] = expression
                exp['noDataValue'] = 0.0

                if len(refarea) > 1:
                    bm_tc.source = bm_tc.source + [math.id]
                else:
                    bm_tc = parse_node('BandMerge')
                    workflow.insert_node(bm_tc, before=[tf.id, math.id], resetSuccessorSource=False)
                    tc.source = bm_tc.id

                # modify Terrain-Correction source bands
                tc_bands = tc.parameters['sourceBands'] + ',' + ratio
                tc.parameters['sourceBands'] = tc_bands

                # add scattering Area to list of band directly written from Terrain-Correction
                tc_selection.append(ratio)
            else:
                raise RuntimeError("ID '{}' not valid for argument 'export_extra'".format(item))
        # directly write export_extra layers to avoid dB scaling
        if scaling in ['db', 'dB'] and len(tc_selection) > 0:
            tc_write = parse_node('Write')
            workflow.insert_node(tc_write, before=tc.id, resetSuccessorSource=False)
            tc_write.parameters['file'] = outname
            tc_write.parameters['formatName'] = 'ENVI'
            tc_select = parse_node('BandSelect')
            workflow.insert_node(tc_select, after=tc_write.id)
            tc_select.parameters['sourceBands'] = tc_selection
    ############################################
    ############################################
    # DEM handling
    dem_parametrize(workflow=workflow, demName=demName,
                    externalDEMFile=externalDEMFile,
                    externalDEMNoDataValue=externalDEMNoDataValue,
                    externalDEMApplyEGM=externalDEMApplyEGM)
    ############################################
    ############################################
    # configure the resampling methods

    options_img = ['NEAREST_NEIGHBOUR',
                   'BILINEAR_INTERPOLATION',
                   'CUBIC_CONVOLUTION',
                   'BISINC_5_POINT_INTERPOLATION',
                   'BISINC_11_POINT_INTERPOLATION',
                   'BISINC_21_POINT_INTERPOLATION',
                   'BICUBIC_INTERPOLATION']
    options_dem = options_img + ['DELAUNAY_INTERPOLATION']

    message = '{0} must be one of the following:\n- {1}'
    if demResamplingMethod not in options_dem:
        raise ValueError(message.format('demResamplingMethod', '\n- '.join(options_dem)))
    if imgResamplingMethod not in options_img:
        raise ValueError(message.format('imgResamplingMethod', '\n- '.join(options_img)))

    workflow.set_par('demResamplingMethod', demResamplingMethod)
    workflow.set_par('imgResamplingMethod', imgResamplingMethod,
                     exceptions=resampling_exceptions)
    ############################################
    ############################################
    # additional parameter settings applied to the whole workflow

    workflow.set_par('nodataValueAtSea', nodataValueAtSea)
    ############################################
    ############################################
    # write workflow to file and optionally execute it
    log.debug('writing workflow to file')

    wf_name = outname.replace(tmpdir, outdir) + '_proc.xml'
    workflow.write(wf_name)

    # execute the newly written workflow
    if not test:
        try:
            groups = groupbyWorkers(wf_name, groupsize)
            gpt(wf_name, groups=groups, cleanup=cleanup, tmpdir=outname,
                gpt_exceptions=gpt_exceptions, gpt_args=gpt_args,
                removeS1BorderNoiseMethod=removeS1BorderNoiseMethod)
            writer(xmlfile=wf_name, outdir=outdir, basename_extensions=basename_extensions,
                   clean_edges=clean_edges, clean_edges_npixels=clean_edges_npixels)
        except Exception as e:

            with open(wf_name.replace('_proc.xml', '_error.log'), 'w') as logfile:
                logfile.write(str(e))
        finally:
            if cleanup and os.path.isdir(outname):
                log.info('deleting temporary files')
                shutil.rmtree(outname, onerror=windows_fileprefix)
        log.info('done')
    if returnWF:
        return wf_name
