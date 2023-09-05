import asf_search as asf
import shutil
import pyroSAR.snap
from pyrasta.raster import Raster
from cdtec.base import single_change_detection, multi_change_detection, single_change_detection_tcs
from cdtec.main import cdtec_run
from cusum.utils import cdtec_run_tcs
from cusum.post_processing_tools import *


class CuSum_:

    def __init__(self, path, date_start, date_end, study_site):
        """

        Parameters
        ----------
        path: str
            Path to [...]
        date_start: str
        date_end: str
        study_site: str
        """

        self.path = path

        try:
            os.mkdir(os.path.join(self.path, "Downloaded_img"))
        except FileExistsError:
            pass
        finally:
            self.out_dl = os.path.join(self.path, "Downloaded_img")

        try:
            os.mkdir(os.path.join(self.path, "zone_files"))
        except FileExistsError:
            pass
        finally:
            self.out_zone = os.path.join(self.path, "zone_files")

        if '.shp' not in study_site:
            self.wkt = study_site
            zone_shp = gpd.GeoSeries.from_wkt([self.wkt], crs="EPSG:4326")
            zone_shp.to_file(os.path.join(self.out_zone, "zone.shp"))

            zone_filename = os.path.join(self.out_zone, "zone.shp")
            self.zone_files = zone_filename

        if '.shp' in study_site:
            geo = gpd.read_file(study_site).to_crs("EPSG:4326")
            wkt_zone = gpd.GeoDataFrame.to_wkt(geo)
            self.wkt = str(list(wkt_zone.geometry)[0])
            self.zone_files = study_site

        try:
            os.mkdir(os.path.join(path, "Pre_processed_img"))
        except FileExistsError:
            pass
        finally:
            self.out_pre = os.path.join(path, "Pre_processed_img")

        date1 = "".join(date_start.split('-'))
        date2 = "".join(date_end.split('-'))
        try:
            os.mkdir(os.path.join(path, "Output_CuSum" + "_" + date1 + "_" + date2))
        except FileExistsError:
            pass
        finally:
            self.out_algo = os.path.join(path, "Output_CuSum" + "_" + date1 + "_" + date2)

    def sentinel1_download(self, date_start, date_end, user, pw, flight_direction, platform, relativeOrbit=None,
                           mode="IW", level="GRD_HD", polarizations=None):

        """ Download S1 data

                                Description
                                -----------
                                Download S1 data between your start date and end date

                                Parameters
                                ----------
                                date_start : str
                                    First date of your time period
                                date_end : str
                                    Last date of your time period
                                user : str
                                    ASF username
                                pw : str
                                    ASF password
                                flight_direction : str
                                    Satellite orbit direction during data acquisition
                                relativeOrbit : int
                                    Path or track of satellite during data acquisition. For UAVSAR it is the Line ID
                                mode : str
                                    The mode used to acquire the data.
                                level : str
                                    Level to which the data has been processed
                                polarizations : str
                                    A property of SAR electromagnetic waves that can be used to extract meaningful
                                    information about surface properties of the earth.
                                platform : str
                                    Sentinel-1A or Sentinel-1B
                                    

                                """

        results = asf.search(start=date_start, end=date_end, processingLevel=level, flightDirection=flight_direction,
                             beamMode=mode, maxResults=200, intersectsWith=self.wkt, relativeOrbit=relativeOrbit,
                             polarization=polarizations, platform=platform)
        print(results)
        session = asf.ASFSession().auth_with_creds(user, pw)
        results.download(path=self.out_dl, session=session, processes=8)

    def preprocess(self, spacing=10, scaling='dB', t_srs=4326, allow_RES_OSV=True, polarizations='all',
                   removeS1BorderNoiseMethod='ESA', removeS1ThermalNoise=True,
                   demResamplingMethod='BILINEAR_INTERPOLATION',
                   imgResamplingMethod='BILINEAR_INTERPOLATION', DEM='SRTM 1Sec HGT', speckleFilter='Lee Sigma'):
        """ preprocess

                        Description
                        -----------
                        Preprocessing calculation

                        Parameters
                        ----------
                        spacing: int or float, optional
                            The target pixel spacing in meters. Default is 20

                        scaling : {'dB', 'db', 'linear'}, optional
                            Should the output be in linear or decibel scaling? Default is 'dB'.

                        allow_RES_OSV: bool

                            (only applies to Sentinel-1) Also allow the less accurate RES orbit files to be used?
                            The function first tries to download a POE file for the scene.
                            If these fails and RES files are allowed, it will download the RES file.
                            The selected OSV type is written to the workflow XML file.
                            Processing is aborted if the correction fails (Apply-Orbit-File parameter continueOnFail set
                            to false).

                        removeS1BorderNoiseMethod : The border noise removal method to be applied if
                            `removeS1BorderNoise` is True. See :func:`pyroSAR.S1.removeGRDBorderNoise` for details.
                            One of the following: - 'ESA': the pure implementation as described by ESA
                            - 'pyroSAR': the ESA method plus the custom pyroSAR refinement (default)

                        t_srs: int or str or osgeo.osr.SpatialReference
                            A target geographic reference system in WKT, EPSG, PROJ4 or
                            OPENGIS format. See function :func:`spatialist.auxil.crsConvert()` for details. Default:
                            `4326 <https://spatialreference.org/ref/epsg/4326/>`_.

                        scaling: {'dB', 'db', 'linear'}, optional
                            Should the output be in linear or decibel scaling? Default is 'dB'.

                        polarizations: list[str] or str
                            The polarizations to be processed; can be a string for a single polarization, e.g. 'VV',
                            or a list of several polarizations, e.g. ['VV', 'VH']. With the special value 'all'
                            (default) all available polarizations are processed.

                        removeS1ThermalNoise: bool, optional
                            Enables removal of S1 thermal noise (default).

                        removeS1BorderNoiseMethod: str, optional
                            The border noise removal method to be applied if `removeS1BorderNoise` is True.
                            See :func:`pyroSAR.S1.removeGRDBorderNoise` for details. One of the following:

                            - 'ESA': the pure implementation as described by ESA
                            - 'pyroSAR': the ESA method plus the custom pyroSAR refinement (default)

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
                            The resampling method for geocoding the SAR image; the options are identical to
                            demResamplingMethod.

                        DEM: str
                            The name of the auto-download DEM. Default is 'SRTM 1Sec HGT'. Is ignored when
                             `externalDEMFile` is not None.
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

                        speckleFilter: str
                            One of the following:
                            - 'Boxcar'
                            - 'Median'
                            - 'Frost'
                            - 'Gamma Map'
                            - 'Refined Lee'
                            - 'Lee'
                            - 'Lee Sigma'

                        """

        for filename in tqdm(os.listdir(self.out_dl), desc="Pre_processing"):
            pyroSAR.snap.util.geocode(os.path.join(self.out_dl, filename), self.out_pre, t_srs,
                                      spacing=spacing, polarizations=polarizations, shapefile=self.zone_files,
                                      scaling=scaling, removeS1BorderNoise=True, allow_RES_OSV=allow_RES_OSV,
                                      removeS1BorderNoiseMethod=removeS1BorderNoiseMethod,
                                      removeS1ThermalNoise=removeS1ThermalNoise, speckleFilter=speckleFilter,
                                      terrainFlattening=True, export_extra=None, groupsize=8, cleanup=True,
                                      gpt_args=None, returnWF=False, demResamplingMethod=demResamplingMethod,
                                      demName=DEM, imgResamplingMethod=imgResamplingMethod, alignToStandardGrid=False,
                                      refarea='gamma0', clean_edges=False, clean_edges_npixels=1)

    def run(self, c_levels, max_samples=500, nb_cores=8, method='single'):
        """ Raster expression calculation

                Description
                -----------
                Run CuSum algorithm

                Parameters
                ----------
                c_levels: List or float
                    List of thresholds

                max_samples: int
                    Maximum samples in bootstrap analysis, 500 is enough
                nb_cores: int
                    Numbers of CPU used to compute changes in images
                method: str
                    One of the following : 'single', 'single_tcs', 'multi'

                    Single computes changes occurred regardless to the number of them on each pixel. Multi compute the
                    number of changes on each pixel. Single_tcs computes changes as single but for two different
                    thresholds in the same time
                """

        if not os.path.exists(os.path.join(self.out_pre, "VV")):
            os.mkdir(os.path.join(self.out_pre, "VV"))
            in_vv = os.path.join(self.out_pre, "VV")

        else:
            in_vv = os.path.join(self.out_pre, "VV")

        if not os.path.exists(os.path.join(self.out_pre, "VH")):
            os.mkdir(os.path.join(self.out_pre, "VH"))
            in_vh = os.path.join(self.out_pre, "VH")

        else:
            in_vh = os.path.join(self.out_pre, "VH")

        for filename in os.listdir(self.out_pre):
            if ('VH' in filename) & ('.tif' in filename):
                shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.out_pre, "VH", filename))
            if ('VV' in filename) & ('.tif' in filename):
                shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.out_pre, "VV", filename))

        def clip():

            x_min_list = []
            x_max_list = []
            y_min_list = []
            y_max_list = []
            for file_vv, file_vh in zip(os.listdir(in_vv), os.listdir(in_vh)):
                raster_vv, raster_vh = Raster(os.path.join(in_vv, file_vv)), Raster(os.path.join(in_vh, file_vh))
                x_min, y_max, x_max, y_min = raster_vv.bounds
                x_min_list.append(x_min)
                x_max_list.append(x_max)
                y_min_list.append(y_min)
                y_max_list.append(y_max)

            for file_vv, file_vh in zip(os.listdir(in_vv), os.listdir(in_vh)):
                raster_vv, raster_vh = Raster(os.path.join(in_vv, file_vv)), Raster(os.path.join(in_vh, file_vh))
                new_vv = raster_vv.clip(bounds=(np.max(x_min_list), np.min(y_max_list),
                                                np.min(x_max_list), np.max(y_min_list)))
                new_vh = raster_vh.clip(bounds=(np.max(x_min_list), np.min(y_max_list),
                                                np.min(x_max_list), np.max(y_min_list)))
                new_vv.to_file(os.path.join(in_vv, file_vv))
                new_vh.to_file(os.path.join(in_vh, file_vh))

        clip()

        FHANDLE = dict(single=single_change_detection,
                       multi=multi_change_detection,
                       single_tcs=single_change_detection_tcs)

        for in_path in [in_vh, in_vv]:
            images = []
            dates = []

            pol = re.search("VV|VH", in_path)[0]

            for _, _, file in os.walk(in_path):

                for name in file:
                    if not name.startswith("."):
                        if '.tif' in name:
                            images.append(os.path.join(in_path, name))

                            dates.append(np.array(int(re.search(r"\d{8}", name)[0])))

            images.sort()

            sources = [Raster(img) for img in images]
            if method == 'single_tcs':
                result = cdtec_run_tcs(FHANDLE[method], sources, dates, c_levels, max_samples=max_samples,
                                       nb_processes=nb_cores, window_size=50, chunksize=1, data_type="float32",
                                       no_data=0,
                                       progress_bar=True)

                result1 = result.extract_bands([1])
                result2 = result.extract_bands([2])

                result2.to_file(os.path.join(self.out_algo, str(method) + "_" + str(dates[0]) + r"_" + str(dates[-1])
                                             + '_' + pol + "_" + str(int(c_levels[0] * 100)) + '_.tif'))
                result1.to_file(os.path.join(self.out_algo, str(method) + "_" + str(dates[0]) + r"_" + str(dates[-1])
                                             + '_' + pol + "_" + str(int(c_levels[1] * 100)) + '_.tif'))
            else:
                for c_level in c_levels:
                    result = cdtec_run(FHANDLE[method], sources, dates, c_level, max_samples=max_samples,
                                       nb_processes=nb_cores, window_size=50, chunksize=1, data_type="float32",
                                       no_data=0, progress_bar=True)

                    result.to_file(os.path.join(self.out_algo, str(method) + "_" + str(dates[0]) + r"_" + str(dates[-1])
                                                + '_' + pol + "_" + str(int(c_level * 100)) + '_.tif'))

    def post_processing(self, th, tl, method, mf_shp=None, area_th=300, area_tl=1000, nb_max=10):

        """ Post processing

                        Description
                        -----------
                        Post processing included :
                            Intersection of VV and VH changes,
                            Intersection with forest mask,
                            Remove of high number of changes (in case of 'multi' run)
                            CrossTc computation


                        Parameters
                        ----------
                        th: float
                            Threshold high
                        tl: float
                            Threshold low
                        method: str
                            One of the following : 'single', single_tcs, 'multi'
                            Single computes changes occurred regardless to the number of them on each pixel
                            Multi compute the number of changes on each pixel
                            Single_tcs computes changes as single but for two different thresholds in the same time

                        mf_shp: str
                            Filename of forest mask shapefile to be removed
                            Set up None if you do not need to mask

                        area_th: int
                            Minimum mapping in square meters for threshold high changes
                        area_tl: int
                            Minimum mapping in quare meters for threshold low changes (final minimum mapping)
                        nb_max: int
                            Maximum number of changes accepted

                        """
        th = int(th * 100)
        tl = int(tl * 100)

        filename = [file for file in os.listdir(self.out_algo) if (str(th) in file) & ('VH' in file)
                    & (method in file)][0]
        period = re.search(r"\d{8}_\d{8}", filename)[0]

        make_paths(self.out_algo)

        high_vh_list = [file for file in os.listdir(self.out_algo) if (str(th) in file) & ('VH' in file)
                        & (period in file) & (method in file)]
        high_vv_list = [file for file in os.listdir(self.out_algo) if (str(th) in file) & ('VV' in file)
                        & (period in file) & (method in file)]
        low_vh_list = [file for file in os.listdir(self.out_algo) if (str(tl) in file) & ('VH' in file)
                       & (period in file) & (method in file)]
        low_vv_list = [file for file in os.listdir(self.out_algo) if (str(tl) in file) & ('VV' in file)
                       & (period in file) & (method in file)]

        for vh, vv in zip(high_vh_list, high_vv_list):

            if method == 'multi':

                remove_high_nb_change(cross_nb(self.out_algo, vh, vv, mf_shp),
                                      nb_change(self.out_algo, "multi_" + period + "_VH_95_.tif", period), nb_max)
            else:
                cross_nb(self.out_algo, vh, vv, mf_shp)

        for vh, vv in zip(low_vh_list, low_vv_list):

            if method == 'multi':

                remove_high_nb_change(cross_nb(self.out_algo, vh, vv, mf_shp),
                                      nb_change(self.out_algo, "multi_" + period + "_VH_95_.tif", period), nb_max)
            else:
                cross_nb(self.out_algo, vh, vv, mf_shp)

        if mf_shp is not None:
            low_file = [ras for ras in os.listdir(os.path.join(self.out_algo, "Post_raster")) if (str(tl) in ras)
                        & (period in ras) & (method in ras) & ("rmv" in ras) & ("crossP" in ras)]
            high_file = [ras for ras in os.listdir(os.path.join(self.out_algo, "Post_raster")) if (str(th) in ras)
                         & (period in ras) & (method in ras) & ("rmv" in ras) & ("crossP" in ras)]
        else:

            low_file = [ras for ras in os.listdir(os.path.join(self.out_algo, "Post_raster")) if (str(tl) in ras)
                        & (period in ras) & (method in ras) & ("crossP" in ras)]
            high_file = [ras for ras in os.listdir(os.path.join(self.out_algo, "Post_raster")) if (str(th) in ras)
                         & (period in ras) & (method in ras) & ("crossP" in ras)]

        for ras_low, ras_high in zip(low_file, high_file):
            polygonize_raster(self.out_algo, ras_low)
            polygonize_raster(self.out_algo, ras_high)

            if mf_shp is not None:

                shp_low = [shp for shp in os.listdir(os.path.join(self.out_algo, "Post_shp")) if (str(tl) in shp)
                           & (period in shp) & (method in shp) & ("rmv" in shp) & ("crossP" in shp) & (
                               shp.endswith(".shp"))][0]

                shp_high = [shp for shp in os.listdir(os.path.join(self.out_algo, "Post_shp")) if (str(th) in shp)
                            & (period in shp) & (method in shp) & ("rmv" in shp) & ("crossP" in shp) & (
                                shp.endswith(".shp"))][0]
            else:
                shp_low = [shp for shp in os.listdir(os.path.join(self.out_algo, "Post_shp")) if (str(tl) in shp)
                           & (period in shp) & (method in shp) & ("crossP" in shp) & (
                               shp.endswith(".shp"))][0]

                shp_high = [shp for shp in os.listdir(os.path.join(self.out_algo, "Post_shp")) if (str(th) in shp)
                            & (period in shp) & (method in shp) & ("crossP" in shp) & (
                                shp.endswith(".shp"))][0]

            repair_shapefile(os.path.join(self.out_algo, "Post_shp", shp_low))
            repair_shapefile(os.path.join(self.out_algo, "Post_shp", shp_high))

            cross_tc(self.out_algo, shp_high, shp_low, tl, th, area_th, area_tl)

            results_algo = [f for f in os.listdir(os.path.join(self.out_algo, "Post_shp", "crossTc")) if
                            (period in f) & ('.shp' in f) & (str(th) + '_' + str(tl) in f)]
            erosion_dilation(self.out_algo, results_algo[0], os.path.join(self.out_algo, "Post_shp", "crossTc"),
                             'EPSG:3857', period, 20, th, tl, method)
