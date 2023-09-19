import asf_search as asf
from asf_search import ASFSearchResults
import shutil
from pyrasta.raster import Raster
from cdtec.base import single_change_detection, multi_change_detection, single_change_detection_tcs
from cdtec.main import cdtec_run
from cusum.utils import cdtec_run_tcs, geocode
from cusum.post_processing_tools import *
from osgeo import gdal
from dateutil.parser import parse
import pytz


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
            os.mkdir(os.path.join(path, "pre_processed_img"))
        except FileExistsError:
            pass
        finally:
            self.out_pre = os.path.join(path, "pre_processed_img")

        if not os.path.exists(os.path.join(self.out_pre, "VV")):
            os.mkdir(os.path.join(self.out_pre, "VV"))
            self.vv_path = os.path.join(self.out_pre, "VV")

        else:
            self.vv_path = os.path.join(self.out_pre, "VV")

        if not os.path.exists(os.path.join(self.out_pre, "VH")):
            os.mkdir(os.path.join(self.out_pre, "VH"))
            self.vh_path = os.path.join(self.out_pre, "VH")

        else:
            self.vh_path = os.path.join(self.out_pre, "VH")

        date1 = "".join(date_start.split('-'))
        date2 = "".join(date_end.split('-'))
        try:
            os.mkdir(os.path.join(path, "Output_CuSum" + "_" + date1 + "_" + date2))
        except FileExistsError:
            pass
        finally:
            self.out_algo = os.path.join(path, "Output_CuSum" + "_" + date1 + "_" + date2)

    def sentinel1_download(self, date_start, date_end, user, pw, flight_direction, platform, relativeOrbit=None,
                           mode="IW", level="GRD_HD", polarizations=None, exclude_period=None):

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
                                exclude_period :
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
                                exclude_period : str
                                    Define a period if you need to exclude some period in your time series.
                                    (rainy season for example)
                                    

                                """

        results = asf.search(start=date_start, end=date_end, processingLevel=level, flightDirection=flight_direction,
                             beamMode=mode, maxResults=200, intersectsWith=self.wkt, relativeOrbit=relativeOrbit,
                             polarization=polarizations, platform=platform)
        if exclude_period is not None:

            date_start_exclude = parse(exclude_period.split('_')[0]).replace(tzinfo=pytz.UTC)
            date_end_exclude = parse(exclude_period.split('_')[1]).replace(tzinfo=pytz.UTC)
            session = asf.ASFSession().auth_with_creds(user, pw)

            filtered_results = ASFSearchResults()
            for product in results:
                scene_start_time = parse(product.properties["startTime"])

                if not date_start_exclude <= scene_start_time <= date_end_exclude:
                    filtered_results.append(product)

            filtered_results.download(path=self.out_dl, session=session, processes=8)
            session.close()

        else:

            session = asf.ASFSession().auth_with_creds(user, pw)
            results.download(path=self.out_dl, session=session, processes=8)

    def preprocess(self, spacing=10, scaling='dB', t_srs=4326, allow_RES_OSV=True, polarizations='all',
                   removeS1BorderNoiseMethod='ESA', removeS1ThermalNoise=True,
                   demResamplingMethod='BILINEAR_INTERPOLATION',
                   imgResamplingMethod='BILINEAR_INTERPOLATION', DEM='SRTM 1Sec HGT', speckleFilter='Lee Sigma',
                   merge=False, interp=False, remove_temp_files=False):
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
                        merge : bool
                            Set up True if you need to merge same date and same polarization tif. It appears when your
                            study site cover two S1 tiles.
                        interp : bool
                            If you had to merge S1 tiles, it could appear some no data lines between the two S1 tiles
                            merged result. And it could be necessary to fill those no data lines by interpolation.
                            Set up True to
                            process this step. Default is False.
                        remove_temp_files : bool, optional
                            Set to True to conserve disk space by removing all temporary pre-processed files after the
                             preprocessing is complete.
                            Default is False, which retains the temporary files.

                        """
        zip_files = [filename for filename in os.listdir(self.out_dl) if filename.endswith(".zip")]

        for filename in tqdm(zip_files, desc=f"Pre-processing zip files"):
            geocode(os.path.join(self.out_dl, filename), self.out_pre, t_srs, spacing=spacing,
                    polarizations=polarizations, shapefile=self.zone_files,
                    scaling=scaling, removeS1BorderNoise=True, allow_RES_OSV=allow_RES_OSV,
                    removeS1BorderNoiseMethod=removeS1BorderNoiseMethod,
                    removeS1ThermalNoise=removeS1ThermalNoise,
                    speckleFilter=speckleFilter, terrainFlattening=True, export_extra=None, groupsize=8,
                    cleanup=True,
                    gpt_args=None, returnWF=False, demResamplingMethod=demResamplingMethod, demName=DEM,
                    imgResamplingMethod=imgResamplingMethod, alignToStandardGrid=False, refarea='gamma0',
                    clean_edges=False, clean_edges_npixels=1)

        def crop_rasters_to_max_bounds(path):

            max_x_min = float('inf')
            max_x_max = -float('inf')
            max_y_min = float('inf')
            max_y_max = -float('inf')
            for file_name in os.listdir(path):
                if file_name.endswith(".tif"):
                    file_name_path = os.path.join(path, file_name)
                    ds = gdal.Open(file_name_path)

                    x_min = ds.GetGeoTransform()[0]
                    x_max = ds.GetGeoTransform()[0] + ds.GetGeoTransform()[1] * ds.RasterXSize
                    y_min = ds.GetGeoTransform()[3] + ds.GetGeoTransform()[5] * ds.RasterYSize
                    y_max = ds.GetGeoTransform()[3]

                    max_x_min = min(max_x_min, x_min)
                    max_x_max = max(max_x_max, x_max)
                    max_y_min = min(max_y_min, y_min)
                    max_y_max = max(max_y_max, y_max)

                    ds = None  # Close dataset

            for file_name in tqdm(os.listdir(path), desc="Cropping rasters", position=0):
                if file_name.endswith(".tif"):
                    file_crop_path = os.path.join(path, file_name)
                    temp_output_path = os.path.join(path, "temp_cropped_" + file_name)  # Temporary output path

                    ds = gdal.Open(file_crop_path)

                    gdal.Translate(temp_output_path, ds, projWin=[max_x_min, max_y_max, max_x_max, max_y_min])

                    ds = None  # Close dataset

                    os.remove(file_crop_path)
                    os.rename(temp_output_path, file_crop_path)
        if merge:
            file_groups = {}
            list_file_tif = [filename for filename in os.listdir(self.out_pre) if ('.tif' in filename) &
                             ('merged' not in filename)]
            for filename in list_file_tif:
                date_polarization = (re.search(r"\d{8}", filename)[0], re.search("VV|VH", filename)[0])
                if date_polarization in file_groups:
                    file_groups[date_polarization].append(filename)
                else:
                    file_groups[date_polarization] = [filename]

            for date_polarization, filenames in tqdm(file_groups.items(),
                                                     desc=f"Merging {date_polarization} tif files"):
                for i in range(len(filenames)):
                    filename1 = filenames[i]
                    for j in range(i + 1, len(filenames)):
                        filename2 = filenames[j]

                        vh_merged_path = os.path.join(self.vh_path, filename1[:-4] + "_merged.tif")
                        vv_merged_path = os.path.join(self.vv_path, filename1[:-4] + "_merged.tif")
                        out_pre_merged_path = os.path.join(self.out_pre, filename1[:-4] + "_merged.tif")

                        if not (os.path.exists(vh_merged_path) or os.path.exists(vv_merged_path) or os.path.exists(
                                out_pre_merged_path)):
                            ras1 = Raster(os.path.join(self.out_pre, filename1))
                            ras2 = Raster(os.path.join(self.out_pre, filename2))
                            ras_merged = Raster.merge([ras1, ras2], input_no_data=[0], output_no_data=0)
                            ras_merged_path = os.path.join(self.out_pre, filename1[:-4] + "_merged.tif")
                            ras_merged.to_file(ras_merged_path)

            for filename in os.listdir(self.out_pre):
                if ('VH' in filename) & ('.tif' in filename) & ('merged' in filename):
                    shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.vh_path, filename))
                if ('VV' in filename) & ('.tif' in filename) & ('merged' in filename):
                    shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.vv_path, filename))

            if interp:
                vv_files = os.listdir(self.vv_path)
                vh_files = os.listdir(self.vh_path)

                if len(vv_files) != len(vh_files):
                    raise ValueError("The number of files in vv and vh directories does not match.")

                p_interp_bar = tqdm(total=len(vv_files), desc="Filling no data lines in merged files")

                for file_vv, file_vh in zip(vv_files, vh_files):
                    if "merged" in file_vv:
                        for img_path in [os.path.join(self.vv_path, file_vv), os.path.join(self.vh_path, file_vh)]:
                            img = gdal.Open(img_path, gdal.GA_Update)
                            band = img.GetRasterBand(1)
                            gdal.FillNodata(band, None, 3, 0)
                            img = None  # Close dataset

                        p_interp_bar.set_description(f"Filling no data lines in : {file_vv}, {file_vh}")
                        p_interp_bar.update(1)

            crop_rasters_to_max_bounds(self.vh_path)
            crop_rasters_to_max_bounds(self.vv_path)

        else:
            for filename in os.listdir(self.out_pre):
                if ('VH' in filename) & ('.tif' in filename):
                    shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.vh_path, filename))
                if ('VV' in filename) & ('.tif' in filename):
                    shutil.move(os.path.join(self.out_pre, filename), os.path.join(self.vv_path, filename))
            crop_rasters_to_max_bounds(self.vv_path)
            crop_rasters_to_max_bounds(self.vh_path)

        if remove_temp_files:
            # Remove temporary files that are not in self.vv_path and self.vh_path
            for filename in os.listdir(self.out_pre):
                if filename not in os.listdir(self.vv_path) and filename not in os.listdir(self.vh_path):
                    file_path = os.path.join(self.out_pre, filename)
                    try:
                        os.remove(file_path)
                    except PermissionError:
                        pass

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

        FHANDLE = dict(single=single_change_detection,
                       multi=multi_change_detection,
                       single_tcs=single_change_detection_tcs)

        for in_path in [self.vh_path, self.vv_path]:
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

    def post_processing(self, th, tl, method, mf_shp=None, area_th=300, area_tl=1000, nb_max=10, cpu_cores=None,
                        remove_temp_files=None):

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
                        cpu_cores : int
                            Number of cores to multiprocess post_processing cross tc step
                        remove_temp_files : str, optional
                            Choose which output directory to remove: "chunks", "Post_raster", "Post_shp", or "all".
                            Default is None.



                        """
        th = int(th * 100)
        tl = int(tl * 100)

        filename = [file for file in os.listdir(self.out_algo) if (str(th) in file) & ('VH' in file)
                    & (method in file)][0]
        period = re.search(r"\d{8}_\d{8}", filename)[0]

        make_paths(self.out_algo)

        def generate_file_list(folder, threshold, period, method, polarization):
            file_list = [
                file for file in os.listdir(folder) if
                (threshold in file) and
                (polarization in file) and
                (period in file) and
                (method in file)
            ]

            if method == "single":
                file_list = [file for file in file_list if not ('tcs' in file)]

            return file_list

        high_vh_list = generate_file_list(self.out_algo, str(th), period, method, 'VH')
        high_vv_list = generate_file_list(self.out_algo, str(th), period, method, 'VV')
        low_vh_list = generate_file_list(self.out_algo, str(tl), period, method, 'VH')
        low_vv_list = generate_file_list(self.out_algo, str(tl), period, method, 'VV')

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
                        & (period in ras) & (method in ras) & ("crossP" in ras) & ("rmv" not in ras)]
            high_file = [ras for ras in os.listdir(os.path.join(self.out_algo, "Post_raster")) if (str(th) in ras)
                         & (period in ras) & (method in ras) & ("crossP" in ras) & ("rmv" not in ras)]

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

            try:
                os.mkdir(os.path.join(os.path.join(self.out_algo, "Post_shp", "crossTc")))
            except FileExistsError:
                pass
            finally:
                out_cross_tc = os.path.join(os.path.join(self.out_algo, "Post_shp", "crossTc"))

            if cpu_cores is not None:
                try:
                    os.mkdir(os.path.join(self.out_algo, "chunks"))

                except FileExistsError:
                    pass
                finally:
                    out_chunks = os.path.join(self.out_algo, "chunks")

                try:
                    os.mkdir(os.path.join(self.out_algo, "chunks", "crossTc"))
                except FileExistsError:
                    pass

                grid_shapefile = create_custom_grid(os.path.join(self.out_algo, "Post_shp", shp_low),
                                                    out_chunks, cpu_cores)

                spatially_split_shapefile(os.path.join(self.out_algo, "Post_shp", shp_low), grid_shapefile, out_chunks,
                                          str(tl))
                spatially_split_shapefile(os.path.join(self.out_algo, "Post_shp", shp_high), grid_shapefile, out_chunks,
                                          str(th))

                if mf_shp:

                    chunks_high_list = [chunks_file for chunks_file in os.listdir(out_chunks) if ('.shp' in chunks_file)
                                        & ('_' + str(th) + '_' in chunks_file) & (method in chunks_file)
                                        & ("rmv" in chunks_file)]

                    chunks_low_list = [chunks_file for chunks_file in os.listdir(out_chunks) if ('.shp' in chunks_file)
                                       & ('_' + str(tl) + '_' in chunks_file) & (method in chunks_file)
                                       & ("rmv" in chunks_file)]
                else:
                    chunks_high_list = [chunks_file for chunks_file in os.listdir(out_chunks) if ('.shp' in chunks_file)
                                        & ('_' + str(th) + '_' in chunks_file) & (method in chunks_file)
                                        & ("rmv" not in chunks_file)]

                    chunks_low_list = [chunks_file for chunks_file in os.listdir(out_chunks) if ('.shp' in chunks_file)
                                       & ('_' + str(tl) + '_' in chunks_file) & (method in chunks_file)
                                       & ("rmv" not in chunks_file)]

                def extract_chunk_number(filename):
                    match = re.search(r'_chunk_(\d+)_', filename)
                    if match:
                        return int(match.group(1))
                    return -1  # Return -1 if no match is found

                chunks_low_list = sorted(chunks_low_list, key=extract_chunk_number)
                chunks_high_list = sorted(chunks_high_list, key=extract_chunk_number)
                print(chunks_high_list)

                run_parallel_cross_tc(out_chunks, chunks_high_list, chunks_low_list, tl, th, area_th, area_tl,
                                      cpu_cores)

                chunks_list_cross_tc = [chunks for chunks in os.listdir(os.path.join(out_chunks, 'crossTc'))
                                        if chunks.endswith(".shp") & (method in chunks)]

                merge_shapefiles(os.path.join(out_chunks, 'crossTc'),
                                 chunks_list_cross_tc, out_cross_tc)

            else:

                cross_tc(os.path.join(self.out_algo, "Post_shp"), shp_high, shp_low, tl, th, area_th, area_tl)

            results_algo = [f for f in os.listdir(os.path.join(self.out_algo, "Post_shp", "crossTc")) if
                            (period in f) & ('.shp' in f) & (str(th) + '_' + str(tl) in f)]
            erosion_dilation(self.out_algo, results_algo[0], os.path.join(self.out_algo, "Post_shp", "crossTc"),
                             'EPSG:3857', period, 20, th, tl, method)

            def remove_temp_cross_tc_files(directory_path, extensions):
                try:
                    for file_temp in os.listdir(directory_path):
                        file_path = os.path.join(directory_path, file_temp)
                        if any(file_temp.endswith(ext) for ext in extensions) & ('temp' in file_temp):
                            os.remove(file_path)
                except Exception as e:
                    print(f"An error occurred: {str(e)}")

            directory_to_clean = os.getcwd()
            file_extensions_to_remove = ['.shp', '.cpg', '.dbf', '.prj', '.shx']
            remove_temp_cross_tc_files(directory_to_clean, file_extensions_to_remove)

            if remove_temp_files:
                if remove_temp_files == "all":
                    output_dirs = ["chunks", "Post_raster", "Post_shp"]
                elif remove_temp_files in ["chunks", "Post_raster", "Post_shp"]:
                    output_dirs = [remove_temp_files]
                else:
                    raise ValueError(
                        "Invalid value for 'remove_output'. Choose from 'chunks', 'Post_raster', 'Post_shp', or 'all'.")

                for output_dir in output_dirs:
                    output_path = os.path.join(self.out_algo, output_dir)
                    if os.path.exists(output_path):
                        shutil.rmtree(output_path)
