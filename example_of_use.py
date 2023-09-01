from cusum.cusum_cls import CuSum_

USER = r"your_username"  # id account for asf https://search.asf.alaska.edu/#/
PASS = r"your_password"

study_site = "POLYGON((46.2726 -15.9429,46.3306 -15.9429,46.3306 -15.8863,46.2726 -15.8863,46.2726 -15.9429))"

forest_mask_shp = r"path_to_your_shapefile_forest_mask"

start_date = "04-01-2022"  # MMDDYYYY format
last_date = "11-15-2023"

path = r"C:\Users\pfefer\guyane"

if __name__ == '__main__':
    obj = CuSum_(path, start_date, last_date, study_site)

    obj.sentinel1_download(start_date, last_date, USER, PASS, "DESCENDING", relativeOrbit=39, platform="sentinel-1A")
    # relativeOrbit should be found by searching images manually on https://search.asf.alaska.edu/#/

    obj.preprocess()

    obj.run([1, 0.95], method='single', max_samples=500)

    obj.post_processing(method='single', th=1, tl=0.95, mf_shp=forest_mask_shp, nb_max=20, area_th=100, area_tl=300)
