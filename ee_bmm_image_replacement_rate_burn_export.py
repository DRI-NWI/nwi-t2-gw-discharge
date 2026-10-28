#--------------------------------
# Name:    ee_beamerminor_image_replacement_rate_burn_export.py
# Desc:    Export Beamer Minor Method long-term median ET*, NDVI, PPT, ETo, ETa, or ETg image
#            to Google Drive for a phreatophyte shapefile
# Python:  3.11
#--------------------------------

import argparse
import logging
import ee
from google.cloud import storage
# import multiprocessing
# import copy
# from retry import retry
import pprint
import sys
import json

import os
import pandas as pd
import geopandas as gpd
import numpy as np
import datetime
import requests

import ee_bmm.inputs as inputs
import ee_bmm.utils as utils


"""
This tool prepares a Beamer-Minor long-term median image and burns 
replacement ETg rates into the image before exporting.

It exports a long-term median NDVI, ET*, WY PPT, WY ETo, ETa, and ETg 
image for a phreatophyte boundary in Nevada and exports them to Google Drive. 
This is the simplified Beamer-Minor Method based on mid-summer median
NDVI composites from Landsat Collection 2 Surface Reflectance data

"""

pp = pprint.PrettyPrinter(indent=4)


def main(ini_path=None):
    """ Earth Engine Image Export

    Parameters
    ----------
    ini_path : str

    """

    logging.info('\nEarth Engine image export for a phreatophyte boundary')

    # Read config file
    ini = inputs.read(ini_path)
    inputs.parse_section(ini, section='INPUTS')
    inputs.parse_section(ini, section='ZONAL_STATS')
    inputs.parse_section(ini, section='IMAGES')

    # Initialize Earth Engine API key
    logging.info('\nInitializing Earth Engine')
    
    # high volume endpoint
    ee.Initialize(project=ini['INPUTS']['gcloud_project_id'])

    # project folder to output data in
    #proj_fol = ini['ZONAL_STATS']['output_ws'].split('\\')[-1]

    # input phreatophyte shapefile
    phreat = ini['INPUTS']['zone_shp_path']
    phreat_nm = phreat.split("\\")[-1]

    shp_nm = phreat_nm.split('.')[0]

    # local output path for getDownloadUrl export
    # local output path/locatin where zonal stats are directly exported to
    local_o_path = ini['INPUTS']['local_output_path']
    out_stats_path = os.path.join(local_o_path, 'stats')
    if not os.path.isdir(out_stats_path):
        os.mkdir(out_stats_path)
    
    # HA number
    ha = ini['INPUTS']['hyd_area']

    # start and end years
    start_yr = ini['INPUTS']['start_year']
    end_yr = ini['INPUTS']['end_year']
    year_filter_list = list(range(start_yr, end_yr+1))
    
    # Google Drive folder
    gdrive = ini['ZONAL_STATS']['output_folder']

    # image storage location (earth engine asset or COGs in cloud storage)
    img_storage_location = ini['IMAGES']['image_storage_location']
    
    # asset ID location or cloud storage location for image export
    output_path = ini['IMAGES']['export_dest']

    # landsat NDVI harmonization type
    harmonization_typ = ini['IMAGES']['harmonization_type']
    
    # image band list
    band_nm_list = ini['IMAGES']['band_list']

    # water index image collection id
    ndwi_id = ini['IMAGES']['water_index_id']

    # minimum ETg threshold to use for pixels whose NDVI is too low to predict ET*
    min_etg_thresh = ini['IMAGES']['min_etg_threshold']
    min_etg_thresh_mm = min_etg_thresh * 304.8

    # asset ID location or cloud storage location
    output_path_n = f"{output_path}/{ha}_{harmonization_typ}"

    # function to calc area of feature
    def calcArea(ftr):
        return ftr.set({
            'area_ee_acres': ftr.geometry().area().divide(4046.86)
        })

    # Beamer Minor image collection
    bmm_col = ee.ImageCollection(output_path_n)
    # if img_storage_location == 'ee_asset':
        
    #     bmm_col = ee.ImageCollection(output_path_n)
        
    # elif img_storage_location == 'cloud_storage':

    #     # output location path split
    #     out_split = output_path_n.split('/')
        
    #     # Instantiate a Cloud Storage client
    #     storage_client = storage.Client()
        
    #     # cloud storage bucket prefix (path)
    #     prefix = f"{out_split[1]}/{out_split[2]}/"
        
    #     # List blobs (objects) in the specified bucket, with optional prefix filtering
    #     blobs = storage_client.list_blobs(out_split[0], prefix=prefix)

    #     # list of uri to all images on the bucket
    #     cog_uris = [f"gs://{out_split[0]}/{blob.name}" for blob in blobs]

    #     # list of loaded geotiff images 
    #     bmm_img_list = [ee.Image.loadGeoTIFF(uri).set({'system:index': f"{ha}_{uri.split('/')[-1].split('.')[0][-4:]}", 'year': uri.split('/')[-1].split('.')[0][-4:]}) for uri in cog_uris]

    #     bmm_col = ee.ImageCollection(bmm_img_list)

    # replacement rate field name dictionary
    burn_nm_dict = {
        'etg_adj': 'rplc_rt',
        'etg_lci_adj': 'rplc_rt_lc',
        'etg_uci_adj': 'rplc_et_uc',
    }
    
    # CRS and Transform for export
    projection = bmm_col.first().projection()
    ee_json = ee.Dictionary(ee.Algorithms.Describe(projection))
    ee_crs = ee_json.get('crs')
    ee_transform = ee.List(ee_json.get('transform'))
    
    # Beamer Minor long term mean image (collection)
    bmm_med = (
        ee.ImageCollection(bmm_col
            .select(band_nm_list)
            .filter(ee.Filter.inList('year', year_filter_list))
            .median()
            .set('year', 'long_term_median'))
    )
    # carson desert boundary too large/complex to read in from dataframe, using preloaded ee asset
    # fish lake valley as well
    # try: 
    #     # phreat_col_pre = ee.FeatureCollection('projects/ee-bminor/assets/NWI/HA101_CarsonDesert_phreat')
    #     phreat_col_pre = ee.FeatureCollection('users/bminor-dri/NWI/Task2_groundwater_discharge/shapefiles/HA117_Fish_Lake_Valley')
    # except Exception as e:
    #     print(e)
    
    # read phreatophyte shapefile geometry into EE feature collection
    df = gpd.read_file(phreat).to_crs(epsg=4326)
    df = df.loc[df['HA_number'] == ha]
    if len(df) == 0:
        logging.error(f'      No phreatophyte boundary found for HA {ha}, make sure there is a feature available in the shapefile')
    geo_json = df.to_json()
    phreat_col_pre = ee.FeatureCollection(json.loads(geo_json))
    
    phreat_col = phreat_col_pre.map(calcArea)

    logging.info(f"\nExporting Beamer-Minor long-term median image(s) ({band_nm_list}) to Google Drive folder {gdrive}")
    logging.info(f"\nPhreatophyte boundary: {phreat_nm}")
    logging.info(f"\nHA number: {ha}")
    logging.info(f"\nStart Year: {start_yr}")
    logging.info(f"\nEnd Year: {end_yr}")
    if img_storage_location == 'ee_asset':
        logging.info(f"\nAsset ID: {output_path_n}")
    elif img_storage_location == 'cloud_storage':
        logging.info(f"\nCOGs location: {output_path_n}")

    for band_nm in band_nm_list:
        
        if not (band_nm == 'etg' or band_nm == 'etg_lci' or band_nm == 'etg_uci' or band_nm == 'etg_adj' or band_nm == 'etg_lci_adj' or band_nm == 'etg_uci_adj' or band_nm == 'etg_lpi_adj' or band_nm == 'etg_upi_adj'):
            print(f'Replacement rates only apply to ETg images, skipping {band_nm}')
            continue

        # shapefile replacement ETg rate burn field
        burn_field = burn_nm_dict[band_nm]
        
        # convert image to ft
        out_img = bmm_med.select([band_nm]).first().divide(304.8)
        filename_o = f'HA{ha}_{shp_nm}_{band_nm}_long_term_median_ft_burned_rplc_rates_{harmonization_typ}_{start_yr}_{end_yr}'

        # filter phreatophyte feature collection to fields/polygons we are using scaling factors for
        scl_ftr_coll = (
            phreat_col
                .filter(ee.Filter.gt('scale_fctr', 0))
        )

        # filter phreatophyte feature collection to fields/polygons we are replacing ETg rates for (place-based)
        plc_ftr_coll = (
            phreat_col
                .filter(ee.Filter.gt(burn_field, 0))
        )

        # filter phreatophyte feature collection to fields/polygons we are calculating the basin-wide average ETg from
        bsnAv_non_ftr_coll = (
            phreat_col
                .filter(ee.Filter.eq('bsnAv_flag', 0))
        )

        # filter phreatophyte feature collection to fields/polygons we are replacing ETg rates for with the basin-wide average ETg
        bsnAv_burn_ftr_coll = (
            phreat_col
                .filter(ee.Filter.eq('bsnAv_flag', 1))
        )

        # check the size of the filtered feature collections to see if place-based rates, scaling factors, and/or basin-wide average flags are provided
        scl_ftr_coll_size = scl_ftr_coll.size().getInfo() 
        plc_ftr_coll_size = plc_ftr_coll.size().getInfo()
        bsnAv_burn_ftr_coll_size = bsnAv_burn_ftr_coll.size().getInfo()

        # check to see if replacement rates and scaling factors have been found in the shapefile and proceed accordingly
        if (scl_ftr_coll_size == 0 and plc_ftr_coll_size == 0 and bsnAv_burn_ftr_coll_size == 0):
            print('no ETg scaling factors, place-based replacement ETg rates, or basin-wide average replacement ETg rate flags found in the phreatophyte boundary. Please check the phreatophyte shapefile')
            continue
        elif (scl_ftr_coll_size == 0 and plc_ftr_coll_size > 0 and bsnAv_burn_ftr_coll_size == 0):
            print('no ETg scaling factors or basin-wide average replacement ET rate flags were found, but place-based replacement ETg rates were found in the phreatophyte boundary. Proceeding with just place-based replacement ETg rates')
    
            # convert the feature collection to an image with the replacement ETg rates as the pixel values
            plc_burn_img = (
                plc_ftr_coll
                    .reduceToImage([burn_field], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )
    
            # mosaic the ETg image and the place-based image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_img, plc_burn_img]).mosaic()  

        elif (scl_ftr_coll_size > 0 and plc_ftr_coll_size == 0 and bsnAv_burn_ftr_coll_size == 0):
            print('ETg scaling factors were found, but no place-based replacement ETg rates or basin-wide average replacement ETg rate flags were found in the phreatophyte boundary. Proceeding with just ETg scaling factors')

            # build a scaling image to use to adjust ETg values
            scl_burn_img = (
                scl_ftr_coll
                    .reduceToImage(['scale_fctr'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # multiply the scaling image by the ETg image
            out_scl_burn_img = out_img.multiply(scl_burn_img)

            # mosaic the ETg image and the scaled image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_img, out_scl_burn_img]).mosaic()

        elif (scl_ftr_coll_size == 0 and plc_ftr_coll_size == 0 and bsnAv_burn_ftr_coll_size > 0):
            
            # calculate the average ETg for the filtered features to get at a basin-wide average
            Av_stats = (
                out_img
                    .reduceRegions(
                        collection=bsnAv_non_ftr_coll,
                        reducer=ee.Reducer.mean(),
                        crs=ee_crs,
                        crsTransform=ee_transform,
                    )
            )
            
            # reduceColumns to calculate the basin-wide average from each feature's average ETg
            bsnAv_stat = (
                Av_stats
                    .reduceColumns(
                        reducer=ee.Reducer.mean(),
                        selectors=['mean'],
                        weightSelectors=['area_ee_acres'],
                    )
            )
              
            # filter phreatophyte feature collection fields/polygons we need to assign the basin-wide average ETg to
            # also set the basin-wide average ETg value as a new property
            def setBsnAv(ftr):
                return (
                    ftr.set({
                        'bsnAv_val': bsnAv_stat.get('mean')
                    })
                )

            # map the function
            bsnAv_burn_ftr_coll = bsnAv_burn_ftr_coll.map(setBsnAv)

            # convert the feature collection to an image with the basin-wide average replacement ETg rate as the pixel values
            bsnAv_burn_img = (
                bsnAv_burn_ftr_coll
                    .reduceToImage(['bsnAv_val'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # mosaic the scaled and place-based image with the basin-wide average ETg image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_img, bsnAv_burn_img]).mosaic()
        
        elif (scl_ftr_coll_size == 0 and plc_ftr_coll_size > 0 and bsnAv_burn_ftr_coll_size > 0):
            
            # convert the feature collection to an image with the replacement ETg rates as the pixel values
            plc_burn_img = (
                plc_ftr_coll
                    .reduceToImage([burn_field], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )
    
            # mosaic the ETg image, the scaled image, and the place-based image to replace the ETg values for irrigated areas
            out_burned_img_pre = ee.ImageCollection([out_img, plc_burn_img]).mosaic()
            
            # calculate the average ETg for the filtered features to get at a basin-wide average
            Av_stats = (
                out_burned_img_pre
                    .reduceRegions(
                        collection=bsnAv_non_ftr_coll,
                        reducer=ee.Reducer.mean(),
                        crs=ee_crs,
                        crsTransform=ee_transform,
                    )
            )
            
            # reduceColumns to calculate the basin-wide average from each feature's average ETg
            bsnAv_stat = (
                Av_stats
                    .reduceColumns(
                        reducer=ee.Reducer.mean(),
                        selectors=['mean'],
                        weightSelectors=['area_ee_acres'],
                    )
            )
              
            # filter phreatophyte feature collection fields/polygons we need to assign the basin-wide average ETg to
            # also set the basin-wide average ETg value as a new property
            def setBsnAv(ftr):
                return (
                    ftr.set({
                        'bsnAv_val': bsnAv_stat.get('mean')
                    })
                )

            # map the function
            bsnAv_burn_ftr_coll = bsnAv_burn_ftr_coll.map(setBsnAv)

            # convert the feature collection to an image with the basin-wide average replacement ETg rate as the pixel values
            bsnAv_burn_img = (
                bsnAv_burn_ftr_coll
                    .reduceToImage(['bsnAv_val'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # mosaic the scaled and place-based image with the basin-wide average ETg image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_burned_img_pre, bsnAv_burn_img]).mosaic()

        elif (scl_ftr_coll_size > 0 and plc_ftr_coll_size > 0 and bsnAv_burn_ftr_coll_size == 0):

            # build a scaling image to use to adjust ETg values
            scl_burn_img = (
                scl_ftr_coll
                    .reduceToImage(['scale_fctr'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # multiply the scaling image by the ETg image
            out_scl_burn_img = out_img.multiply(scl_burn_img)
            
            # convert the feature collection to an image with the replacement ETg rates as the pixel values
            plc_burn_img = (
                plc_ftr_coll
                    .reduceToImage([burn_field], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )
    
            # mosaic the ETg image, the scaled image, and the place-based image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_img, out_scl_burn_img, plc_burn_img]).mosaic()
        
        elif (scl_ftr_coll_size > 0 and plc_ftr_coll_size == 0 and bsnAv_burn_ftr_coll_size > 0):

            # build a scaling image to use to adjust ETg values
            scl_burn_img = (
                scl_ftr_coll
                    .reduceToImage(['scale_fctr'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # multiply the scaling image by the ETg image
            out_scl_burn_img = out_img.multiply(scl_burn_img)
    
            # mosaic the ETg image and the scaled image to replace the ETg values for irrigated areas
            out_burned_img_pre = ee.ImageCollection([out_img, out_scl_burn_img]).mosaic()
            
            # calculate the average ETg for the filtered features to get at a basin-wide average
            Av_stats = (
                out_burned_img_pre
                    .reduceRegions(
                        collection=bsnAv_non_ftr_coll,
                        reducer=ee.Reducer.mean(),
                        crs=ee_crs,
                        crsTransform=ee_transform,
                    )
            )
            
            # reduceColumns to calculate the basin-wide average from each feature's average ETg
            bsnAv_stat = (
                Av_stats
                    .reduceColumns(
                        reducer=ee.Reducer.mean(),
                        selectors=['mean'],
                        weightSelectors=['area_ee_acres'],
                    )
            )
              
            # filter phreatophyte feature collection fields/polygons we need to assign the basin-wide average ETg to
            # also set the basin-wide average ETg value as a new property
            def setBsnAv(ftr):
                return (
                    ftr.set({
                        'bsnAv_val': bsnAv_stat.get('mean')
                    })
                )

            # map the function
            bsnAv_burn_ftr_coll = bsnAv_burn_ftr_coll.map(setBsnAv)

            # convert the feature collection to an image with the basin-wide average replacement ETg rate as the pixel values
            bsnAv_burn_img = (
                bsnAv_burn_ftr_coll
                    .reduceToImage(['bsnAv_val'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # mosaic the scaled and place-based image with the basin-wide average ETg image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_burned_img_pre, bsnAv_burn_img]).mosaic()
        
        elif (scl_ftr_coll_size > 0 and plc_ftr_coll_size > 0 and bsnAv_burn_ftr_coll_size > 0):

            # build a scaling image to use to adjust ETg values
            scl_burn_img = (
                scl_ftr_coll
                    .reduceToImage(['scale_fctr'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # multiply the scaling image by the ETg image
            out_scl_burn_img = out_img.multiply(scl_burn_img)
            
            # convert the feature collection to an image with the replacement ETg rates as the pixel values
            plc_burn_img = (
                plc_ftr_coll
                    .reduceToImage([burn_field], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )
    
            # mosaic the ETg image, the scaled image, and the place-based image to replace the ETg values for irrigated areas
            out_burned_img_pre = ee.ImageCollection([out_img, out_scl_burn_img, plc_burn_img]).mosaic()
            
            # calculate the average ETg for the filtered features to get at a basin-wide average
            Av_stats = (
                out_burned_img_pre
                    .reduceRegions(
                        collection=bsnAv_non_ftr_coll,
                        reducer=ee.Reducer.mean(),
                        crs=ee_crs,
                        crsTransform=ee_transform,
                    )
            )
            
            # reduceColumns to calculate the basin-wide average from each feature's average ETg
            bsnAv_stat = (
                Av_stats
                    .reduceColumns(
                        reducer=ee.Reducer.mean(),
                        selectors=['mean'],
                        weightSelectors=['area_ee_acres'],
                    )
            )
              
            # filter phreatophyte feature collection fields/polygons we need to assign the basin-wide average ETg to
            # also set the basin-wide average ETg value as a new property
            def setBsnAv(ftr):
                return (
                    ftr.set({
                        'bsnAv_val': bsnAv_stat.get('mean')
                    })
                )

            # map the function
            bsnAv_burn_ftr_coll = bsnAv_burn_ftr_coll.map(setBsnAv)

            # convert the feature collection to an image with the basin-wide average replacement ETg rate as the pixel values
            bsnAv_burn_img = (
                bsnAv_burn_ftr_coll
                    .reduceToImage(['bsnAv_val'], ee.Reducer.first())
                    .select(['first'], [band_nm])
            )

            # mosaic the scaled and place-based image with the basin-wide average ETg image to replace the ETg values for irrigated areas
            out_burned_img = ee.ImageCollection([out_burned_img_pre, bsnAv_burn_img]).mosaic()

        else:
            logging.error(f"\nerror finding the ETg scaling factors, place-based replacement ETg values, and/or basin-wide average replacement ETg rate flags in the shapefile, stop")
            continue

        if band_nm == 'etg_adj':
            # use a constant image to burn in values for missing pixels (most often associated saturated pixels in playas/bare soil areas, e.g., Carson Desert)
            fill_img = ee.Image.constant(min_etg_thresh).rename([band_nm]).updateMask(out_burned_img.mask().Not())
            out_burned_img_filled = ee.ImageCollection([out_burned_img, fill_img]).mosaic()
    
            out_burned_img_final = out_burned_img_filled.where(out_burned_img_filled.select(band_nm).lt(min_etg_thresh), min_etg_thresh)

        # download/export image directly to local folder (fails for larger image regions)
        # downloadUrlTiff = out_burned_img_final.getDownloadURL({
        #     'bands': [band_nm],
        #     'region': phreat_col.geometry(),
        #     'crs': ee_crs,
        #     'crsTransform': ee_transform,
        #     'format': 'GEO_TIFF',
        #     'formatOptions':{
        #         'noData': -9999,
        #         'cloudOptimized': True,
        #     },
        #     'maxPixels': 1e13,
        # })
        # responseImg = requests.get(downloadUrlTiff)
        # if responseImg.status_code == 200:
        #     with open(os.path.join(out_stats_path, f'{filename_o}.tif'), 'wb') as f:
        #         f.write(responseImg.content)
        #     print("image downloaded successfully to the 'stats' folder")
        # else:
        #     print(f"Failed to download image. Status code: {responseImg.status_code}")
        # logging.info(f"\nimage exported for HA: {ha}, phreatophyte shapefile: {phreat_nm}, location: {out_stats_path}")
        
        # export task
        out_task = ee.batch.Export.image.toDrive(**{
            'image': out_burned_img_final,
            'description': f'{ha}_{band_nm}_{harmonization_typ}_LongTermMedianBurnedImageExport',
            'folder': gdrive,
            'fileNamePrefix': filename_o,
            'fileFormat': 'GeoTIFF',
            'formatOptions':{
                'noData': -9999,
                'cloudOptimized': True,
            },
            'crs': ee_crs,
            'crsTransform': ee_transform,
            # 'scale': 30,
            'region': phreat_col.geometry(),
            'maxPixels': 1e13
        })
        out_task.start()      
        logging.info(f"\ntask started for HA: {ha}")
    

def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Earth Engine Beamer-Minor Method long-term median replacement ETg rate burned image export',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        '-i', '--ini', type=utils.arg_valid_file,
        help='Input file', metavar='FILE')
    
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    
    args = parser.parse_args()

    if args.ini and os.path.isfile(os.path.abspath(args.ini)):
        
        args.ini = os.path.abspath(args.ini)
    
    else:
        
        args.ini = utils.get_ini_path(os.getcwd())
    
    return args


if __name__ == '__main__':
    
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)
    logging.info('\n{}'.format('#' * 80))
    log_f = '{:<20s} {}'
    logging.info(log_f.format(
        'Start Time:', datetime.datetime.now().isoformat(' ')))
    logging.info(log_f.format('Current Directory:', os.getcwd()))
    logging.info(log_f.format('Script:', os.path.basename(sys.argv[0])))
    
    main(ini_path=args.ini)
