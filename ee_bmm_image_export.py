#--------------------------------
# Name:    ee_beamerminor_image_export.py
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

import ee_bmm.inputs as inputs
import ee_bmm.utils as utils


"""
This tool exports a Beamer-Minor long-term median image to Google Drive.

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
    
    # HA number
    ha = ini['INPUTS']['hyd_area']

    # Google Drive folder
    gdrive = ini['ZONAL_STATS']['output_folder']

    # image storage location (earth engine asset or COGs in cloud storage)
    img_storage_location = ini['IMAGES']['image_storage_location']
    
    # asset ID location or cloud storage location
    output_path = ini['IMAGES']['export_dest']

    # image band list
    band_nm_list = ini['IMAGES']['band_list']

    # water index image collection id
    ndwi_id = ini['IMAGES']['water_index_id']

    # asset ID location or cloud storage location
    output_path_n = f"{output_path}/{ha}"

    # Beamer Minor image collection
    bmm_col = ee.ImageCollection(output_path_n)
    # if img_storage_location == 'ee_asset':
        
        # bmm_col = ee.ImageCollection(output_path_n)
        
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

    #     # list of loaded geotiff images and set a couple of properties we need
    #     bmm_img_list = [ee.Image.loadGeoTIFF(uri).set({'system:index': f"{ha}_{uri.split('/')[-1].split('.')[0][-4:]}", 'year': uri.split('/')[-1].split('.')[0][-4:]}) for uri in cog_uris]

    #     bmm_col = ee.ImageCollection(bmm_img_list)
        
    # list of years from the collection
    yr_lst = bmm_col.aggregate_array('year').getInfo()
    start_yr = yr_lst[0]
    end_yr = yr_lst[-1]
    
    # CRS and Transform for export
    projection = bmm_col.first().projection()
    ee_json = ee.Dictionary(ee.Algorithms.Describe(projection))
    ee_crs = ee_json.get('crs')
    ee_transform = ee.List(ee_json.get('transform'))
    
    # Beamer Minor long term mean image (collection)
    bmm_med = (
        ee.ImageCollection(bmm_col
            .select(band_nm_list)
            .median()
            .set('year', 'long_term_median'))
    )
        
    # read phreatophyte shapefile geometry into EE feature collection
    df = gpd.read_file(phreat).to_crs(epsg=4326)
    geo_json = df.to_json()
    phreat_col = ee.FeatureCollection(json.loads(geo_json))

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
        
        # convert image to ft
        if not (band_nm == 'ndvi' or band_nm == 'et_star'):
            out_img = bmm_med.select([band_nm]).first().divide(304.8)
            filename_o = f'{ha}_{shp_nm}_{band_nm}_long_term_median_ft'
        else:
            out_img = bmm_med.select([band_nm]).first()
            filename_o = f'{ha}_{shp_nm}_{band_nm}_long_term_median'

        # export task
        out_task = ee.batch.Export.image.toDrive(**{
            'image': out_img,
            'description': f'{ha}_{band_nm}_LongTermMedianImageExport',
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
        description='Earth Engine Beamer-Minor Method long-term median image export',
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
