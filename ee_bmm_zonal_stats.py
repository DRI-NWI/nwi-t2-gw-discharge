#--------------------------------
# Name:    ee_beamerminor_zonal_stats.py
# Desc:    Zonal stats of Beamer Minor Method ET*, NDVI, PPT, ETo, ETa, and ETg long-term median images 
#            to Earth Engine asset for a phreatophyte boundary
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
This tool exports a summary table of zonal stats for Beamer-Minor long-term images to Google Drive.

It calculates the spatial average NDVI, ET*, WY PPT, WY ETo, ETa, and ETg long-term image
values for each feature in the phreatophyte shapefile.

"""

pp = pprint.PrettyPrinter(indent=4)


def main(ini_path=None):
    """ Earth Engine Image Export

    Parameters
    ----------
    ini_path : str

    """

    logging.info('\nEarth Engine zonal stats for phreatophyte boundary')

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
    shp_nm = phreat_nm.split(".")[0]

    # landsat NDVI harmonization type
    harmonization_typ = ini['IMAGES']['harmonization_type']
    
    # local output path for getDownloadUrl export
    # local output path/locatin where zonal stats are directly exported to
    local_o_path = ini['INPUTS']['local_output_path']
    out_stats_path = os.path.join(local_o_path, f'stats_{harmonization_typ}')
    if not os.path.isdir(out_stats_path):
        os.mkdir(out_stats_path)

    # attribute to use as the unique identifier, default FID
    zone_nm = ini['INPUTS']['zone_field']
    
    # HA boundary
    ha = ini['INPUTS']['hyd_area']

    # start and end years
    start_yr = ini['INPUTS']['start_year']
    end_yr = ini['INPUTS']['end_year']
    year_filter_list = list(range(start_yr, end_yr+1))
    
    # Google Drive folder
    gdrive = ini['ZONAL_STATS']['output_folder']

    # Place-based ETg rate flag to use replacement ETg rates burned into an ETg image for zonal stats
    replacement_etg_flag = ini['ZONAL_STATS']['replacement_rate_etg_flag']

    # image storage location (earth engine asset or COGs in cloud storage)
    img_storage_location = ini['IMAGES']['image_storage_location']
    
    # asset ID location or cloud storage location
    output_path = ini['IMAGES']['export_dest']
    # output_path = 'projects/dri-nwi-425400/assets/bmm' # hard code image collection (and registered COGs, from cloud storage) paths
    
    # water index image collection id
    ndwi_id = ini['IMAGES']['water_index_id']

    # minimum ETg threshold to use for pixels whose NDVI is too low to predict ET*
    min_etg_thresh = ini['IMAGES']['min_etg_threshold']
    min_etg_thresh_mm = min_etg_thresh * 304.8

    # image band list
    band_nm_list = ini['IMAGES']['band_list']

    # asset ID location or cloud storage location
    output_path_n = f"{output_path}/{ha}_{harmonization_typ}"

    # function to calc area of feature
    def calcArea(ftr):
        return ftr.set({
            'area_ee_acres': ftr.geometry().area().divide(4046.86)
        })

    # need to convert units on replacement rate column
    def converUnit(ftr):
        return ftr.set({
            'rplc_rt': ee.Number(ftr.get('rplc_rt_old')).multiply(304.8)
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

    # if shapefile geometries are too complex, creation of an feature collection from manual upload is required
    # try: 
    #     phreat_col_pre = (
    #         ee.FeatureCollection('users/bminor-dri/NWI/Task2_groundwater_discharge/shapefiles/HA117_Fish_Lake_Valley')
    #         # ee.FeatureCollection('projects/ee-bminor/assets/NWI/HA101_CarsonDesert_phreat')
    #         # ee.FeatureCollection('projects/ee-bminor/assets/NWI/HA173B_Railroad_Valley_North')
    #             .select(['FID2', 'rplc_rt', 'scale_fctr', 'bsnAv_flag'], ['FID2', 'rplc_rt_old', 'scale_fctr', 'bsnAv_flag'])
    #             .map(converUnit)
    #     )
    # except Exception as e:
    #     print(e)

    # read phreatophyte shapefile geometry into EE feature collection and modify as needed (above method for complex geometries if it fails)
    try:
        df = gpd.read_file(phreat).to_crs(epsg=4326)
        df = df.loc[df['HA_number'] == ha]
        if len(df) == 0:
            logging.error(f'      No phreatophyte boundary found for HA {ha}, make sure there is a feature available in the shapefile')

        df[zone_nm] = df.index
    
        # Need to first convert replacement rate field units to millimeters before making the feature collection
        if replacement_etg_flag:
    
            for band_nm in band_nm_list:
                if not (band_nm == 'etg' or band_nm == 'etg_lci' or band_nm == 'etg_uci' or band_nm == 'etg_adj' or band_nm == 'etg_lci_adj' or band_nm == 'etg_uci_adj'):
                    print(f'Replacement rates only apply to ETg images, skipping {band_nm}')
                    continue
                
                # shapefile replacement ETg rate burn field
                burn_field = burn_nm_dict[band_nm]
    
                # convert the units from feet to millimeters
                df[burn_field] = df[burn_field] * 304.8
        
        # convert dataframe to json for earth engine feature collection
        geo_json = df.to_json()
        phreat_col_pre = ee.FeatureCollection(json.loads(geo_json))
    except Exception as e:
        print(e)


    phreat_col = phreat_col_pre.map(calcArea)

    
    # Beamer Minor long term median image
    bmm_med = (
        bmm_col
            .filter(ee.Filter.inList('year', year_filter_list))
            .median()
            .set('year', 'long_term_median')
    )
    
    if replacement_etg_flag:

        # list of all the possible band names in the collection
        all_band_nm_list = ['ndvi', 'et_star', 'eto', 'pr', 'eta', 'eta_lci', 'eta_uci', 'etg', 'etg_lci', 'etg_uci', 'etg_adj', 'etg_lci_adj', 'etg_uci_adj']

        # non place based band list, we can't process these bands so need to keep them away from further processing if they were provided in the INI list
        non_plc_band_nm_list = ['ndvi', 'et_star', 'eto', 'pr', 'eta', 'eta_lci', 'eta_uci']

        # place based band list that are able to be processed
        plc_band_nm_list = ['etg', 'etg_lci', 'etg_uci', 'etg_adj', 'etg_lci_adj', 'etg_uci_adj']
        
        # create a list of band names without the bands we are modifiying with replacement rates
        band_names = [x for x in all_band_nm_list if x in non_plc_band_nm_list]
        mis_band_names = [x for x in all_band_nm_list if x in plc_band_nm_list and x not in band_nm_list]
        band_names = band_names + mis_band_names
        
        bmm_med_sub = (
            bmm_med
                .select(band_names)
        )
        
        for band_nm in band_nm_list:

            if not (band_nm == 'etg' or band_nm == 'etg_lci' or band_nm == 'etg_uci' or band_nm == 'etg_adj' or band_nm == 'etg_lci_adj' or band_nm == 'etg_uci_adj'):
                print(f'Replacement rates only apply to ETg images, skipping {band_nm}')
                continue 
            
            # shapefile replacement ETg rate burn field
            burn_field = burn_nm_dict[band_nm]

            # select band from the current loop
            out_img = bmm_med.select([band_nm])

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

            # use a constant image to burn in values for missing pixels (most often associated saturated pixels in playas/bare soil areas, e.g., Carson Desert)
            fill_img = ee.Image.constant(min_etg_thresh_mm).rename([band_nm]).updateMask(out_burned_img.mask().Not())
            out_burned_img_filled = ee.ImageCollection([out_burned_img, fill_img]).mosaic()
            
            out_burned_img_filled = out_burned_img_filled.where(out_burned_img_filled.select(band_nm).lt(min_etg_thresh_mm), min_etg_thresh_mm)
            
            # add the burned image band back into the stack
            bmm_med_sub = (
                bmm_med_sub
                    .addBands(out_burned_img_filled)
            )
    else:

        bmm_med_sub = bmm_med
    
    def removGeom(ftr):
        return ftr.setGeometry(None)

    logging.info(f"\nExporting zonal stats of Beamer-Minor long-term ETg images {harmonization_typ} to Google Drive")
    logging.info(f"\nReplacement ETg flag used: {replacement_etg_flag}")
    logging.info(f"\nLandsat Harmonization Type: {harmonization_typ}")
    logging.info(f"\nPhreatophyte boundary: {phreat_nm}")
    logging.info(f"\nPhreatophyte zone field: {zone_nm}")
    logging.info(f"\nHA number: {ha}")

    # reduceRegions
    stats_coll = (
        bmm_med_sub
            .reduceRegions(
                collection=phreat_col,
                reducer=(
                    ee.Reducer.mean()
                        .combine(ee.Reducer.count(), '', True)
                ),
                crs=ee_crs,
                crsTransform=ee_transform,
                tileScale=8
            )
            .map(removGeom)
    )

    # output filename with start/end years
    out_filename = f'BMM_etg_{shp_nm}_regions_{harmonization_typ}_{start_yr}_{end_yr}'
    
    # download/export stats directly to local folder
    try:
        downloadUrlStats = stats_coll.getDownloadURL(**{
            'filetype': 'CSV',
            'filename': out_filename,
            'selectors': [zone_nm, 'area_ee_acres', 'ndvi_count', 'ndvi_mean', 'et_star_mean', 'eto_mean', 'pr_mean', 'eta_mean', 'eta_lci_mean', 
                                   'eta_uci_mean', 'etg_mean', 'etg_lci_mean', 'etg_uci_mean', 'etg_adj_mean', 'etg_lci_adj_mean', 'etg_uci_adj_mean'],
            # 'selectors': ['FID2', 'area_ee_acres', 'ndvi_count', 'ndvi_mean', 'et_star_mean', 'eto_mean', 'pr_mean', 'eta_mean', 'eta_lci_mean', 
            #                       'eta_uci_mean', 'etg_mean', 'etg_lci_mean', 'etg_uci_mean', 'etg_adj_mean', 'etg_lci_adj_mean', 'etg_uci_adj_mean'], # these selectors were for CarsonDesert101 run
        })
        responseStats = requests.get(downloadUrlStats)
        if responseStats.status_code == 200:
            with open(os.path.join(out_stats_path, f'{out_filename}.csv'), 'wb') as f:
                f.write(responseStats.content)
            print(f"zonal stats table downloaded successfully to the 'stats_{harmonization_typ}' folder")
        else:
            print(f"Failed to download zonal stats table. Status code: {responseStats.status_code}")
        logging.info(f"\ntable exported for HA: {ha}, phreatophyte shapefile: {phreat_nm}, location: {out_stats_path}")

    except Exception as e:
        print('attempt to use getDownloadURL export failed, submitting batch table export instead')
        # batch table export
        out_table = ee.batch.Export.table.toDrive(**{
          'collection': stats_coll,
          'description': f'{ha}_BMM_TableExportRegions_{harmonization_typ}',
          'folder': gdrive,
          'fileNamePrefix': out_filename,
          'fileFormat': 'CSV',
          'selectors': [zone_nm, 'area_ee_acres', 'ndvi_count', 'ndvi_mean', 'et_star_mean', 'eto_mean', 'pr_mean', 'eta_mean', 'eta_lci_mean', 
                                'eta_uci_mean', 'etg_mean', 'etg_lci_mean', 'etg_uci_mean', 'etg_adj_mean', 'etg_lci_adj_mean', 'etg_uci_adj_mean'],
          # 'selectors': ['FID2', 'area_ee_acres', 'ndvi_count', 'ndvi_mean', 'et_star_mean', 'eto_mean', 'pr_mean', 'eta_mean', 'eta_lci_mean', 
          #                       'eta_uci_mean', 'etg_mean', 'etg_lci_mean', 'etg_uci_mean', 'etg_adj_mean', 'etg_lci_adj_mean', 'etg_uci_adj_mean'],
        })
        out_table.start() 
        logging.info(f"\ntask started for HA: {ha}, phreatophyte shapefile: {phreat_nm}")
        

def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Earth Engine Beamer-Minor Method zonal stats',
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
