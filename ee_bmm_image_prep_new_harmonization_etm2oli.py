#--------------------------------
# Name:    ee_beamerminor_image_prep_new_harmonization_etm2oli.py
# Desc:    Export/prepare Beamer Minor Method ET*, NDVI, PPT, ETo, ETa, and ETg images 
#            to Earth Engine asset for an HA boundary in Nevada. New Landsat NDVI harmonization ETM 2 OLI polynomial Great Basin corrections
# Python:  3.11
#--------------------------------

import argparse
import logging
import ee
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
This tool exports Beamer-Minor images to Earth Engine assets.

It prepares NDVI, ET*, WY PPT, WY ETo, ETa, and ETg images for a given
HA boundary in Nevada and exports them to Earth Engine assets or cloud storage. This is
the simplified Beamer-Minor Method based on mid-summer NDVI composites
from Landsat Collection 2 Surface Reflectance data

"""

pp = pprint.PrettyPrinter(indent=4)


def main(ini_path=None):
    """ Earth Engine Image Export

    Parameters
    ----------
    ini_path : str

    """

    logging.info('\nEarth Engine image collection preparation for HA boundaries')

    # Read config file
    ini = inputs.read(ini_path)
    inputs.parse_section(ini, section='INPUTS')
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
    
    # HA boundary
    ha = ini['INPUTS']['hyd_area']

    # start and end years
    start_yr = ini['INPUTS']['start_year']
    end_yr = ini['INPUTS']['end_year']
    start_doy = ini['INPUTS']['start_doy']
    end_doy = ini['INPUTS']['end_doy']

    # Landsat path/rows to use
    path_lst = ini['INPUTS']['path_keep_list']
    row_lst = ini['INPUTS']['row_keep_list']

    # gridMET cell buffer flag
    g_buf = ini['INPUTS']['grid_buffer_flag']

    # gridMET cells to skip
    g_skip_lst = ini['INPUTS']['grid_skip_list']

    # image storage location (earth engine asset or COGs in cloud storage)
    img_storage_location = ini['IMAGES']['image_storage_location']

    # BMM regression version for image properties
    bmm_r_vers = ini['IMAGES']['bmm_regress_vers']

    # regression stats/params for calculating confidence intervals (based on the BMM version)
    # B0, B1, A=Tval*RSE, mean of x, SSE, inverse of n, 1+inverse of n
    param_dict = {
        # 'v2p1p3': [-0.035, 1.1385, 0.12253174979396549, 0.21018529657667032, 0.31721084424075263, 0.01639344262295082, 1.0163934426229508] # confidence/prediction interval fix 
        # 'v2p2p1': [], new sites and expansion outside of Great Basin
        'v2p3p1': [-0.0685, 1.1343, 0.11037154434689736, 0.24552620375922582, 0.275380265224164, 0.015384615384615385, 1.0153846153846153] # new landsat harmonization GB ETM to OLI from Alex & more sites: US-LS1 and HIKO
        # 'v2p3p1': [-0.035, 1.1462, 0.1127634360921511, 0.21202746820642543, 0.2874452771158929, 0.015384615384615385, 1.0153846153846153] # new landsat harmonization GB OLI to ETM from Alex & more sites: US-LS1 and HIKO 
    }
    
    # define params/stats based on dictionary above
    B0 = param_dict[bmm_r_vers][0]
    B1 = param_dict[bmm_r_vers][1]
    bmm_eqn = f'ET*={B0} + {B1}NDVI'
    A = param_dict[bmm_r_vers][2]
    u_x = param_dict[bmm_r_vers][3]
    Sxx = param_dict[bmm_r_vers][4]
    inv_n = param_dict[bmm_r_vers][5] # used for confidence intervals
    one_plus_inv_n = param_dict[bmm_r_vers][6] # used for prediction intervals
    logging.info(f'BMM model version: {bmm_r_vers}, ET* = {B0} + {B1}NDVI')
    
    # asset ID location or cloud storage location
    output_path = ini['IMAGES']['export_dest']

    # water index image collection id
    ndwi_id = ini['IMAGES']['water_index_id']

    # minimum ETg threshold to use for pixels whose NDVI is too low to predict ET*
    min_etg_thresh = ini['IMAGES']['min_etg_threshold']
    min_etg_thresh_mm = min_etg_thresh * 304.8

    # asset ID location or cloud storage location
    output_path_n = f"{output_path}/{ha}_etm2oli"

    # earth engine list of years 
    yrs = ee.List.sequence(start_yr, end_yr, 1)

    # python list of years
    yr_lst = list(range(start_yr, end_yr+1))
    
    # earth engine year/doy filters
    YEAR_RANGE = ee.Filter.calendarRange(start_yr, end_yr, 'year')
    DOY_RANGE = ee.Filter.calendarRange(start_doy, end_doy, 'day_of_year')
    
    # flags for masking
    cirrus_flag = True
    dilate_flag = True
    shadow_flag = True
    snow_flag   = True
    
    # Applies scale factors to Landsat Collection 2 Tier 1 image collection
    def applyScaleFactors(img):
        opticalBands = img.select('SR_B.').multiply(0.0000275).add(-0.2)
        thermalBands = img.select('ST_B.*').multiply(0.00341802).add(149.0)
        return (
            img
                .addBands(opticalBands, None, True)
                .addBands(thermalBands, None, True)
        )
    
    # Cloud Mask SR using bits
    def cloudMask(img):
        qa_img = img.select(['QA_PIXEL'])
        cloud_mask = qa_img.rightShift(3).bitwiseAnd(1).neq(0)
        #     .And(qa_img.rightShift(6).bitwiseAnd(3).gte(cloud_confidence)) issue with bit 6
        if cirrus_flag:
            cloud_mask = cloud_mask.Or(qa_img.rightShift(2).bitwiseAnd(1).neq(0))
        if dilate_flag:
            cloud_mask = cloud_mask.Or(qa_img.rightShift(1).bitwiseAnd(1).neq(0))
        if shadow_flag:
            cloud_mask = cloud_mask.Or(qa_img.rightShift(4).bitwiseAnd(1).neq(0))
        if snow_flag:
            cloud_mask = cloud_mask.Or(qa_img.rightShift(5).bitwiseAnd(1).neq(0))
    
        # Flip to set cloudy pixels to 0 and clear to 1
        return img.updateMask(cloud_mask.Not())
    
    def dateProperties(img):
        return (
            img
                .set('year', ee.Number.parse(ee.String(img.get('DATE_ACQUIRED')).slice(0, 4)).toInt())
                .copyProperties(img)
        )
    
    def srNDVI(img):
        return (
            ee.Image(img
                .select(['red','nir'])
                .expression('((b(1) - b(0)) / (b(1) + b(0)))')
                .rename(['ndvi'])
                .set('year', ee.Number.parse(ee.String(img.get('DATE_ACQUIRED')).slice(0, 4)).toInt())
                .copyProperties(img))
        )

    # Alex updated NDVI adjustments for cross-sensor calibration (June, 2026)
    def ndvi_etm_2_oli_func(img):
        return (
            ee.Image(img
                .select(['ndvi'])
                # .expression('0.0055 + (1.2289 * b(0)) - (0.3197 * (b(0)**2)) + (0.0575 * (b(0)**3))') # 3rd order polynomial CONUS NDVI adjustments
                .expression('0.0164 + (1.1656 * b(0)) - (0.2799 * (b(0)**2)) + (0.0896 * (b(0)**3))') # 3rd order polynomial Great Basin NDVI adjustments
                .rename(['ndvi'])
                .set('year',  ee.Number.parse(ee.String(img.get('DATE_ACQUIRED')).slice(0, 4)).toInt())
                .copyProperties(img))
        )
    
    # build collection of median NDVI composites
    def buildNDVIComp(yr):
        return (
            ndvi_col
                .filter(ee.Filter.eq('year', yr))
                .median()
                .set({
                    'start_doy': start_doy,
                    'end_doy': end_doy,
                    'year': yr,
                    'date_ingested': datetime.datetime.today().strftime('%Y-%m-%d'),
                })
        )

    # version history
    # ET Star images from NDVI ET* = -0.035 + 1.1097(NDVI), LCI, UCI, LPI, and UPI 2024-05-24
    # ET Star images from NDVI ET* = -0.035 + 1.1316(NDVI), LCI, UCI, LPI, and UPI 2024-10-22 (BMM v2p1p0)
    # ET Star images from NDVI ET* = -0.035 + 1.1343(NDVI), LCI, UCI, LPI, and UPI 2025-01-02 (BMM v2p1p1) added EGBRX and PKL sites from DeMeo et al. (2003), added CC mesquite site, LONS and ERG-2 spring sites
    # ET Star images from NDVI ET* = -0.035 + 1.1385(NDVI), LCI, UCI, LPI, and UPI 2025-02-06 (BMM v2p1p2) implemented QC flags for EC/BR flux data following AmeriFlux guidelines and Garcia et al. (2015)
    # ET Star images from NDVI ET* = -0.035 + 1.1385(NDVI), LCI, UCI, LPI, and UPI 2026-03-05 (BMM v2p1p3) update/fix for confidence and prediction interval estimation
    # ET Star images from NDVI ET* = -0.0685 + 1.1343(NDVI), LCI, UCI, LPI, and UPI 2026-07-21 (BMM v2p3p1) new landsat harmonization updates GB ETM to OLI from Alex and additional sites HIKO, US-LS1
    # ET Star images from NDVI ET* = -0.035 + 1.1462(NDVI), LCI, UCI, LPI, and UPI 2026-07-21 (BMM v2p3p1) new landsat harmonization updates GB OLI to ETM from Alex and additional sites HIKO, US-LS1
    def buildETStar(img):

        # apply predictive equation to ndvi image
        yhat = (
            img
                .select('ndvi')
                .multiply(ee.Number(B1))
                .add(ee.Number(B0))
        )

        # confidence interval term
        term_ci = (
            img
                .select('ndvi')
                .subtract(ee.Number(u_x))
                .pow(2)
                .divide(ee.Number(Sxx))
                .add(ee.Number(inv_n))
                .sqrt()
        )

        # prediction interval term
        # term_pi = (
        #     img
        #         .select('ndvi')
        #         .subtract(ee.Number(u_x))
        #         .pow(2)
        #         .divide(ee.Number(Sxx))
        #         .add(one_plus_inv_n)
        #         .sqrt()
        
        return (
            img
                .addBands(yhat)
                .addBands(yhat.subtract(term_ci.multiply(ee.Number(A)))) # lower 90% confidence interval
                .addBands(yhat.add(term_ci.multiply(ee.Number(A)))) # upper 90% confidence interval 
                .rename(['ndvi', 'et_star', 'et_star_lci', 'et_star_uci'])
        )

    # build image collection of water year total ETo and PPT
    def buildGridMetWY(yr):
    
        # start of water year
        start_date = ee.Date.fromYMD(ee.Number(yr).subtract(1), 10, 1)
    
        # end of water year, exclusive
        end_date = ee.Date.fromYMD(yr, 10, 1)
    
        # calculate the water year sum
        grid_wy_img = (
            gridmet_col
                .filterDate(start_date, end_date)
                .sum()
        )
    
        return (
            ee.Image(grid_wy_img)
                .set('year', yr)
        )
    
    # calculate the average water year ETo and PPT for the selected/filtered gridMET cells
    def calcAvgGridMet(img):
    
        year = img.get('year')
        
        # spatial average gridMET WY ETo and PPT and counts
        grid_stat = (
            img
                .reduceRegion(
                    reducer=(
                        ee.Reducer.mean().unweighted()
                            .combine(ee.Reducer.count().unweighted(), '', True)
                    ),
                    geometry=grids_final.geometry(),
                    scale=4638.312116386398
                )
        )
    
        grid_stat = grid_stat.combine(img.toDictionary())
    
        return ee.Feature(None, grid_stat)
        
    # set centroid property on gridMET cell feature collection for adjusting the grid selection
    def setCent(ftr):
        return ftr.set({
                'centroid': ftr.geometry().centroid()
            })
    
    
    # read filtered phreatophyte shapefile geometry into EE feature collection
    df = gpd.read_file(phreat).to_crs(epsg=4326)
    df = df.loc[df['HA_number'] == ha]
    if len(df) == 0:
        logging.error(f'      No phreatophyte boundary found for HA {ha}, make sure there is a feature available in the shapefile')

    # dataframe to geojson
    geo_json = df.to_json()

    # feature collection filtered to the basin of interest
    phreat_col = ee.FeatureCollection(json.loads(geo_json))
    # phreat_col = ee.FeatureCollection('projects/ee-bminor/assets/NWI/HA101_CarsonDesert_phreat')

    # gridMET grid feature collection and set the centroid property
    grids = ee.FeatureCollection('users/bminor/gridmet_grids').map(setCent)
    
    # filter gridMET cell feature collection to intersection grids with phreatophyte boundary
    grids_sub = (
        grids
            .filter(
                ee.Filter.intersects(
                    leftValue=phreat_col.geometry(), 
                    rightField=".geo"
                )
            )
    )
    
    # if buffer flag is true, buffer gridMET cells inwards by 1/2 cell (won't work with really small or linear features)
    if g_buf:
        # buffer gridmet cell extent/selection inwards by 1/2 cell
        grids_perim = (
            ee.Feature(grids_sub
                .union()
                .first())
                .buffer(-2000)
        )
    
    else:
        # dont buffer gridmet cell extent/selection inwards, keep at phreatophyte intersection
        grids_perim = (
            ee.Feature(grids_sub
                .union()
                .first())
        )
    
    # final gridMET cells to mask to
    if g_skip_lst:
        grids_final = (
            grids
                .filter(
                    ee.Filter.intersects(
                        leftValue=grids_perim.geometry(), 
                        rightField="centroid"
                    )
                )
                .filter(ee.Filter.inList('GRIDMET_ID', g_skip_lst).Not())
        )

    else:
        grids_final = (
            grids
                .filter(
                    ee.Filter.intersects(
                        leftValue=grids_perim.geometry(), 
                        rightField="centroid"
                    )
                )
                # .filter(ee.Filter.inList('GRIDMET_ID', [438180, 438181, 438182, 439566, 439568, 440953, 440954, 443724, 445109, 445110, 446495, 446496]))
        )

    # UNUSED list and mask
    # list of gridMET cells to include in image properties
    # grid_lst = grids_final.aggregate_array('GRIDMET_ID').getInfo()

    # final gridMET mask to mask ETo and PPT bands to
    # gridmet_mask = ee.FeatureCollection(grids_final).reduceToImage(['GRIDMET_ID'], ee.Reducer.first())
    
    # prepare gridMET ETo and Precip data
    gridmet_col = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET").select(['eto', 'pr'])
        
    # check if image collection exists for the HA, overwrite if it does exist
    if img_storage_location == 'ee_asset':
        try:
            ee.data.createAsset({'type': 'ImageCollection'}, output_path_n)
        except Exception as e:
            print(f'{output_path_n} already exists, attempting to upload images from {start_yr} to {end_yr}')

    # HA boundary to use for image export
    # ha_ftr = (
    #     ee.Feature(ee.FeatureCollection('users/bminor/NWI/shapefiles/HA_admin_boundaries')
    #         .filter(ee.Filter.eq('HYD_AREA', ha))
    #         .select(['HYD_AREA', 'HYD_AREA_N'])
    #         .first())
    # )
    ha_ftr = (
        ee.Feature(ee.FeatureCollection('projects/dri-nwi-425400/assets/shapefiles/ha_boundaries_cleaned_20260819')
            .filter(ee.Filter.eq('BasinID', ha))
            .select(['BasinID', 'BasinName'], ['HYD_AREA', 'HYD_AREA_N'])
            .first())
    )

    # NDWI or MNDWI water pixels (long-term median) in the final image
    ndwi_img = ee.Image(ndwi_id)
    ndwi_mask = ndwi_img.mask().Not()
    
    # prep Landsat image collections
    landsat5 = (
        ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
            .filter(ee.Filter.inList('WRS_PATH', path_lst))
            .filter(ee.Filter.inList('WRS_ROW', row_lst))
            .filter(YEAR_RANGE)
            .filter(DOY_RANGE)
            .map(applyScaleFactors)
            .select(['SR_B3', 'SR_B4', 'QA_PIXEL'], ['red', 'nir', 'QA_PIXEL'])
            .filter(ee.Filter.lt('CLOUD_COVER', 70))
    )

    if start_yr < 2021:
        landsat7 = (
            ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
                .filter(ee.Filter.inList('WRS_PATH', path_lst))
                .filter(ee.Filter.inList('WRS_ROW', row_lst))
                .filter(ee.Filter.calendarRange(start_yr, 2020, 'year')) # Landsat 7 removal starting in 2021 for platform/sensor drift?
                .filter(DOY_RANGE)
                .map(applyScaleFactors)
                .select(['SR_B3', 'SR_B4', 'QA_PIXEL'], ['red', 'nir', 'QA_PIXEL'])
                .filter(ee.Filter.lt('CLOUD_COVER', 70))
        )
        
    landsat8 = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
            .filter(ee.Filter.inList('WRS_PATH', path_lst))
            .filter(ee.Filter.inList('WRS_ROW', row_lst))
            .filter(YEAR_RANGE)
            .filter(DOY_RANGE)
            .map(applyScaleFactors)
            .select(['SR_B4', 'SR_B5', 'QA_PIXEL'], ['red', 'nir', 'QA_PIXEL'])
            .filter(ee.Filter.lt('CLOUD_COVER', 70))
    )

    landsat9 = (
        ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            .filter(ee.Filter.inList('WRS_PATH', path_lst))
            .filter(ee.Filter.inList('WRS_ROW', row_lst))
            .filter(YEAR_RANGE)
            .filter(DOY_RANGE)
            .map(applyScaleFactors)
            .select(['SR_B4', 'SR_B5', 'QA_PIXEL'], ['red', 'nir', 'QA_PIXEL'])
            .filter(ee.Filter.lt('CLOUD_COVER', 70))
    )
    
    
    # merge L5, L7, L8, and L9 collections
    landsat_col = ee.ImageCollection([])
    # not using L7 after 2020 b/c of sensor drift
    if start_yr < 2021:
        landsat_col = ee.ImageCollection(landsat_col.merge(landsat5.merge(landsat7.merge(landsat8.merge(landsat9)))))
    else:
        landsat_col = ee.ImageCollection(landsat_col.merge(landsat5.merge(landsat8.merge(landsat9))))
        
    # cloud masking, ndvi calc, set date properties mapping over the merged collection
    ndvi_col_l5_l7 = (
        landsat_col
            .filter(ee.Filter.inList('SPACECRAFT_ID', ['LANDSAT_5', 'LANDSAT_7']))
            .map(cloudMask)
            .map(srNDVI)
            .map(dateProperties)
            .map(ndvi_etm_2_oli_func)
    )
    ndvi_col_l8_l9 = (
        landsat_col
            .filter(ee.Filter.inList('SPACECRAFT_ID', ['LANDSAT_8', 'LANDSAT_9']))
            .map(cloudMask)
            .map(srNDVI)
            .map(dateProperties)
    )
    ndvi_col = ee.ImageCollection(ndvi_col_l5_l7.merge(ndvi_col_l8_l9))
    
    # CRS and Transform for export
    projection = ndvi_col.first().projection()
    ee_json = ee.Dictionary(ee.Algorithms.Describe(projection))
    ee_crs = ee_json.get('crs')
    ee_transform = ee.List(ee_json.get('transform'))
    
    # final collection of NDVI composites and ET* images 
    final_col = (
        ee.ImageCollection.fromImages(yrs                       
            .map(buildNDVIComp))
            .map(buildETStar)
    )

    # gridMET feature collection w/ average WY values
    grid_wy_col = ee.ImageCollection.fromImages(yrs.map(buildGridMetWY))

    grid_stats = ee.FeatureCollection(grid_wy_col.map(calcAvgGridMet))

    if img_storage_location == 'ee_asset':
        logging.info(f"\nExporting annual Beamer-Minor images to Earth Engine asset {output_path_n}")
    elif img_storage_location == 'cloud_storage':
        logging.info(f"\nExporting annual Beamer-Minor images to Google Cloud Storage Bucket: {output_path_n}")
    logging.info(f"\nPhreatophyte boundary for gridMET cell selection: {phreat_nm}")
    logging.info(f"\nHA boundary: {ha}")
    logging.info(f"\nStart Year: {start_yr}")
    logging.info(f"\nEnd Year: {end_yr}")
    logging.info(f"\nStart DOY: {start_doy}")
    logging.info(f"\nEnd DOY: {end_doy}")
    logging.info(f"\nPath List: {', '.join([str(p) for p in path_lst])}")
    logging.info(f"\nRow List: {', '.join([str(r) for r in row_lst])}")
    if img_storage_location == 'ee_asset':
        logging.info(f"\nAsset ID: {output_path_n}")
    elif img_storage_location == 'cloud_storage':
        logging.info(f"\nCOGs location: {output_path_n}")
    
    # loop through years and export final image w/ all beamer minor bands
    for yr in yr_lst:

        # ndvi composite and ET* image bands for the year
        ndvi_et_img = (
            ee.Image(final_col
                .filter(ee.Filter.eq('year', yr))
                .first())
        )

        # filter to the current year avg value
        grid_c = (
            grid_stats
                .filter(ee.Filter.eq('year', yr))
                .first()
        )

        # convert the average gridMET values into an image (separate bands, constants)
        grid_avg_img = (
            grid_c
                .toDictionary()
                .toImage(['eto_mean', 'pr_mean'])
                .toShort() # signed 16-bit integer
                .updateMask(ndvi_et_img.select(['ndvi']))
        )

        # add gridMET images to ndvi and et_star image
        merged_img = (
            ndvi_et_img
                .addBands(grid_avg_img)
                .rename(['ndvi', 'et_star', 'et_star_lci', 'et_star_uci', 'eto', 'pr'])
        )

        # calculate ETa and ETg bands, along with all CI bands
        out_img = (
            merged_img
                .addBands((merged_img.select('pr').add(merged_img.select('et_star').multiply(merged_img.select('eto').subtract(merged_img.select('pr'))))).toShort()) # signed 16-bit integer
                .addBands((merged_img.select('pr').add(merged_img.select('et_star_lci').multiply(merged_img.select('eto').subtract(merged_img.select('pr'))))).toShort()) # signed 16-bit integer
                .addBands((merged_img.select('pr').add(merged_img.select('et_star_uci').multiply(merged_img.select('eto').subtract(merged_img.select('pr'))))).toShort()) # signed 16-bit integer
                .addBands((merged_img.select('et_star').multiply(merged_img.select('eto').subtract(merged_img.select('pr')))).toShort()) # signed 16-bit integer
                .addBands((merged_img.select('et_star_lci').multiply(merged_img.select('eto').subtract(merged_img.select('pr')))).toShort()) # signed 16-bit integer
                .addBands((merged_img.select('et_star_uci').multiply(merged_img.select('eto').subtract(merged_img.select('pr')))).toShort()) # signed 16-bit integer
                .rename(['ndvi', 'et_star', 'et_star_lci', 'et_star_uci', 'eto', 'pr', 'eta', 'eta_lci', 'eta_uci', 'etg', 'etg_lci', 'etg_uci'])
                .select(['ndvi', 'et_star', 'eto', 'pr', 'eta', 'eta_lci', 'eta_uci', 'etg', 'etg_lci', 'etg_uci'])
                .set({
                    'bmm_regression_version': bmm_r_vers,
                    'bmm_regression_equation': bmm_eqn,
                    'date_ingested': datetime.datetime.now().strftime('%Y-%m-%d'),
                    'gridmet_cell_count': grid_c.get('eto_count'),
                    'gridmet_buffer_flag': g_buf,
                    'gridmet_cell_skip_list': g_skip_lst,
                    'mndwi_mask_source': ndwi_id,
                    'minimum_etg_threshold': min_etg_thresh,
                })
        )

        # set floor on etg values using the user-specified threshold. Default is 0.01 ft (approx. 3.048 mm)
        # leaving these adjusted bands as floats b/c of low values
        out_img = (
            out_img
                .addBands(out_img.select('etg').where(out_img.select('etg').lt(min_etg_thresh_mm), min_etg_thresh_mm))
                .addBands(out_img.select('etg_lci').where(out_img.select('etg_lci').lt(min_etg_thresh_mm), min_etg_thresh_mm))
                .addBands(out_img.select('etg_uci').where(out_img.select('etg_uci').lt(min_etg_thresh_mm), min_etg_thresh_mm))
                .rename(['ndvi', 'et_star', 'eto', 'pr', 'eta', 'eta_lci', 'eta_uci', 'etg', 'etg_lci', 'etg_uci', 'etg_adj', 'etg_lci_adj', 'etg_uci_adj'])
        )

        # mask water pixels out of the final ETg image 
        out_img_masked = out_img.updateMask(ndwi_mask)

        if img_storage_location == 'ee_asset':
            # export task
            out_task = ee.batch.Export.image.toAsset(**{
                'image': out_img_masked,
                'description': f'beamerminor{ha}{yr}',
                'assetId': f'{output_path_n}/{ha}_{yr}',
                'crs': ee_crs,
                'crsTransform': ee_transform,
                # 'scale': 30,
                # 'region': phreat_col.geometry(),
                'region': ha_ftr.geometry(),
                'maxPixels': 1e13,
            })
            out_task.start()
            
        elif img_storage_location == 'cloud_storage':

            # output location path split
            out_split = output_path_n.split('/')
            
            # export task
            out_task = ee.batch.Export.image.toCloudStorage(**{
                'image': out_img_masked.toFloat(),
                'description': f'beamerminor{ha}{yr}_etm2oli',
                'crs': ee_crs,
                'crsTransform': ee_transform,
                'bucket': out_split[0],
                'fileNamePrefix': f'{out_split[1]}/{out_split[2]}/bmm_HA_{out_split[2]}_etg_{yr}',
                'fileFormat': 'GeoTIFF',
                'formatOptions': {
                    'noData': -9999,
                    'cloudOptimized': True,
                },
                # 'scale': 30,
                # 'region': phreat_col.geometry(),
                'region': ha_ftr.geometry(),
                'maxPixels': 1e13,
            })
            out_task.start()

        else:
            logging.info(f"\nincorrect image storage location defined in the config file, cannot export images")
        
        logging.info(f"\ntask started for HA: {ha}, {yr}")
    

def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Earth Engine Beamer-Minor Method image collection prep',
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
