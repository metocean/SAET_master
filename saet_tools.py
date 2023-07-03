'''
SAET (Shoreline Analysis and Extraction Tool). V2.0
Author: Geo-Environmental Cartography and Remote Sensing Group (CGAT).
Polytechnic University of Valencia (Spain)

This module contains all functions needed to run the shoreline extraction
algorithm from Landsat-8/9 and Sentinel-2 satellite imagery.

The module is structured into twelve sections:

- Imports
- Parameters functions.
- High level run functions.
- Download functions.
- Image processing functions.
- Helper functions
- Folder management functions
- Time functions
- Shoreline extraction functions
- Average points functions
- Shapefile creation functions
- Point cleaning functions

'''

# *******************************************************************************
# SECTION: IMPORTS
# *******************************************************************************

# common imports
import json
import shutil
import zipfile
import tarfile
import os
import sys
import pathlib
import time
import datetime
import fnmatch
import numpy as np
from math import exp, sqrt, degrees, atan2, ceil, sin
import gc
import csv
from tqdm import tqdm
import xml.etree.ElementTree as ET
import requests
import re
from collections import OrderedDict
import webbrowser
import logging
import warnings

# specific imports
import shapefile
from osgeo import gdal, osr, ogr
from sentinelsat import read_geojson, geojson_to_wkt
from sentinelsat import SentinelAPI
from landsatxplore2 import api
from landsatxplore2.earthexplorer import EarthExplorer
# from landsatxplore import api
# from landsatxplore.earthexplorer import EarthExplorer
import skimage.morphology as morphology
from skimage.transform import rescale
from skimage.filters import threshold_multiotsu, threshold_otsu
from sklearn.cluster import KMeans
from skimage import exposure
from polynomial import polyfit2d, polyval2d
from matplotlib.pyplot import clf, contour
from scipy.spatial import Delaunay
import networkx as nx

# some warnings related with maptplotlib.contours function will be deactivated
warnings.filterwarnings('ignore')

# *******************************************************************************
# SECTION: PARAMETERS FUNCTION
# *******************************************************************************


def checkCommandArgumments(args):
    '''
    Description:
    ------------
    Check the suitable command line parameters format to avoid input errors.
    Converts footprint parameter into suitable geometry in wkt format (needed to
    searching for images).

    Arguments:
    ------------
    - argparse object


    Returns:
    ------------
    - Parameters in suitable format (string or int)

    '''

    # check integrity for parameters
    if '=' in args.lp:
        l8_product = args.lp.split('=')[1]
    else:
        l8_product = args.lp

    if '=' in args.sp:
        s2_product = args.sp.split('=')[1]
    else:
        s2_product = args.sp

    try:
        if '=' in str(args.sd):
            start_date = str(args.sd).split('=')[1]
        else:
            start_date = str(args.sd)
        datetime.datetime.strptime(start_date, '%Y%m%d')
        if '=' in str(args.ed):
            end_date = str(args.ed).split('=')[1]
        else:
            end_date = str(args.ed)
        datetime.datetime.strptime(end_date, '%Y%m%d')
        if '=' in str(args.cd):
            central_date = str(args.cd).split('=')[1]
        else:
            central_date = str(args.cd)
        datetime.datetime.strptime(central_date, '%Y%m%d')
    except:
        # print('Invalid format for start or end date (YYYYMMDD).')
        logging.warning('Invalid format for start or end date (YYYYMMDD).')
        sys.exit(1)

    if not start_date < central_date < end_date:
        # print('Central date must be between start date and end date.')
        logging.warning(
            'Central date must be between start date and end date.')
        sys.exit(1)

    max_cloud_cover = args.mc

    # processing the filter by beach code
    if '=' in args.bc:
        bc = args.sl.split('=')[1]
    else:
        bc = args.bc
    if bc != 'NONE':
        if ',' in bc:
            bc = bc.replace(" ", "")
            res1 = re.findall(r'^\d+(?:[ \t]*,[ \t]*\d+)+$', bc)
            if len(res1) == 1:
                bc = res1[0].split(',')
            else:
                logging.warning('Invalid format for beach code filter.')
                sys.exit(1)
        else:
            res2 = re.findall(r'^\d+$', bc)
            if len(res2) == 1:
                bc = res2
            else:
                logging.warning('Invalid format for beach code filter.')
                sys.exit(1)
    else:
        bc = 'NONE'

    # convert beach code list in string with tuple format
    query = '('
    for i in bc:
        query = query+i+','
    query = query[:-1]
    query = query+')'

    beach_code_filter = query

    # processing footprint for image searching. The final format is different depending on L8 or S2
    scene_landast_list = []
    scene_sentinel_list = []
    if '=' in args.fp:
        footprint_path = args.fp.split('=')[1]
    else:
        footprint_path = args.fp
    if footprint_path != 'NONE':
        if not os.path.isfile(footprint_path):
            res = re.findall(r'-?\d+\.?\d*', footprint_path)
            if len(res) != 4:
                # print('Invalid footprint path or coordinate format.')
                logging.warning('Invalid footprint path or coordinate format.')
                sys.exit(1)
            else:
                xmin = res[0]
                ymin = res[1]
                xmax = res[2]
                ymax = res[3]
                footprint = 'GEOMETRYCOLLECTION(POLYGON(('+xmin+' '+ymin+','+xmax + \
                    ' '+ymin+','+xmax+' '+ymax+','+xmin+' '+ymax+','+xmin+' '+ymin+')))'
        else:
            # reading roi footprint
            footprint = geojson_to_wkt(read_geojson(footprint_path))
    else:
        footprint = 'NONE'
        if '=' in args.ll:
            ll = args.ll.split('=')[1]
        else:
            ll = args.ll
        if ll != 'NONE':
            for id_ll in args.ll.split(','):
                res = re.findall('[0-9]{6}', id_ll)
                if len(res) == 0:
                    # print('Invalid Landsat 8 scene code.')
                    logging.warning('Invalid Landsat scene code.')
                    sys.exit(1)

            ll = args.ll.split(',')  # list of Landsat 8 scenes
        else:
            ll = 'NONE'

        if '=' in args.sl:
            sl = args.sl.split('=')[1]
        else:
            sl = args.sl
        if sl != 'NONE':
            for id_sl in sl.split(','):
                res = re.findall('[0-9]{2}[A-Z]{3}', id_sl)
                if len(res) == 0:
                    # print('Invalid Sentinel 2 scene code.')
                    logging.warning('Invalid Sentinel 2 scene code.')
                    sys.exit(1)

            sl = args.sl.split(',')  # list of sentinel 2 scenes
        else:
            sl = 'NONE'

        scene_landast_list = ll
        scene_sentinel_list = sl

    if '=' in args.rm:
        run_mode = args.rm.split('=')[1]
    else:
        run_mode = args.rm
    if not run_mode in ['os', 'dp', 'od', 'op', 'or']:
        # print('Invalid run mode.')
        logging.warning('Invalid run mode.')
        sys.exit(1)

    # output dictionary
    run_parameters = {}
    run_parameters['l8_product'] = l8_product
    run_parameters['s2_product'] = s2_product
    run_parameters['start_date'] = start_date
    run_parameters['central_date'] = central_date
    run_parameters['end_date'] = end_date
    run_parameters['max_cloud_cover'] = max_cloud_cover
    run_parameters['footprint'] = footprint
    run_parameters['run_mode'] = run_mode
    run_parameters['scene_landsat_list'] = scene_landast_list
    run_parameters['scene_sentinel_list'] = scene_sentinel_list
    run_parameters['beach_code_filter'] = beach_code_filter

    return run_parameters


def parseNumberProducts(np, scenes):
    '''
    Description:
    ------------
    The function filters the scenes to be downloaded and procesed 
    acording with a list of numbers (identifiers)

    Arguments:
    ------------
    - np (list of integers): list with the number of products
    - scenes (list of dictionaries): list of scenes to be filtered

    Returns:
    ------------
    - list of dictionaries: list of filtered scenes

    '''
    #lista = list(range(0, len(scenes)))
    lista = []
    for key in scenes.keys():
        lista.append(int(key))
    filtered_scenes = []
    if np != 'NONE':
        if len(np) == 1 and np != '*':
            try:
                # return [scenes[int(np)]]
                return [scenes[str(np)]]
            except:
                # logging.warning(f'Incorrect number of product: {np}'+'\n')
                # sys.exit(1)
                pass
        if np == '*':
            # return scenes
            return scenes.values()
        if len(np) > 1:
            res = re.findall(r'^\d+(?:[ \t]*,[ \t]*\d+)+$', np)
            res2 = re.findall(r'^[0-9]\d*-[0-9]\d*', np)
            if len(res) != 0:
                n_lista = []
                for n in res[0].split(','):
                    if int(n) in lista:
                        n_lista.append(int(n))
                    else:
                        # logging.warning(
                        #     f'Incorrect number of product: {n}'+'\n')
                        # sys.exit(1)
                        pass
                n_lista_sorted = sorted(list(set(n_lista)))
                for n in n_lista_sorted:
                    # filtered_scenes.append(scenes[n])
                    filtered_scenes.append(scenes[str(n)])
                return filtered_scenes
            if len(res2) != 0:
                lim_inf = int(res2[0].split('-')[0])
                lim_sup = int(res2[0].split('-')[1])
                if lim_inf > lim_sup:
                    logging.warning(
                        f'Incorrect range of products'+'\n')
                    sys.exit(1)
                else:
                    if lim_inf not in lista or lim_sup not in lista:
                        logging.warning(
                            f'Incorrect range of products'+'\n')
                        sys.exit(1)
                    else:
                        n_lista = list(range(lim_inf, lim_sup+1))
                        for n in n_lista:
                            # filtered_scenes.append(scenes[n])
                            filtered_scenes.append(scenes[str(n)])
                        return filtered_scenes
            logging.warning(
                f'Incorrect format for list or range of products'+'\n')
            sys.exit(1)
    return []


def parseNumberProductsActivation(np, scenes):
    '''
    Description:
    ------------
    The function filters the S2 scenes to be activated 
    acording with a list of numbers (identifiers)

    Arguments:
    ------------
    - np (list of integers): list with the number of products
    - scenes (list of dictionaries): list of scenes to be filtered

    Returns:
    ------------
    - list of dictionaries: list of filtered scenes

    '''
    lista = list(range(0, len(scenes)))
    filtered_scenes = []
    if np != 'NONE':
        if len(np) == 1 and np != '*':
            try:
                return [scenes[np]]
            except:
                logging.warning(f'Incorrect number of product: {np}'+'\n')
                sys.exit(1)
        if np == '*':
            return scenes
        if len(np) > 1:
            res = re.findall(r'^\d+(?:[ \t]*,[ \t]*\d+)+$', np)
            res2 = re.findall(r'^[0-9]\d*-[0-9]\d*', np)
            if len(res) != 0:
                n_lista = []
                for n in res[0].split(','):
                    if int(n) in lista:
                        n_lista.append(int(n))
                    else:
                        logging.warning(
                            f'Incorrect number of product: {n}'+'\n')
                        sys.exit(1)
                n_lista_sorted = sorted(list(set(n_lista)))
                for n in n_lista_sorted:
                    filtered_scenes.append(scenes[str(n)])
                return filtered_scenes
            if len(res2) != 0:
                lim_inf = int(res2[0].split('-')[0])
                lim_sup = int(res2[0].split('-')[1])
                if lim_inf > lim_sup:
                    logging.warning(
                        f'Incorrect range of products'+'\n')
                    sys.exit(1)
                else:
                    if lim_inf not in lista or lim_sup not in lista:
                        logging.warning(
                            f'Incorrect range of products'+'\n')
                        sys.exit(1)
                    else:
                        n_lista = list(range(lim_inf, lim_sup+1))
                        for n in n_lista:
                            filtered_scenes.append(scenes[str(n)])
                        if len(filtered_scenes) > 20:
                            logging.warning(
                                f'Maximum number of offline images activation at the same time is 20. You have ordered {len(filtered_scenes)}'+'\n')
                            sys.exit(1)
                        return filtered_scenes
            logging.warning(
                f'Incorrect format for list or range of products'+'\n')
            sys.exit(1)
    return []


def parseNumberProductsFolder(np, scenes):
    '''
    Description:
    ------------
    The function filters the scenes from the output data to be processed 
    acording with a list of numbers (identifiers)

    Arguments:
    ------------
    - np (list of integers): list with the number of products
    - scenes (list of dictionaries): list of scenes to be filtered

    Returns:
    ------------
    - list of dictionaries: list of filtered scenes

    '''
    lista = list(range(0, len(scenes)))
    filtered_scenes = []
    if np != 'NONE':
        if len(np) == 1 and np != '*':
            try:
                return [scenes[int(np)]]
            except:
                logging.warning(f'Incorrect number of product: {np}'+'\n')
                sys.exit(1)
        if np == '*':
            return scenes
        if len(np) > 1:
            res = re.findall(r'^\d+(?:[ \t]*,[ \t]*\d+)+$', np)
            res2 = re.findall(r'^[0-9]\d*-[0-9]\d*', np)
            if len(res) != 0:
                n_lista = []
                for n in res[0].split(','):
                    if int(n) in lista:
                        n_lista.append(int(n))
                    else:
                        logging.warning(
                            f'Incorrect number of product: {n}'+'\n')
                        sys.exit(1)
                n_lista_sorted = sorted(list(set(n_lista)))
                for n in n_lista_sorted:
                    filtered_scenes.append(scenes[n])
                return filtered_scenes
            if len(res2) != 0:
                lim_inf = int(res2[0].split('-')[0])
                lim_sup = int(res2[0].split('-')[1])
                if lim_inf > lim_sup:
                    logging.warning(
                        f'Incorrect range of products'+'\n')
                    sys.exit(1)
                else:
                    if lim_inf not in lista or lim_sup not in lista:
                        logging.warning(
                            f'Incorrect range of products'+'\n')
                        sys.exit(1)
                    else:
                        n_lista = list(range(lim_inf, lim_sup+1))
                        for n in n_lista:
                            filtered_scenes.append(scenes[n])
                        return filtered_scenes
            logging.warning(
                f'Incorrect format for list or range of products'+'\n')
            sys.exit(1)
    return []


def readAuthenticationParameters():
    '''
    Description:
    ------------
    Read authentication parameters for ESA scihub and USGS servers.
    The file in json format containing access credentials must be in the root folder and its name is 'auth.json'

    Arguments:
    ------------
    None

    Returns:
    ------------
    - ESA and USGS users and passwords to access to their servers (strings)

    '''
    with open('auth.json') as json_file:
        data = json.load(json_file)

        for p in data['authentication']:
            user_esa = p['user_esa']
            password_esa = p['password_esa']
            user_usgs = p['user_usgs']
            password_usgs = p['password_usgs']

    return user_esa, password_esa, user_usgs, password_usgs


def esaAuthentication(user_esa, pass_esa):
    '''
    Description:
    ------------
    Access to the ESA scihub api using credentials from
    readAuthenticationParameters() function

    Arguments:
    ------------
    None

    Returns:
    ------------
    - ESA scihub api object

    '''

    try:
        ss_api = SentinelAPI(user_esa, pass_esa)
        return ss_api
    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def usgsAuthentication(user_usgs, pass_usgs):
    '''
    Description:
    ------------
    Access to the USGS apis (seraching and download) using credentials from
    readAuthenticationParameters() function

    Arguments:
    ------------
    - user_usgs(string): user credential
    - pass_usgs(string): password

    Returns:
    ------------
    - USGS search and download api objects

    '''

    try:
        usgs_api_search = api.API(user_usgs, pass_usgs)
        usgs_api_download = EarthExplorer(user_usgs, pass_usgs)
        return usgs_api_search, usgs_api_download
    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


# *******************************************************************************
# SECTION: HIGH LEVEL RUN FUNCTIONS
# *******************************************************************************

def processLandsat8Scenes(l8_scenes, data_path, sds_path, shp_path, user_usgs, pass_usgs, run_mode, wi_type, thr_method, cloud_mask_level, morphology_method, kernel_size, bc):
    '''
    Description:
    ------------
    Run the whole workflow for downloading and processing Landsat-8 scenes.
    For every scene, a zip file with all bands is downloaded.

    Arguments:
    ------------
    - l8_scenes: dictionary of scenes provided for  startSearchForLandsat8() function
    - data_path (string): folder for store every downloaded and processing image. Obtained from createFolderTree() function
    - sds_path (string): folder for store every extracted shoreline. Obtained from createFolderTree() function
    - shp_path (string): path to shapefile containing beaches areas
    - user_usgs(string): usgs user
    - pass_usgs (string): usgs password
    - run_mode (string): run mode
    - wi_type (string): type of water index
    - thr_method (string); type of thresholding method
    - cloud_mask_level (string): cloud masking severity
    - morphology_method (string): type of morphology method
    - kernel_size (int): kernel size for the extraction algorithm
    - bc (list of string): list of beach polygons to be processed

    Returns:
    ------------
    None

    '''
    # download and processing every scene
    for l8_scene in l8_scenes:
        scene_path, scene_id = downloadLandsat8_2(
            l8_scene, data_path, user_usgs, pass_usgs)

        # print('Processing '+scene_path+' ...')
        logging.info('Processing '+scene_id+' ...')
        processing_path = str(pathlib.Path(
            os.path.join(scene_path, 'temp')))
        # print('Computing water index band...')
        logging.info('Computing water index band...')
        if wi_type == 'aweinsh':
            aweinshLandsat(scene_path)
        if wi_type == 'aweish':
            aweishLandsat(scene_path)
        if wi_type == 'mndwi':
            mndwiLandsat(scene_path)
        if wi_type == 'kmeans':
            computeKmeansLandsat(scene_path)
        # print('Computing cloud mask...')
        logging.info('Computing cloud mask...')
        createCloudMaskL8(scene_path, cloud_mask_level)
        logging.info('Computing water index mask...')
        if wi_type != 'kmeans':
            wi_path = getBandPath(processing_path, 'wi.tif')
            wi_mask = getIndexMask(wi_path, thr_method)
        else:
            wi_path = getBandPath(processing_path, 'kmeans_mask.tif')
            wi_mask = getBandData(wi_path)
        cmask_band = getBandPath(processing_path, 'cmask')
        # print('Computing rough pixel line...')
        logging.info('Computing rough pixel line...')
        pixel_line = createPixelLine(
            morphology_method, wi_mask, cmask_band)
        saveMask(pixel_line, str(pathlib.Path(os.path.join(
            processing_path, 'pl.tif'))), cmask_band)
        source_epsg = getSourceEpsg()
        target_epsg = getTargetEpsg(scene_path, 'B6')
        # print('Reprojecting shp of beaches...')
        logging.info('Reprojecting shp of beaches...')
        reprojectShp(shp_path, str(pathlib.Path(os.path.join(processing_path,
                                                             'bb300_r.shp'))), source_epsg, target_epsg)
        # print('Computing footprint band...')
        logging.info('Computing footprint band...')
        createShapefileFromRasterFootprint(getBandPath(scene_path, 'B6'), str(pathlib.Path(os.path.join(
            processing_path, 'scene_footprint.shp'))), target_epsg, geom_type='polygon')
        # print('Clipping shp of beaches by scene footprint...')
        logging.info('Clipping shp of beaches by scene footprint...')
        clipShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
            processing_path, 'clip_bb300_r.shp'))), str(pathlib.Path(os.path.join(processing_path, 'scene_footprint.shp'))))
        # print('Rasterizing beaches subset...')
        logging.info('Rasterizing beaches subset...')
        rasterizeShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
            processing_path, 'bb300_r.tif'))), getBandPath(scene_path, 'B6'), bc)
        # print('Masking rough pixel line with beaches subset...')
        logging.info('Masking rough pixel line with beaches subset...')
        maskPixelLine(str(pathlib.Path(os.path.join(processing_path, 'pl.tif'))),
                      str(pathlib.Path(os.path.join(processing_path, 'bb300_r.tif'))))
        # print('Extracting points...')
        logging.info('Extracting points...')
        res = extractPoints(getBandPath(scene_path, 'B6'), str(pathlib.Path(os.path.join(
            processing_path, 'pl.tif'))), processing_path, int(kernel_size), 4, 3)
        if res:
            # print('Computing average points...')
            logging.info('Computing average points...')
            averagePoints(getBandPath(scene_path, 'B6'),
                          processing_path, 50, 3)
            # print('Making point shp...')
            logging.info('Making point shp...')
            shp_path_average = createShpFromAverageFile(
                getBandPath(scene_path, 'B6'), processing_path)
            # print('Transfering beaches identifiers...')
            logging.info('Transfering beaches identifiers...')
            copyShpIdentifiers(str(pathlib.Path(os.path.join(
                processing_path, 'clip_bb300_r.shp'))), shp_path_average)
            # print('Cleaning points and making final shoreline in line vector format...')
            logging.info(
                'Cleaning points and making final shoreline in line vector format...')
            cleanPoints2(shp_path_average, 150, 1)
            # print('Export final shoreline shapefiles to SDS folder...')
            logging.info(
                'Export final shoreline shapefiles to SDS folder...')
            copyShpToFolder(processing_path, sds_path, target_epsg)
        else:
            logging.warning('No results in extraction points process.')
            sys.exit(1)


def processLandsatSceneFromPath(scene_path, scene_id, sds_path, shp_path, wi_type, thr_method, cloud_mask_level, morphology_method, kernel_size, bc):
    '''
    Description:
    ------------
    Run the whole workflow reprocessing Landsat-8 scenes from the data folder.

    Arguments:
    ------------
    - scene_path (string): folder that contains the image to be processed
    - scene_id (string): name of the image (id)
    - sds_path (string): folder for store every extracted shoreline. Obtained from createFolderTree() function
    - shp_path (string): path to shapefile containing beaches areas
    - wi_type (string): type of water index
    - thr_method (string); type of thresholding method
    - cloud_mask_level (string): cloud masking severity
    - morphology_method (string): type of morphology method
    - kernel_size (int): kernel size for the extraction algorithm
    - bc (list of string): list of beach polygons to be processed

    Returns:
    ------------
    None

    '''

    # print('Processing '+scene_path+' ...')
    logging.info('Processing '+scene_id+' ...')
    processing_path = str(pathlib.Path(
        os.path.join(scene_path, 'temp')))
    # print('Computing water index band...')
    logging.info('Computing water index band...')
    if wi_type == 'aweinsh':
        aweinshLandsat(scene_path)
    if wi_type == 'aweish':
        aweishLandsat(scene_path)
    if wi_type == 'mndwi':
        mndwiLandsat(scene_path)
    if wi_type == 'kmeans':
        computeKmeansLandsat(scene_path)
    # print('Computing cloud mask...')
    logging.info('Computing cloud mask...')
    createCloudMaskL8(scene_path, cloud_mask_level)
    logging.info('Computing water index mask...')
    if wi_type != 'kmeans':
        wi_path = getBandPath(processing_path, 'wi.tif')
        wi_mask = getIndexMask(wi_path, thr_method)
    else:
        wi_path = getBandPath(processing_path, 'kmeans_mask.tif')
        wi_mask = getBandData(wi_path)
    cmask_band = getBandPath(processing_path, 'cmask')
    # print('Computing rough pixel line...')
    logging.info('Computing rough pixel line...')
    pixel_line = createPixelLine(
        morphology_method, wi_mask, cmask_band)
    saveMask(pixel_line, str(pathlib.Path(os.path.join(
        processing_path, 'pl.tif'))), cmask_band)
    source_epsg = getSourceEpsg()
    target_epsg = getTargetEpsg(scene_path, 'B6')
    # print('Reprojecting shp of beaches...')
    logging.info('Reprojecting shp of beaches...')
    reprojectShp(shp_path, str(pathlib.Path(os.path.join(processing_path,
                                                         'bb300_r.shp'))), source_epsg, target_epsg)
    # print('Computing footprint band...')
    logging.info('Computing footprint band...')
    createShapefileFromRasterFootprint(getBandPath(scene_path, 'B6'), str(pathlib.Path(os.path.join(
        processing_path, 'scene_footprint.shp'))), target_epsg, geom_type='polygon')
    # print('Clipping shp of beaches by scene footprint...')
    logging.info('Clipping shp of beaches by scene footprint...')
    clipShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
        processing_path, 'clip_bb300_r.shp'))), str(pathlib.Path(os.path.join(processing_path, 'scene_footprint.shp'))))
    # print('Rasterizing beaches subset...')
    logging.info('Rasterizing beaches subset...')
    rasterizeShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
        processing_path, 'bb300_r.tif'))), getBandPath(scene_path, 'B6'), bc)
    # print('Masking rough pixel line with beaches subset...')
    logging.info('Masking rough pixel line with beaches subset...')
    maskPixelLine(str(pathlib.Path(os.path.join(processing_path, 'pl.tif'))),
                  str(pathlib.Path(os.path.join(processing_path, 'bb300_r.tif'))))
    # print('Extracting points...')
    logging.info('Extracting points...')
    res = extractPoints(getBandPath(scene_path, 'B6'), str(pathlib.Path(os.path.join(
        processing_path, 'pl.tif'))), processing_path, int(kernel_size), 4, 3)
    if res:
        # print('Computing average points...')
        logging.info('Computing average points...')
        averagePoints(getBandPath(scene_path, 'B6'),
                      processing_path, 50, 3)
        # print('Making point shp...')
        logging.info('Making point shp...')
        shp_path_average = createShpFromAverageFile(
            getBandPath(scene_path, 'B6'), processing_path)
        # print('Transfering beaches identifiers...')
        logging.info('Transfering beaches identifiers...')
        copyShpIdentifiers(str(pathlib.Path(os.path.join(
            processing_path, 'clip_bb300_r.shp'))), shp_path_average)
        # print('Cleaning points and making final shoreline in line vector format...')
        logging.info(
            'Cleaning points and making final shoreline in line vector format...')
        cleanPoints2(shp_path_average, 150, 1)
        # print('Export final shoreline shapefiles to SDS folder...')
        logging.info(
            'Export final shoreline shapefiles to SDS folder...')
        copyShpToFolder(processing_path, sds_path, target_epsg)
    else:
        logging.warning('No results in extraction points process.')
        sys.exit(1)


def offlineS2Activation(scenes, run_parameters):
    '''
    Description:
    ------------
    Triggers the retrieval of offline Sentinel-2 images. The user will select the number of
    images to be retrieved. In order to check if the image is online after triggering, use the
    parameter --oa=check

    Arguments:
    ------------
    - scenes: dictionary of the S2 scenes provided for  startSearchSentinel2() function
    - run_parameters (list of strings): list of parameters.

    Returns:
    ------------
    None

    '''
    user_esa = run_parameters['user_esa']
    pass_esa = run_parameters['pass_esa']
    esa_api = esaAuthentication(user_esa, pass_esa)
    print('')
    np = input('Number of S2 images to be activated?: ')
    filtered_scenes = parseNumberProductsActivation(np, scenes)
    if len(filtered_scenes) > 0:
        for scene in filtered_scenes:
            product_info = esa_api.get_product_odata(scene, full=True)
            scene_title = product_info['title']
            img_id = product_info['id']
            time.sleep(2)
            try:
                esa_api.trigger_offline_retrieval(img_id)
                logging.info(scene_title+' -> Offline retrieval successful...')
                #print(scene_title, 'Offline retrieval successful...')
            except:
                logging.info(scene_title+' -> Offline retrieval failed...')
                #print(scene_title, 'Offline retrieval failed...')
    else:
        logging.warning('There are no products S2 to be activated.')
        sys.exit(1)


def processSentinel2Scenes(product_type, s2_scenes, data_path, sds_path, shp_path, user_esa, pass_esa, run_mode, s2_titles, wi_type, thr_method, cloud_mask_level, morphology_method, kernel_size, bc):
    '''
    Description:
    ------------
    Run the whole workflow for downloading and processing Sentinel-2 scenes.
    For every scene, just only needed bands, according with the type of product, are downloaded.

    Arguments:
    ------------
    - product_type (string): type of product for S2 images. Can be S2MSI1C or S2MSI2A
    - s2_scenes: dictionary of scenes provided for  startSearchSentinel2() function
    - data_path (string): folder for store every downloaded and processing image. Obtained from createFolderTree() function
    - sds_path (string): folder for store every extracted shoreline. Obtained from createFolderTree() function
    - shp_path (string): path to shapefile containing beaches areas
    - user_esa(string): esa user
    - pass_esa (string): esa password
    - run_mode (string): run mode
    - s2_titles (list of strings): list with the titles for each S2 image
    - wi_type (string): type of water index
    - thr_method (string); type of thresholding method
    - cloud_mask_level (string): cloud masking severity
    - morphology_method (string): type of morphology method
    - kernel_size (int): kernel size for the extraction algorithm
    - bc (list of string): list of beach polygons to be processed

    Returns:
    ------------
    None

    '''

    # selection of suitable bands to be downloaded
    if product_type == 'S2MSI2A':
        bands = ['B02', 'B03', 'B08', 'B11', 'B12', 'SCL']

    if product_type == 'S2MSI1C':
        bands = ['B02', 'B03', 'B08', 'B11', 'B12', 'QA60', 'CPM']

    # processing of every scene
    for s2_scene in s2_scenes:
        scene_path, title = downloadS2_2(
            product_type, s2_scene, bands, data_path, user_esa, pass_esa)

        # print('Processing '+title+' ...')
        logging.info('Processing '+title+' ...')
        processing_path = str(pathlib.Path(
            os.path.join(scene_path, 'temp')))
        # print('Computing water index band...')
        logging.info('Computing water index band...')
        if wi_type == 'aweinsh':
            aweinshS2(scene_path)
        if wi_type == 'aweish':
            aweishS2(scene_path)
        if wi_type == 'mndwi':
            mndwiS2(scene_path)
        if wi_type == 'kmeans':
            computeKmeansS2(scene_path)
        # print('Computing cloud mask...')
        logging.info('Computing cloud mask...')
        createCloudMaskS2(scene_path, cloud_mask_level)
        # print('Computing water index mask...')
        logging.info('Computing water index mask...')
        if wi_type != 'kmeans':
            wi_path = getBandPath(processing_path, 'wi.tif')
            wi_mask = getIndexMask(wi_path, thr_method)
        else:
            wi_path = getBandPath(processing_path, 'kmeans_mask.tif')
            wi_mask = getBandData(wi_path)
        cmask_band = getBandPath(processing_path, 'cmask.tif')
        # print('Computing rough pixel line...')
        logging.info('Computing rough pixel line...')
        pixel_line = createPixelLine(
            morphology_method, wi_mask, cmask_band)
        saveMask(pixel_line, str(pathlib.Path(os.path.join(
            processing_path, 'pl.tif'))), cmask_band)
        source_epsg = getSourceEpsg()
        target_epsg = getTargetEpsg(scene_path, 'B11')
        # print('Reprojecting shp of beaches...')
        logging.info('Reprojecting shp of beaches...')
        reprojectShp(shp_path, str(pathlib.Path(os.path.join(processing_path,
                                                             'bb300_r.shp'))), source_epsg, target_epsg)
        # print('Computing footprint band...')
        logging.info('Computing footprint band...')
        createShapefileFromRasterFootprint(getBandPath(scene_path, 'B11'), str(pathlib.Path(os.path.join(
            processing_path, 'scene_footprint.shp'))), target_epsg, geom_type='polygon')
        # print('Clipping shp of beaches by scene footprint...')
        logging.info('Clipping shp of beaches by scene footprint...')
        clipShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
            processing_path, 'clip_bb300_r.shp'))), str(pathlib.Path(os.path.join(processing_path, 'scene_footprint.shp'))))
        # print('Rasterizing beaches subset...')
        logging.info('Rasterizing beaches subset...')
        rasterizeShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
            processing_path, 'bb300_r.tif'))), getBandPath(scene_path, 'B11'), bc)
        # print('Masking rough pixel line with beaches subset...')
        logging.info('Masking rough pixel line with beaches subset...')
        maskPixelLine(str(pathlib.Path(os.path.join(processing_path, 'pl.tif'))),
                      str(pathlib.Path(os.path.join(processing_path, 'bb300_r.tif'))))
        # print('Extracting points...')
        logging.info('Extracting points...')
        res = extractPoints(getBandPath(scene_path, 'B11'), str(pathlib.Path(os.path.join(
            processing_path, 'pl.tif'))), processing_path, int(kernel_size), 4, 3)
        if res:
            # print('Computing average points...')
            logging.info('Computing average points...')
            averagePoints(getBandPath(scene_path, 'B11'),
                          processing_path, 50, 3)
            # print('Making point shp...')
            logging.info('Making point shp...')
            shp_path_average = createShpFromAverageFile(
                getBandPath(scene_path, 'B11'), processing_path)
            # print('Transfering beaches identifiers...')
            logging.info('Transfering beaches identifiers...')
            copyShpIdentifiers(str(pathlib.Path(os.path.join(
                processing_path, 'clip_bb300_r.shp'))), shp_path_average)
            # print('Cleaning points and making final shoreline in line vector format...')
            logging.info(
                'Cleaning points and making final shoreline in line vector format...')
            cleanPoints2(shp_path_average, 150, 1)
            # print('Export final shoreline shapefiles to SDS folder...')
            logging.info(
                'Export final shoreline shapefiles to SDS folder...')
            copyShpToFolder(processing_path, sds_path, target_epsg)
        else:
            logging.warning('No results in extraction points process.')
            sys.exit(1)


def processSentinel2SceneFromPath(scene_path, title, sds_path, shp_path, wi_type, thr_method, cloud_mask_level, morphology_method, kernel_size, bc):
    '''
    Description:
    ------------
    Run the whole workflow reprocessing Sentinel-2 scenes from a path.

    Arguments:
    ------------

    - scene_path (string): path to the folder of the Sentinel-2 image
    - title (string): identifier of the image (coincident with the name of the folder)
    - shp_path (string): path to shapefile containing beaches areas
    - wi_type (string): type of water index
    - thr_method (string); type of thresholding method
    - cloud_mask_level (string): cloud masking severity
    - morphology_method (string): type of morphology method
    - kernel_size (int): kernel size for the extraction algorithm
    - bc (list of string): list of beach polygons to be processed

    Returns:
    ------------
    None

    '''
    #print('Processing '+title+' ...')
    logging.info('Processing '+title+' ...')
    processing_path = str(pathlib.Path(
        os.path.join(scene_path, 'temp')))
    #print('Computing water index band...')
    logging.info('Computing water index band...')
    if wi_type == 'aweinsh':
        aweinshS2(scene_path)
    if wi_type == 'aweish':
        aweishS2(scene_path)
    if wi_type == 'mndwi':
        mndwiS2(scene_path)
    if wi_type == 'kmeans':
        computeKmeansS2(scene_path)
    #print('Computing cloud mask...')
    logging.info('Computing cloud mask...')
    createCloudMaskS2(scene_path, cloud_mask_level)
    #print('Computing water index mask...')
    logging.info('Computing water index mask...')
    if wi_type != 'kmeans':
        wi_path = getBandPath(processing_path, 'wi.tif')
        wi_mask = getIndexMask(wi_path, thr_method)
    else:
        wi_path = getBandPath(processing_path, 'kmeans_mask.tif')
        wi_mask = getBandData(wi_path)
    cmask_band = getBandPath(processing_path, 'cmask.tif')
    #print('Computing rough pixel line...')
    logging.info('Computing rough pixel line...')
    pixel_line = createPixelLine(
        morphology_method, wi_mask, cmask_band)
    saveMask(pixel_line, str(pathlib.Path(os.path.join(
        processing_path, 'pl.tif'))), cmask_band)
    source_epsg = getSourceEpsg()
    target_epsg = getTargetEpsg(scene_path, 'B11')
    #print('Reprojecting shp of beaches...')
    logging.info('Reprojecting shp of beaches...')
    reprojectShp(shp_path, str(pathlib.Path(os.path.join(processing_path,
                                                         'bb300_r.shp'))), source_epsg, target_epsg)
    #print('Computing footprint band...')
    logging.info('Computing footprint band...')
    createShapefileFromRasterFootprint(getBandPath(scene_path, 'B11'), str(pathlib.Path(os.path.join(
        processing_path, 'scene_footprint.shp'))), target_epsg, geom_type='polygon')
    #print('Clipping shp of beaches by scene footprint...')
    logging.info('Clipping shp of beaches by scene footprint...')
    clipShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
        processing_path, 'clip_bb300_r.shp'))), str(pathlib.Path(os.path.join(processing_path, 'scene_footprint.shp'))))
    #print('Rasterizing beaches subset...')
    logging.info('Rasterizing beaches subset...')
    rasterizeShapefile(str(pathlib.Path(os.path.join(processing_path, 'bb300_r.shp'))), str(pathlib.Path(os.path.join(
        processing_path, 'bb300_r.tif'))), getBandPath(scene_path, 'B11'), bc)
    #print('Masking rough pixel line with beaches subset...')
    logging.info('Masking rough pixel line with beaches subset...')
    maskPixelLine(str(pathlib.Path(os.path.join(processing_path, 'pl.tif'))),
                  str(pathlib.Path(os.path.join(processing_path, 'bb300_r.tif'))))
    #print('Extracting points...')
    logging.info('Extracting points...')
    res = extractPoints(getBandPath(scene_path, 'B11'), str(pathlib.Path(os.path.join(
        processing_path, 'pl.tif'))), processing_path, int(kernel_size), 4, 3)
    if res:
        #print('Computing average points...')
        logging.info('Computing average points...')
        averagePoints(getBandPath(scene_path, 'B11'),
                      processing_path, 50, 3)
        #print('Making point shp...')
        logging.info('Making point shp...')
        shp_path_average = createShpFromAverageFile(
            getBandPath(scene_path, 'B11'), processing_path)
        #print('Transfering beaches identifiers...')
        logging.info('Transfering beaches identifiers...')
        copyShpIdentifiers(str(pathlib.Path(os.path.join(
            processing_path, 'clip_bb300_r.shp'))), shp_path_average)
        #print('Cleaning points and making final shoreline in line vector format...')
        logging.info(
            'Cleaning points and making final shoreline in line vector format...')
        cleanPoints2(shp_path_average, 150, 1)
        #print('Export final shoreline shapefiles to SDS folder...')
        logging.info(
            'Export final shoreline shapefiles to SDS folder...')
        copyShpToFolder(processing_path, sds_path, target_epsg)
    else:
        logging.warning('No results in extraction points process.')
        sys.exit(1)


# *******************************************************************************
# SECTION: DOWNLOAD FUNCTIONS
# *******************************************************************************


def downloadS2_2(product_type, scene, bands, base_folder, user_esa, pass_esa):
    '''
    Description:
    ------------
    Download specific bands from one single S2 scene.

    Arguments:
    ------------
    - product_type (string): type of product for S2 images. Can be S2MSI1C or S2MSI2A
    - scene: dictionary of scenes provided for  startSearchSentinel2() function
    - bands: list of bands to be downloaded (list of strings: ['B11','B02','QA60',...])
    - base_folder (string): folder for store every downloaded band. Obtained from createFolderTree() function
    - user_esa (string): user credential for the ESA server
    - pass_esa (string): password credential for the ESA server

    Returns:
    ------------
    - scene_path (string): path to the folder where the bands of the scene will be download

    '''

    # authentication for ESA api
    api = esaAuthentication(user_esa, pass_esa)

    # authentication for url requests with credentials
    api_session = requests.Session()
    api_session.auth = (user_esa, pass_esa)

    # create output folder if it is needed
    createFolderCheck(str(pathlib.Path(os.path.join(base_folder, 's2'))))
    output_folder_s2 = str(pathlib.Path(os.path.join(base_folder, 's2')))

    # band id matching for manifest file in every type of product
    if product_type == 'S2MSI1C':
        # multiple resolution bands
        band_dict = {
            'B01': 'IMG_DATA_Band_60m_1_Tile1_Data',
            'B02': 'IMG_DATA_Band_10m_1_Tile1_Data',
            'B03': 'IMG_DATA_Band_10m_2_Tile1_Data',
            'B04': 'IMG_DATA_Band_10m_3_Tile1_Data',
            'B05': 'IMG_DATA_Band_20m_1_Tile1_Data',
            'B06': 'IMG_DATA_Band_20m_2_Tile1_Data',
            'B07': 'IMG_DATA_Band_20m_3_Tile1_Data',
            'B08': 'IMG_DATA_Band_10m_4_Tile1_Data',
            'B09': 'IMG_DATA_Band_60m_2_Tile1_Data',
            'B10': 'IMG_DATA_Band_60m_3_Tile1_Data',
            'B11': 'IMG_DATA_Band_20m_5_Tile1_Data',
            'B12': 'IMG_DATA_Band_20m_6_Tile1_Data',
            'B8A': 'IMG_DATA_Band_20m_4_Tile1_Data',
            'TCI': 'IMG_DATA_Band_TCI_Tile1_Data',
            'QA60': 'FineCloudMask_Tile1_Data',
            'CPM':  'ClassiPixelsMask_Band_00m_0_Tile1_Data'
        }

    if product_type == 'S2MSI2A':
        # 20 m resolution bands
        band_dict = {
            'B01': 'IMG_DATA_Band_AOT_20m_Tile1_Data',
            'B02': 'IMG_DATA_Band_B02_20m_Tile1_Data',
            'B03': 'IMG_DATA_Band_B03_20m_Tile1_Data',
            'B04': 'IMG_DATA_Band_B04_20m_Tile1_Data',
            'B05': 'IMG_DATA_Band_B05_20m_Tile1_Data',
            'B06': 'IMG_DATA_Band_B06_20m_Tile1_Data',
            'B07': 'IMG_DATA_Band_B07_20m_Tile1_Data',
            'B08': 'IMG_DATA_Band_B08_10m_Tile1_Data',
            'B8A': 'IMG_DATA_Band_B8A_20m_Tile1_Data',
            'B11': 'IMG_DATA_Band_B11_20m_Tile1_Data',
            'B12': 'IMG_DATA_Band_B12_20m_Tile1_Data',
            'SCL': 'IMG_DATA_Band_SCL_20m_Tile1_Data',
            'TCI': 'IMG_DATA_Band_TCI_20m_Tile1_Data'
        }

    # product metadata
    product_info = api.get_product_odata(scene)
    title = product_info['title']
    availability = product_info['Online']
    if availability == False:
        title = product_info['title']
        logging.info(
            title + ' is offline. You must activate the image first'+'\n')
        sys.exit(0)
    # get manifest file
    manifest = getManifest(product_info, api_session)
    if product_type == 'S2MSI1C':
        for band in bands:
            # url generation for single band
            band_id = band_dict[band]

            file_info = [
                file_info for file_info in manifest if file_info['id'] == band_id]
            if len(file_info) == 1:
                file_info = file_info[0]
                file_info['url'] = '/'.join(product_info['url'].split('/')
                                            [:-1]) + "/Nodes('{}.SAFE')/".format(title)
                file_info['url'] += '/'.join(["Nodes('{}')".format(token)
                                              for token in file_info['href'].split('/')[1:]]) + '/$value'
                # print('Downloading '+file_info['url']+' ...')
                logging.info('Downloading ' + title + ' band: '+band)
                # download single band
                res = downloadS2Band(
                    product_type, file_info['url'], band, api_session, title, output_folder_s2)
                if res == False:
                    logging.error('Error Downloading ' +
                                  title + 'band: '+band_id)
                    sys.exit(1)

    if product_type == 'S2MSI2A':
        for band in bands:
            # url generation for single band
            band_id = band_dict[band]
            file_info = [
                file_info for file_info in manifest if file_info['id'] == band_id]
            if len(file_info) == 1:
                file_info = file_info[0]
                title = product_info['title']
                url = product_info['url'].split(
                    '/$')[0]+"/Nodes('"+title+".SAFE')/Nodes('GRANULE')/Nodes('"
                href = file_info['href'].split("/GRANULE/")[1]
                id_granule = href.split('/')[0]
                id_file = href.split('/')[-1]
                if band == 'B08':
                    url = url+id_granule + \
                        "')/Nodes('IMG_DATA')/Nodes('R10m')/Nodes('" + \
                        id_file+"')/$value"
                else:
                    url = url+id_granule + \
                        "')/Nodes('IMG_DATA')/Nodes('R20m')/Nodes('" + \
                        id_file+"')/$value"
                # print('Downloading '+url+' ...')
                logging.info('Downloading ' + title + 'band: '+band_id)
                # download single band
                res = downloadS2Band(product_type, url, band,
                                     api_session, title, output_folder_s2)
                if res == False:
                    logging.error('Error Downloading ' +
                                  title + 'band: '+band_id)
                    sys.exit(1)

    scene_path = str(pathlib.Path(os.path.join(output_folder_s2, title)))
    return scene_path, title


def searchS2(footprint, product_type, start_date, end_date, cloud_cover, s2grid_path, central_date, quicklook, user_esa, pass_esa, output_res, run_mode):
    '''
    Description:
    ------------
    Searches for Sentinel 2 scenes using ESA api through SentinelSat module.
    Scenes with the same date will be removed (case of areas of interest that
    overlap the same area in different coordinate reference systems)

    WARNING: currently, the Copernicus server only provide access to the online products during one month.
             Products older than one month are identified as 'offline' products and have to be ordered firstly
             to the ESA to become 'offline' to 'online'. This can takes a long time (more than one hour in some cases).
             For more information google 'Long Term Access for sentinel products'.

    Arguments:
    ------------
    - footprint (string): region of interest in wkt format.
    Example: GEOMETRYCOLLECTION(POLYGON((-0.4532 39.3975,-0.2856 39.3975,-0.2856 39.5295,-0.4532 39.5295,-0.4532 39.3975)
    - product_type (string): type of product for S2 images. Can be S2MSI1C or S2MSI2A
    - start_date (string): start date in format YYYYMMDD
    - end_date (string): end date in format YYYYMMDD
    - cloud_cover (int): cloud coverage. Number between 0 and 100.
    - s2grid_path (string): path to the shapefile of Sentinel-2 grid
    - central_date (string): central date in format YYYYMMDD
    - quicklook (string): quicklook activation (0: deactivated, 1: activated)
    - user_esa (string): user credential for the ESA server
    - pass_esa (string): password credential for the ESA server
    - output_res (string): tipe of output for the searching results (CONSOLE:0, TXT FILE:1, JSON FILE:2)
    - run_mode (string): run mode (os, dp, od, op)

    Returns:
    ------------
    - filtered_scenes (list of dictionaries): list of metadata for every found scene.
    - filtered_titles (list of strings): list of ids for every found scene.

    '''

    # authentication
    esa_api = esaAuthentication(user_esa, pass_esa)

    # arguments for query
    try:
        scenes = esa_api.query(footprint,
                               date=(start_date, end_date),
                               platformname='Sentinel-2',
                               producttype=product_type,
                               cloudcoverpercentage=(0, cloud_cover))

        filtered_scenes = []
        filtered_titles = []
        filtered_clouds = []
        filtered_quicklooks = []
        filtered_availability = []
        txt_string = ''
        if len(scenes) > 0:
            for scene in scenes:
                # get scene metadata
                product_info = esa_api.get_product_odata(scene, full=True)
                # print (product_info)
                scene_title = product_info['title']
                scene_cloud = str(
                    round(float(product_info['Cloud cover percentage']), 2))
                scene_ql = product_info['quicklook_url']
                # avoiding 'offline' products and removing scenes with the same date
                if run_mode != 'op':
                    if product_info['Online'] == True:
                        if output_res == '0':
                            print(
                                f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online')
                        else:
                            txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online'+'\n'
                        filtered_scenes.append(scene)
                        filtered_titles.append(scene_title)
                        filtered_clouds.append(scene_cloud)
                        filtered_quicklooks.append(scene_ql)
                        filtered_availability.append('Online')
                    else:
                        if output_res == '0':
                            print(
                                f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline')
                        else:
                            txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline'+'\n'
                        filtered_scenes.append(scene)
                        filtered_titles.append(scene_title)
                        filtered_clouds.append(scene_cloud)
                        filtered_quicklooks.append(scene_ql)
                        filtered_availability.append('Offline')
                # run mode p ==> all images, online and offline are included (NEW*****)
                else:
                    if product_info['Online'] == True:
                        filtered_availability.append('Offline')
                        if output_res == '0':
                            print(
                                f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online')
                        else:
                            txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online'+'\n'
                    else:
                        filtered_availability.append('Offline')
                        if output_res == '0':
                            print(
                                f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline')
                        else:
                            txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline'+'\n'
                    filtered_scenes.append(scene)
                    filtered_titles.append(scene_title)
                    filtered_clouds.append(scene_cloud)
                    filtered_quicklooks.append(scene_ql)

            print('\n')

            return filtered_scenes, filtered_titles, filtered_quicklooks, filtered_clouds, txt_string, filtered_availability
        else:
            logging.warning('There are no products S2 to download.')
            sys.exit(1)
            # return [], [], [], [] , txt_string

    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def searchS2ByTiles(tiles, product_type, start_date, end_date, cloud_cover, quicklook, user_esa, pass_esa, output_res):
    '''
    Description:
    ------------
    Searches for Sentinel 2 scenes using ESA api through SentinelSat module.

    WARNING: currently, the Copernicus server only provide access to the online products during one month.
             Products older than one month are identified as 'offline' products and have to be ordered firstly
             to the ESA to become 'offline' to 'online'. This can takes a long time (more than one hour in some cases).
             For more information google 'Long Term Access for sentinel products'.

    Arguments:
    ------------
    - tiles (list): list of Sentinel 2 tiles ('30TYJ','30TFC')
    - product_type (string): type of product for S2 images. Can be S2MSI1C or S2MSI2A
    - start_date (string): start date in format YYYYMMDD
    - end_date (string): end date in format YYYYMMDD
    - cloud_cover (int): cloud coverage. Number between 0 and 100.
    - quicklook (string): quicklook activation (0: deactivated, 1: activated)
    - user_esa (string): user credential for the ESA server
    - pass_esa (string): password credential for the ESA server
    - output_res (string): tipe of output for the searching results (CONSOLE:0, TXT FILE:1, JSON FILE:2)

    Returns:
    ------------
    - filtered_scenes (list of dictionaries): list of metadata for every found scene.
    - filtered_titles (list of strings): list of ids for every found scene.

    '''
    try:
        # authentication
        esa_api = esaAuthentication(user_esa, pass_esa)
        # print(esa_api)

        # arguments for query
        query_kwargs = {
            'platformname': 'Sentinel-2',
            'producttype': product_type,
            'date': (start_date, end_date),
            'cloudcoverpercentage': (0, cloud_cover)}

        scenes = OrderedDict()
        for tile in tiles:
            kw = query_kwargs.copy()
            # kw['tileid'] = tile
            kw['filename'] = '*_T{}_*'.format(tile)
            pp = esa_api.query(**kw)
            scenes.update(pp)

        filtered_scenes = []
        filtered_titles = []
        filtered_quicklooks = []
        filtered_clouds = []
        filtered_availability = []
        txt_string = ''

        if len(scenes) > 0:
            ordered_scenes = []
            img_list = {}
            for scene in scenes:
                product_info = esa_api.get_product_odata(scene, full=True)
                scene_title = product_info['title']
                scene_date = scene_title.split('_')[2][0:8]
                img_list[scene_date] = scene
            scenes_sorted = dict(sorted(img_list.items(), reverse=True))
            for scene in scenes_sorted.values():
                ordered_scenes.append(scene)

            for scene in ordered_scenes:
                # get scene metadata
                product_info = esa_api.get_product_odata(scene, full=True)
                scene_title = product_info['title']
                scene_cloud = str(
                    round(float(product_info['Cloud cover percentage']), 2))
                scene_ql = product_info['quicklook_url']
                # avoiding 'offline' products and removing scenes with the same date
                if product_info['Online'] == True:
                    if output_res == '0':
                        print(
                            f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online')
                    else:
                        txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: online'+'\n'
                    filtered_scenes.append(scene)
                    filtered_titles.append(scene_title)
                    filtered_clouds.append(scene_cloud)
                    filtered_quicklooks.append(scene_ql)
                    filtered_availability.append('Online')
                else:
                    if output_res == '0':
                        print(
                            f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline')
                    else:
                        txt_string += f'Scene: {scene_title} Cloud coverage: {scene_cloud} availability: offline'+'\n'
                    # if quicklook == '1':
                    #    webbrowser.open(scene_ql)
                    filtered_scenes.append(scene)
                    filtered_titles.append(scene_title)
                    filtered_clouds.append(scene_cloud)
                    filtered_quicklooks.append(scene_ql)
                    filtered_availability.append('Offline')

            print('\n')
            return filtered_scenes, filtered_titles, filtered_quicklooks, filtered_clouds, txt_string, filtered_availability
        else:
            logging.warning('There are no products S2 to download.')
            sys.exit(1)
            # return [], [], [], []

    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def filterS2TilesByOverlap(titles, shp_path):
    '''
    Description:
    ------------
    Prevents Sentinel 2 scenes overlapping. Same times two scenes covers similar area
    in the same date. This happens when we have two UTM zones over the same area

    Arguments:
    ------------
    - titles (string): list of sentinel 2 scene names
    - shp_path (string): path to shapefile containing sentinel 2 scene footprints

    Returns:
    ------------
    - filter_tiles (list): list of filtered tiles

    '''

    tiles = []
    date_tiles = []
    for title in titles:
        tile = title.split('_')[5][1:]
        date_tile = title.split('_')[2].split('T')[0]
        tiles.append(tile)
        date_tiles.append(date_tile)

    # opens shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    input_ds = driver.Open(shp_path)
    input_layer = input_ds.GetLayer()

    scene_ids_u = list(set(tiles))
    footprints = {}
    # filter shapefile by scene tile name
    for scene_id in scene_ids_u:
        input_layer.SetAttributeFilter("Name = '"+scene_id+"'")
        for feat in input_layer:
            feat_footprint = feat.GetGeometryRef()
            footprints[scene_id] = feat_footprint.Clone()

    # check if one single scene intersects with other scene more than 50%
    final_ids = []
    for i in range(0, len(scene_ids_u)):
        for j in range(i+1, len(scene_ids_u)):
            geom_inter = footprints[scene_ids_u[i]].Intersection(
                footprints[scene_ids_u[j]])
            if not geom_inter is None:
                area = geom_inter.GetArea()
                if area > 0.5:
                    final_ids.append(scene_ids_u[i])

    # filter final scenes taht not overlap
    unique_ids = list(set(final_ids))
    filter_tiles = []
    for title in titles:
        if not title.split('_')[5][1:] in unique_ids:
            filter_tiles.append(title)

    return filter_tiles


def getNearestDates(filter_tiles, central_date):
    '''
    Description:
    ------------
    Returns the closest tile to a central date in both pre and post date cases.

    Arguments:
    ------------
    - filter_tiles (list): list of scene names
    - central_date (string): central date

    Returns:
    ------------
    - final (list): list of two scenes (pre and post date)

    '''
    dates = []
    for title in filter_tiles:
        tile_date = title.split('_')[2].split('T')[0]
        dates.append(tile_date)
    dates.sort()
    max_post = max_pre = 1000000

    for date in dates:
        if date >= central_date:
            offset_post = (datetime.datetime.strptime(
                date, '%Y%m%d')-datetime.datetime.strptime(central_date, '%Y%m%d')).days
            if offset_post < max_post:
                min_post = date
                max_post = offset_post
        else:
            offset_pre = (datetime.datetime.strptime(
                central_date, '%Y%m%d')-datetime.datetime.strptime(date, '%Y%m%d')).days
            if offset_pre < max_pre:
                min_pre = date
                max_pre = offset_pre

    final = []
    for title in filter_tiles:
        date_tile = title.split('_')[2].split('T')[0]
        if date_tile in [min_pre, min_post]:
            final.append(title)
    return final


def searchLandsat8(footprint, product_type, start_date, end_date, max_cloud_cover, user_usgs, pass_usgs):
    '''
    Description:
    ------------
    Searches for Landsat 8 scenes using USGS api through landsatxplore module.
    Scenes with the same date will be removed (case of areas of interest that
    overlap the same area in different coordinate reference systems)

    Arguments:
    ------------
    - footprint (tuple of 4 strings): region of interest in format (xmin, ymin, xmax, ymax).
    - product_type (string): type of product for L8 images (landsat_8_c1).
    - start_date (string): start date in format YYYYMMDD
    - end_date (string): end date in format YYYYMMDD
    - max_cloud_cover (int): cloud coverage. Number between 0 and 100.
    - user_usgs (string): user credential for the USGS server
    - pass_usgs (string): password credential for the USGS server

    Returns:
    ------------
    - filtered_scenes (list of dictionaries): list of metadata for every found scene.

    '''

    # authentication
    ee_api_search, ee_api_download = usgsAuthentication(user_usgs, pass_usgs)

    try:
        scenes = ee_api_search.search(
            dataset=product_type,
            bbox=footprint,
            start_date=start_date,
            end_date=end_date,
            max_cloud_cover=max_cloud_cover
        )

        ee_api_search.logout()

        if(len(scenes) > 0):
            return scenes
        else:
            return []
    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def searchLandsat8byScene(scene_id, product_type, start_date, end_date, max_cloud_cover, l8grid_path, user_usgs, pass_usgs):
    '''
    Description:
    ------------
    Searches for Landsat 8 scenes using USGS api through landsatxplore module.
    Scenes with the same date will be removed (case of areas of interest that
    overlap the same area in different coordinate reference systems)

    Arguments:
    ------------
    - scene (string): scene code for landsat 8 (199032).
    - product_type (string): type of product for L8 images (landsat_8_c1).
    - start_date (string): start date in format YYYYMMDD
    - end_date (string): end date in format YYYYMMDD
    - max_cloud_cover (int): cloud coverage. Number between 0 and 100.
    - l8grid_path (string): path to the shapefile of Landsat grid
    - user_usgs (string): user credential for the USGS server
    - pass_usgs (string): password credential for the USGS server

    Returns:
    ------------
    - filtered_scenes (list of dictionaries): list of metadata for every found scene.

    '''

    # authentication
    ee_api_search, ee_api_download = usgsAuthentication(user_usgs, pass_usgs)

    try:
        long, lat = getLatLongFromPathRow(l8grid_path, scene_id)
        scenes = ee_api_search.search(
            dataset=product_type,
            longitude=long,
            latitude=lat,
            start_date=start_date,
            end_date=end_date,
            max_cloud_cover=max_cloud_cover
        )

        ee_api_search.logout()

        img_list = []
        for scene in scenes:
            title = scene['display_id']
            # Avoid L1GT products. Just only L1PT
            if not 'L1GT' in title:
                if scene_id in title:
                    img_list.append(scene)

        return img_list

    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def getLatLongFromPathRow(shp_path, scene_id):
    '''
    Description:
    ------------
    Returns long, lat coordinates from path and row code for landsat 8 scenes

    Arguments:
    ------------
    - shp_pat (string): path to the shapefile containing Landsat 8 footprints.
    - scene_id (string): scene code (199032).

    Returns:
    ------------
    - long, lat (float): longitud and latitud of the scene footprint centroid

    '''
    path_id = int(scene_id[0:3])
    row_id = int(scene_id[3:])

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(shp_path)
    layer = ds.GetLayer()
    layer.SetAttributeFilter("PATH = "+str(path_id)+' and ROW = '+str(row_id))
    for feat in layer:
        geom = feat.GetGeometryRef()
        long = geom.Centroid().GetX()
        lat = geom.Centroid().GetY()
    ds = None
    return long, lat


def downloadLandsat8_2(scene, base_folder, user_usgs, pass_usgs):
    '''
    Description:
    ------------
    Downloads one single L8 scene obtained from the searchLandsat8() function.

    Arguments:
    ------------
    - scene (dictionary): metadata for a single L8 scene
    - base_folder (string): folder to download L8 scene
    - user_usgs (string): user credential for the USGS server
    - pass_usgs (string): password credential for the USGS server

    Returns:
    ------------
    - scene_path (string): path to the folder where the L8 scene has been downloaded

    '''

    # USGS authentication
    ee_api_search, ee_api_download = usgsAuthentication(user_usgs, pass_usgs)

    # get scene id
    scene_id = scene['display_id']

    # output folder creation
    if 'LC08' in scene_id:
        createFolderCheck(str(pathlib.Path(os.path.join(base_folder, 'l8'))))
        output_folder_l = str(pathlib.Path(os.path.join(base_folder, 'l8')))
    if 'LC09' in scene_id:
        createFolderCheck(str(pathlib.Path(os.path.join(base_folder, 'l9'))))
        output_folder_l = str(pathlib.Path(os.path.join(base_folder, 'l9')))
    try:
        # print('Downloading... '+scene_id)
        logging.info('Downloading... '+scene_id)
        # download L8 .tar file to the output folder
        ee_api_download.download(scene_id, output_folder_l)
        createFolder(
            str(pathlib.Path(os.path.join(output_folder_l, scene_id))))
        output_folder_date = str(pathlib.Path(
            os.path.join(output_folder_l, scene_id)))
        # uncompress the file
        # print('Unzipping... '+scene_id)
        logging.info('Unzipping... '+scene_id)
        if os.path.exists(str(pathlib.Path(os.path.join(output_folder_l, scene_id+'.tar')))):
            # be careful -> zip file .zip, .tar, .tar.gz
            untarFile(str(pathlib.Path(os.path.join(output_folder_l,
                                                    scene_id+'.tar'))), output_folder_date)
        if os.path.exists(str(pathlib.Path(os.path.join(output_folder_l, scene_id+'.tar.gz')))):
            untarFile(str(pathlib.Path(os.path.join(output_folder_l,
                                                    scene_id+'.tar.gz'))), output_folder_date)
        if os.path.exists(str(pathlib.Path(os.path.join(output_folder_l, scene_id+'.zip')))):
            unzipFile(str(pathlib.Path(os.path.join(output_folder_l,
                                                    scene_id+'.zip'))), output_folder_date)

        ee_api_search.logout()

        scene_path = str(pathlib.Path(os.path.join(output_folder_l, scene_id)))
        return scene_path, scene_id

    except Exception as e:
        logging.error('Exception: %s', e)
        sys.exit(1)


def getManifest(product_info, api_session):
    '''
    Description:
    ------------
    Returns manifest file (xml format) for a S2 scene

    Arguments:
    ------------
    - product_info (dictionary): metadata of the scene
    - api_sesion (object): requests.Session() object. Allows url requests with credential authentication

    Returns:
    ------------
    - _xml (dictionary): manifest file content in dictionary structure

    '''

    # manifest url setting
    _title = product_info['title']
    _url = product_info['url']
    manifest = _url.split('/$')[0]+"/Nodes('" + \
        _title+".SAFE')/Nodes('manifest.safe')/$value"

    # get manifest content
    while True:
        response = api_session.get(manifest)
        if response.status_code == 200:
            _xml = parse_manifest_xml(response.content)
            return _xml


def parse_manifest_xml(xml):
    '''
    Description:
    ------------
    Converts manifest xml content to a dictionary structure

    Arguments:
    ------------
    - xml (string): manifest content in xml format

    Returns:
    ------------
    - outputs (list): list with manifest content in a dictionary structure)

    '''

    outputs = []
    root = ET.fromstring(xml)
    for item in root.findall("./dataObjectSection/dataObject"):
        output = {
            'id': item.get('ID'),
            'mimetype': item.find('./byteStream').get('mimeType'),
            'size': int(item.find('./byteStream').get('size')),
            'href': item.find('./byteStream/fileLocation').get('href'),
            'md5sum': item.find('./byteStream/checksum').text
        }
        outputs.append(output)

    return outputs


def downloadS2Band(product_type, url, band, api_session, title, base_folder):
    '''
    Description:
    ------------
    Donwloads single S2 band from url request.


    Arguments:
    ------------
    - product_type (string): type of product for S2 image (2MSI1C or S2MSI2A).
    - url (string): url to download the S2 scene
    - band (string): band name to create the final output path
    - api_session (object): requests.Session object to make url requests with credential authentication
    - tilte (string): scene id
    - base_folder (string): path to download the scene band

    Returns:
    ------------
    True or False depending on the request success

    '''

    # create output folder if it is needed
    createFolderCheck(str(pathlib.Path(os.path.join(base_folder, title))))
    output_folder = str(pathlib.Path(os.path.join(base_folder, title)))

    # name of the final jp2 band file depending on the product type
    if product_type == 'S2MSI2A':
        if band == 'B08':
            band_name = title+'_'+band+'_10m'+'.jp2'
        else:
            band_name = title+'_'+band+'_20m'+'.jp2'
        img_path = os.path.join(output_folder, band_name)

    if product_type == 'S2MSI1C':
        if '.jp2' in url:
            band_name = title+'_'+band+'.jp2'  # ojo
        else:
            band_name = title+'_'+band+'.gml'
        img_path = os.path.join(output_folder, band_name)

    # url modification to adapt the request to products with different folder structrure
    # ESA uses different folder structure for the same product for older dates.
    parts = url.split('Nodes')
    url2 = parts[0]+'Nodes'+parts[1]+'Nodes'+parts[-1]
    url_base = url

    # img_name = img_path.split('\\')[-1]
    img_name = pathlib.Path(img_path).parts[-1]

    # to prevent wrong answers fron the server, 10 request are performed
    # the more commun issues came from situations for bussy server or for bad url specification
    tries = 0
    while True:
        r = api_session.get(url_base, stream=True)
        # print('bucle 1',r.status_code)
        if r.status_code == 200:
            file_size = int(r.headers.get("Content-Length", 0))
            progress = tqdm(r.iter_content(
                1024), f"Downloading {img_name}", total=file_size, unit="B", unit_scale=True, unit_divisor=1024)
            with open(img_path, "wb") as f:
                for data in progress.iterable:
                    f.write(data)
                    progress.update(len(data))

            return True
        tries += 1
        if tries == 20:
            break

    # url modification to adapt the request to products with different folder structrure
    # ESA uses different folder structure for the same product for older dates.
    url_base = url2

    # to prevent wrong answers fron the server, 10 request are performed
    # the more commun issues came from situations for bussy server or for bad url specification
    tries = 0
    while True:
        r = api_session.get(url_base, stream=True)
        if r.status_code == 200:
            file_size = int(r.headers.get("Content-Length", 0))
            progress = tqdm(r.iter_content(
                1024), f"Downloading {img_name}", total=file_size, unit="B", unit_scale=True, unit_divisor=1024)
            with open(img_path, "wb") as f:
                for data in progress.iterable:
                    f.write(data)
                    progress.update(len(data))

            return True
        tries += 1
        if tries == 20:
            break

    return False  # something wrong happend (more than 10 tries).


def startSearchForLandsat8(run_parameters):
    '''
    Description:
    ------------
    Adapts the format of start_date and end_date for searching L8 scenes.
    Get L8 scene footprint in suitable format for searching L8 scenes.
    Runs the function searchLandsat8 to obtain scenes that matches the initial conditions (footprint, date interval and cloud coverage)

    Arguments:
    ------------
    - run_parameters (list of strings): list of parameters.

    Returns:
    ------------
    l8_scenes (list of dictionaries): list of metadata for every found scene.

    '''
    # get run parameters
    footprint = run_parameters['footprint']
    start_date = run_parameters['start_date']
    central_date = run_parameters['central_date']
    end_date = run_parameters['end_date']
    l8_product = run_parameters['l8_product']
    max_cloud_cover = run_parameters['max_cloud_cover']
    scene_landsat_list = run_parameters['scene_landsat_list']
    l8grid_path = run_parameters['l8grid_path']
    quicklook = run_parameters['quicklook']
    user_usgs = run_parameters['user_usgs']
    pass_usgs = run_parameters['pass_usgs']
    output_res = run_parameters['output_res']
    output_search_folder = run_parameters['output_search_folder']
    run_mode = run_parameters['run_mode']

    # extra settings for donwloading run mode
    if run_mode != 'os':
        output_res = '0'
        quicklook = '0'

    # date adaptation for searching L8 scenes
    start_date_l8 = start_date[0:4]+'-'+start_date[4:6]+'-'+start_date[6:8]
    end_date_l8 = end_date[0:4]+'-'+end_date[4:6]+'-'+end_date[6:8]

    txt_string = ''
    if footprint != 'NONE':  # searchinng by footprint
        # footprint adaptation for searching L8 scenes
        bbox = getLandsatBbox(footprint)
        # print('Searching for L8 images...'+'\n')
        logging.info('Searching for Landsat images...'+'\n')
        l8_scenes = searchLandsat8(
            bbox, l8_product, start_date_l8, end_date_l8, max_cloud_cover, user_usgs, pass_usgs)
        l8_scenes_filtered = {}

        if len(l8_scenes) == 0:
            # print('There are no L8 images to download'+'\n')
            logging.warning('There are no Landsat images to download'+'\n')
            # sys.exit(1)
            return []
        else:
            flag = 0
            ni = 0
            image_url_list = []
            image_title_list = []
            for i in range(0, len(l8_scenes)):
                scene = l8_scenes[i]
                if not 'L1GT' in scene['display_id']:
                    l8_scenes_filtered[str(ni)] = scene
                    scene_date = scene['display_id'].split('_')[3]
                    date_difference = datetime.datetime.strptime(
                        scene_date, "%Y%m%d") - datetime.datetime.strptime(central_date, "%Y%m%d")
                    cloud_cover = str(
                        round(float(scene['scene_cloud_cover']), 2))
                    if central_date < scene_date:
                        if output_res == '0':
                            print(
                                f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                        else:
                            txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'
                    else:
                        if flag == 0:
                            if output_res == '0':
                                print('[*******] Central date: '+central_date)
                                print(
                                    f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                            else:
                                txt_string += '[*******] Central date: ' + \
                                    central_date+'\n'
                                txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'
                            flag = 1
                        else:
                            if output_res == '0':
                                print(
                                    f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                            else:
                                txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'
                    image_url_list.append(getLandsatQuickLookUrl(
                        scene['display_id'], l8_product))
                    image_title_list.append(scene['display_id']+' CC:' +
                                            str(scene['scene_cloud_cover']) + '%  Days: '+str(date_difference.days))
                    ni += 1
            if flag == 0:
                if output_res == '0':
                    print('[*******] Central date: '+central_date)
                else:
                    txt_string += '[*******] Central date: '+central_date+'\n'
            print('\n')
            if output_res == '1':
                txt_filename = str(pathlib.Path(os.path.join(
                    output_search_folder, 'search_result.txt')))
                with open(txt_filename, 'w') as txt_file:
                    txt_file.write(txt_string)
                logging.info(
                    'Writing search results in file '+output_search_folder+'/search_result.txt')
            if output_res == '2':
                results = exportLandsatResultsToJson(txt_string)
                json_path = str(pathlib.Path(os.path.join(
                    output_search_folder, 'search_result.json')))
                with open(str(json_path), 'w') as outfile:
                    json.dump(results, outfile)
                logging.info('Writing search results in file ' +
                             output_search_folder+'/search_result.json')

            html_file = writeHtml(image_url_list, image_title_list,
                                  output_search_folder, l8_product)

            if quicklook == '1':  # quicklook file is created but it only is opened when quicklook variable = 1
                webbrowser.open(html_file)

            return l8_scenes_filtered
    else:
        if scene_landsat_list != 'NONE':
            # print('Searching for L8 images...'+'\n')
            logging.info('Searching for Landsat images...'+'\n')
            l8_scenes = []
            img_list = {}
            ordered_scenes = []
            l8_scenes_filtered = {}
            for scene_id in scene_landsat_list:
                partial_search = searchLandsat8byScene(
                    scene_id, l8_product, start_date_l8, end_date_l8, max_cloud_cover, l8grid_path, user_usgs, pass_usgs)
                l8_scenes = l8_scenes+partial_search
            if len(l8_scenes) == 0:
                # print('There are no L8 images to download'+'\n')
                logging.warning('There are no Landsat images to download'+'\n')
                # exit(1)
                return []
            else:
                txt_string = ''
                for scene in l8_scenes:
                    id_scene = scene['display_id']
                    img_date = id_scene.split('_')[3]
                    if not img_date in img_list:
                        img_list[img_date] = [scene]
                    else:
                        img_list[img_date] += [scene]

                scenes_sorted = dict(sorted(img_list.items(), reverse=True))
                for scenes in scenes_sorted.values():
                    for scene in scenes:
                        ordered_scenes.append(scene)
                flag = 0
                ni = 0
                image_url_list = []
                image_title_list = []
                for i in range(0, len(ordered_scenes)):
                    scene = ordered_scenes[i]
                    l8_scenes_filtered[str(ni)] = scene
                    # filtered_l8_scenes.append(scene)
                    scene_date = scene['display_id'].split('_')[3]
                    date_difference = datetime.datetime.strptime(
                        scene_date, "%Y%m%d") - datetime.datetime.strptime(central_date, "%Y%m%d")
                    cloud_cover = str(
                        round(float(scene['scene_cloud_cover']), 2))
                    if central_date < scene_date:
                        if output_res == '0':
                            print(
                                f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                        else:
                            txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'
                    else:
                        if flag == 0:
                            if output_res == '0':
                                print('[*******] Central date: '+central_date)
                                print(
                                    f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                            else:
                                txt_string += '[*******] Central date: ' + \
                                    central_date+'\n'
                                txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'
                            flag = 1
                        else:
                            if output_res == '0':
                                print(
                                    f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days.")
                            else:
                                txt_string += f"[{ni}] Scene: {scene['display_id']} Cloud cover: {cloud_cover}%  {date_difference.days} days."+'\n'

                    image_url_list.append(getLandsatQuickLookUrl(
                        scene['display_id'], l8_product))
                    image_title_list.append(scene['display_id']+' CC:' +
                                            str(scene['scene_cloud_cover']) + '%  Days: '+str(date_difference.days))
                    ni += 1
                if flag == 0:
                    if output_res == '0':
                        print('[*******] Central date: '+central_date)
                    if output_res == '1':
                        txt_string += '[*******] Central date: ' + \
                            central_date+'\n'
                print('\n')
                if output_res == '1':
                    txt_filename = str(pathlib.Path(os.path.join(
                        os.getcwd(), output_search_folder, 'search_result.txt')))
                    with open(txt_filename, 'w') as txt_file:
                        txt_file.write(txt_string)
                    logging.info(
                        'Writing search results in file '+output_search_folder+'/search_result.txt')
                if output_res == '2':
                    results = exportLandsatResultsToJson(txt_string)
                    json_path = str(pathlib.Path(os.path.join(
                        output_search_folder, 'search_result.json')))
                    with open(str(json_path), 'w') as outfile:
                        json.dump(results, outfile)
                    logging.info(
                        'Writing search results in file '+output_search_folder+'/search_result.json')

                html_file = writeHtml(image_url_list, image_title_list,
                                      output_search_folder, l8_product)
                if quicklook == '1':  # quicklook file is created but it only is opened when quicklook variable = 1
                    webbrowser.open(html_file)

                return l8_scenes_filtered  # ordered_scenes
        else:
            # print('If there is not footprint defined, you must define Landsat scenes list')
            logging.warning(
                'If there is not footprint defined, you must define Landsat scenes list')
            sys.exit(1)


def getPrePostStormDateLandsat(central_date, scenes):
    '''
    Not used function ********************
    '''
    scene_post = ''
    scene_pre = ''
    for i in range(0, len(scenes)-1):
        d_post = scenes[i]['display_id'].split('_')[3]
        if d_post > central_date:
            scene_post = scenes[i]

    for i in range(0, len(scenes)-1):
        d_pre = scenes[i]['display_id'].split('_')[3]
        if d_pre <= central_date:
            scene_pre = scenes[i]
            break
    if scene_post == '' or scene_pre == '':
        print('Pre or post-storm candidate has not been found.')
        sys.exit()
    else:
        return [scene_post, scene_pre]


def getPrePostStormDateSentinel(central_date, s2_scenes, s2_titles, quicklooks):
    '''
    Not used function ********************
    '''
    scene_post = ''
    scene_pre = ''
    for i in range(0, len(s2_titles)-1):
        d_post = s2_titles[i].split('_')[2][0:8]
        if d_post > central_date:
            scene_post = s2_scenes[i]
            quicklook_post = quicklooks[i]
            title_post = s2_titles[i]

    for i in range(0, len(s2_titles)-1):
        d_pre = s2_titles[i].split('_')[2][0:8]
        if d_pre <= central_date:
            scene_pre = s2_scenes[i]
            quicklook_pre = quicklooks[i]
            title_pre = s2_titles[i]
            break
    if scene_post == '' or scene_pre == '':

        print('Pre or post-storm candidate has not been found.')
        sys.exit()
    else:
        return [scene_post, scene_pre], [quicklook_post, quicklook_pre], [title_post, title_pre]


def startSearchForSentinel2(run_parameters, ni):
    '''
    Description:
    ------------
    Runs the function searchS2 to obtain scenes that matches the initial conditions (footprint, date interval and cloud coverage)

    Arguments:
    ------------
    - run_parameters (list of strings): list of parameters.
    - ni (int): stores the whole number of found images (combined with the found Landsat images if necessary).
.
    Returns:
    ------------
    s2_scenes (list of dictionaries): list of metadata for every found scene.

    '''

    # get run parameters
    footprint = run_parameters['footprint']
    start_date = run_parameters['start_date']
    central_date = run_parameters['central_date']
    end_date = run_parameters['end_date']
    s2_product = run_parameters['s2_product']
    max_cloud_cover = run_parameters['max_cloud_cover']
    scene_sentinel_list = run_parameters['scene_sentinel_list']
    s2grid_path = run_parameters['s2grid_path']
    quicklook = run_parameters['quicklook']
    s2_quicklook_server = run_parameters['s2_quicklook_server']
    user_esa = run_parameters['user_esa']
    pass_esa = run_parameters['pass_esa']
    output_res = run_parameters['output_res']
    output_search_folder = run_parameters['output_search_folder']
    run_mode = run_parameters['run_mode']

    # extra settings for donwloading run mode
    if run_mode != 'os':
        output_res = '0'
        quicklook = '0'

    s2_scenes_filtered = {}
    if footprint != 'NONE':
        # print('Searching for S2 images...'+'\n')
        logging.info('Searching for S2 images...'+'\n')
        s2_scenes, s2_titles, s2_quicklooks, s2_clouds, txt_string, s2_availability = searchS2(
            footprint, s2_product, start_date, end_date, max_cloud_cover, s2grid_path, central_date,
            quicklook, user_esa, pass_esa, output_res, run_mode)
        if len(s2_scenes) == 0:
            # print('There are not products S2 to download.'+'\n')
            if txt_string == '':
                logging.warning('There are no products S2 to download.'+'\n')
                sys.exit(1)
            else:
                logging.warning(
                    'There are no products S2 online to download.'+'\n')
                if output_res == '1':
                    txt_filename = str(pathlib.Path(os.path.join(
                        os.getcwd(), output_search_folder, 'search_result.txt')))
                    with open(txt_filename, 'w') as txt_file:
                        txt_file.write(txt_string)
                    logging.info(
                        'Writing search results in file '+output_search_folder+'/search_result.txt')
                if output_res == '2':
                    results = exportSentinelResultsToJson(txt_string)
                    json_path = str(pathlib.Path(os.path.join(
                        output_search_folder, 'search_result.json')))
                    with open(str(json_path), 'w') as outfile:
                        json.dump(results, outfile)
                    logging.info(
                        'Writing search results in file '+output_search_folder+'/search_result.json')
                sys.exit(1)
        else:
            txt_string += '\n'
            flag = 0
            #ni = 0
            image_url_list = []
            image_title_list = []
            for i in range(0, len(s2_scenes)):
                s2_scenes_filtered[str(ni)] = s2_scenes[i]
                title = s2_titles[i]
                scene_date = title.split('_')[2][0:8]
                date_difference = datetime.datetime.strptime(
                    scene_date, "%Y%m%d") - datetime.datetime.strptime(central_date, "%Y%m%d")
                if central_date < scene_date:
                    if output_res == '0':
                        print(
                            f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                    else:
                        txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                else:
                    if flag == 0:
                        if output_res == '0':
                            print('[*******] Central date:'+central_date)
                            print(
                                f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                        else:
                            txt_string += '[*******] Central date:' + \
                                central_date+'\n'
                            txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                        flag = 1
                    else:
                        if output_res == '0':
                            print(
                                f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                        else:
                            txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                image_url_list.append(s2_quicklooks[i])
                image_title_list.append(
                    title+' CC:' + str(s2_clouds[i]) + '  Days: '+str(date_difference.days) + ' '+s2_availability[i])
                ni += 1

            if flag == 0:
                if output_res == '0':
                    print('[*******] Central date: '+central_date)
                else:
                    txt_string += '[*******] Central date: '+central_date+'\n'

            print('\n')
            if output_res == '1':
                txt_filename = str(pathlib.Path(os.path.join(
                    os.getcwd(), output_search_folder, 'search_result.txt')))
                with open(txt_filename, 'w') as txt_file:
                    txt_file.write(txt_string)
                logging.info(
                    'Writing search results in file '+output_search_folder+'/search_result.txt')
            if output_res == '2':
                results = exportSentinelResultsToJson(txt_string)
                json_path = str(pathlib.Path(os.path.join(
                    output_search_folder, 'search_result.json')))
                with open(str(json_path), 'w') as outfile:
                    json.dump(results, outfile)
                logging.info(
                    'Writing search results in file '+output_search_folder+'/search_result.json')

            html_file = writeHtmlS2(image_url_list, image_title_list,
                                    output_search_folder, s2_product, s2_quicklook_server)

            if quicklook == '1':  # quicklook file is created but it only is opened when quicklook variable = 1
                webbrowser.open(html_file)

            return s2_scenes_filtered, s2_titles  # s2_scenes
    else:
        # print('Searching for S2 images...'+'\n')
        logging.info('Searching for S2 images...'+'\n')
        if scene_sentinel_list != 'NONE':
            # scene_sentinel_list = scene_sentinel_list.split(',')
            s2_scenes, s2_titles, s2_quicklooks, s2_clouds, txt_string, s2_availability = searchS2ByTiles(
                scene_sentinel_list, s2_product, start_date, end_date, max_cloud_cover, quicklook, user_esa, pass_esa, output_res)
            if len(s2_scenes) == 0:
                # print('There are not products S2 to download.'+'\n')
                if txt_string == '':
                    logging.warning(
                        'There are no products S2 to download.'+'\n')
                    sys.exit(1)
                else:
                    logging.warning(
                        'There are no products S2 online to download.'+'\n')
                    if output_res == '1':
                        txt_filename = str(pathlib.Path(os.path.join(
                            os.getcwd(), output_search_folder, 'search_result.txt')))
                        with open(txt_filename, 'w') as txt_file:
                            txt_file.write(txt_string)
                        logging.info(
                            'Writing search results in file ${SAET_HOME}/'+output_search_folder+'/search_result.txt')
                    if output_res == '2':
                        results = exportSentinelResultsToJson(txt_string)
                        json_path = str(pathlib.Path(os.path.join(
                            output_search_folder, 'search_result.json')))
                        with open(str(json_path), 'w') as outfile:
                            json.dump(results, outfile)
                        logging.info(
                            'Writing search results in file ${SAET_HOME}/'+output_search_folder+'/search_result.json')
                    sys.exit(1)
            else:
                txt_string += '\n'
                flag = 0
                #ni = 0
                image_url_list = []
                image_title_list = []
                for i in range(0, len(s2_scenes)):
                    s2_scenes_filtered[str(ni)] = s2_scenes[i]
                    title = s2_titles[i]
                    scene_date = title.split('_')[2][0:8]
                    date_difference = datetime.datetime.strptime(
                        scene_date, "%Y%m%d") - datetime.datetime.strptime(central_date, "%Y%m%d")
                    if central_date < scene_date:
                        if output_res == '0':
                            print(
                                f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                        else:
                            txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                    else:
                        if flag == 0:
                            if output_res == '0':
                                print('[*******] Central date:'+central_date)
                                print(
                                    f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                            else:
                                txt_string += '[*******] Central date:' + \
                                    central_date+'\n'
                                txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                            flag = 1
                        else:
                            if output_res == '0':
                                print(
                                    f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}")
                            else:
                                txt_string += f"[{ni}] Scene: {title} Cloud coverage: {s2_clouds[i]}% {date_difference.days} days {s2_availability[i]}"+'\n'
                    image_url_list.append(s2_quicklooks[i])
                    image_title_list.append(
                        title+' CC:' + str(s2_clouds[i]) + ' Days: ' + str(date_difference.days) + ' ' + s2_availability[i])
                    ni += 1

                if flag == 0:
                    if output_res == '0':
                        print('[*******] Central date: '+central_date)
                    else:
                        txt_string += '[*******] Central date: ' + \
                            central_date+'\n'
                print('\n')
                if output_res == '1':
                    txt_filename = str(pathlib.Path(os.path.join(
                        os.getcwd(), output_search_folder, 'search_result.txt')))
                    with open(txt_filename, 'w') as txt_file:
                        txt_file.write(txt_string)
                    logging.info(
                        'Writing search results in file ${SAET_HOME}/'+output_search_folder+'/search_result.txt')
                if output_res == '2':
                    results = exportSentinelResultsToJson(txt_string)
                    json_path = str(pathlib.Path(os.path.join(
                        output_search_folder, 'search_result.json')))
                    with open(str(json_path), 'w') as outfile:
                        json.dump(results, outfile)
                    logging.info(
                        'Writing search results in file ${SAET_HOME}/'+output_search_folder+'/search_result.json')

                html_file = writeHtmlS2(image_url_list, image_title_list,
                                        output_search_folder, s2_product, s2_quicklook_server)

                if quicklook == '1':  # quicklook file is created but it only is opened when quicklook variable = 1
                    webbrowser.open(html_file)

                return s2_scenes_filtered, s2_titles  # s2_scenes
        else:
            # print('If there is not footprint defined, you must define Sentinel tiles list')
            logging.warning(
                'If there is no footprint defined, you must define Sentinel tiles list')
            sys.exit(1)


def unzipFile(filepath, outpath):
    '''
    Description:
    ------------
    Unzips .zip files to an specific output folder.
    Removes the .zip file after uncompress process

    Arguments:
    ------------
    - filepath (string): path to the .zip file
    - outpath (string): path to the output folder

    Returns:
    ------------
    None

    '''

    with zipfile.ZipFile(filepath) as file:
        file.extractall(outpath)
    os.remove(filepath, dir_fd=None)


def untarFile(filepath, outpath):
    '''
    Description:
    ------------
    Untars .tar and .tar.gz files to an specific output folder.
    Removes the .tar or .tar.gz file after uncompress process

    Arguments:
    ------------
    - filepath (string): path to the .zip file
    - outpath (string): path to the output folder

    Returns:
    ------------
    None

    '''
    with tarfile.open(filepath) as file:
        file.extractall(outpath)
    os.remove(filepath, dir_fd=None)


def getLandsatQuickLookUrl(img_id, producttype):
    '''
    Description:
    ------------
    Makes a valid url for dispaying the quiklook image on an standard browser

    Arguments:
    ------------
    - img_id (string): image id
    - producttype (string): product type for landsat (landsat_8_c1 / landsat_ot_c2_l1)

    Returns:
    ------------
    LandsatQuickLookUrl (string): valid url
    '''
    # valid for collections 1 and 2 (2 -> L8 and L9 from Oct 30th 2021)
    if producttype == 'landsat_8_c1':
        url_base = 'https://ims.cr.usgs.gov/browse/'+producttype
        sections = img_id.split('_')
        img_year = sections[3][0:4]
        img_path = sections[2][0:3]
        img_row = sections[2][3:]
        LandsatQuickLookUrl = url_base+'/'+img_year + \
            '/'+img_path+'/'+img_row+'/'+img_id+'.jpg'
    else:  # landsat_ot_c2_l1
        # url_base = 'https://landsatlook.usgs.gov/gen-browse?size=rrb&type=refl&product_id='
        # LandsatQuickLookUrl = url_base+img_id.replace('_01_', '_02_')
        url_base = 'https://earthexplorer.usgs.gov/index/resizeimage?img=https%3A%2F%2Flandsatlook.usgs.gov%2Fgen-browse%3Fsize%3Drrb%26type%3Drefl%26product_id%3Dimg_id&angle=0&size=640'
        #LandsatQuickLookUrl = url_base+img_id
        LandsatQuickLookUrl = url_base.replace('img_id', img_id)
    return LandsatQuickLookUrl


def searchProductInFolder(data_path):
    '''
    Description:
    ------------
    Search for image products inside a folder. Only for mode "reprocessing downloaded images [op]"

    Arguments:
    ------------
    data_path (string): path to the image data folder

    Returns:
    ------------
    - list of paths

    '''

    try:
        list_of_folders = [pathlib.Path(f.path)
                           for f in os.scandir(data_path) if f.is_dir()]
        return list_of_folders
    except:
        return []


def filterScenesInfolder(list_of_paths):
    '''
    Description:
    ------------
    Filter the number of scene from the list of scenes. Only for reprocessing mode

    Arguments:
    ------------
    list_of_paths (list[string]): paths to the products folder

    Returns:
    ------------
    - list of filtered scenes

    '''
    if len(list_of_paths) != 0:
        list_of_folders = []
        print('List of scenes in the data folder: ')
        print('')
        for i in range(0, len(list_of_paths)):
            path = pathlib.Path(list_of_paths[i])
            folder = path.name
            print([i], folder)
            list_of_folders.append(folder)
        print('')
        np = input('Number of images to be reprocessed?: ')
        filtered_scenes = parseNumberProductsFolder(np, list_of_folders)
        return filtered_scenes
    else:
        return []


# *******************************************************************************
# SECTION: IMAGE PROCESSING FUNCTIONS
# *******************************************************************************

def downScaling(input_file):
    '''
    Description:
    ------------
    Reduces the image spatial resolution from 10 m. to 20 m.
    Uses gdal.RegenerateOverviews() function.

    Arguments:
    ------------
    - input_file (string): path to image input file

    Returns:
    ------------
    - final (numpy matrix): output image downscaled

    '''

    factor = 2  # ratio to reduce resolution from 10 m to 20 m.
    input_name = os.path.basename(input_file)
    output_path = os.path.dirname(input_file)
    output_name = input_name.split('.')[0]+'_20.tif'
    output_file = str(pathlib.Path(os.path.join(output_path, output_name)))
    logging.info('Downscaling '+input_name+' ...')
    # print('Downscaling '+input_name+' ...')
    g = gdal.Open(input_file, gdal.GA_ReadOnly)
    total_obs = g.RasterCount
    drv = gdal.GetDriverByName("MEM")
    dst_ds = drv.Create("", g.RasterXSize, g.RasterYSize, 1, gdal.GDT_UInt16)
    dst_ds.SetGeoTransform(g.GetGeoTransform())
    dst_ds.SetProjection(g.GetProjectionRef())
    dst_ds.GetRasterBand(1).WriteArray(g.ReadAsArray())

    geoT = g.GetGeoTransform()
    drv = gdal.GetDriverByName("GTiff")
    resampled = drv.Create(output_file, int(
        g.RasterXSize/factor), int(g.RasterYSize/factor), 1, gdal.GDT_UInt16)

    this_geoT = (geoT[0], geoT[1]*factor, geoT[2],
                 geoT[3], geoT[4], geoT[5]*factor)
    resampled.SetGeoTransform(this_geoT)
    resampled.SetProjection(g.GetProjectionRef())
    resampled.SetMetadata({"TotalNObs": "%d" % total_obs})
    gdal.RegenerateOverviews(dst_ds.GetRasterBand(
        1), [resampled.GetRasterBand(1)], 'average')
    resampled.GetRasterBand(1).SetNoDataValue(0)
    resampled.FlushCache()
    final = resampled.GetRasterBand(1).ReadAsArray()
    resampled = None
    return final


def createCloudMaskL8(scene_path, cloud_mask_level):
    '''
    Description:
    ------------
    Creates binary cloud mask image for L8 according the image of
    cloud classification (band BQA).
    Saves the cloud mask to the processing folder (folder "temp"
    relative to the scene path).

    Arguments:
    ------------
    - scene_path (string): path to scene folder
    - cloud_mask_level (string): 0, 1 or 2. Level of cloud masking

    Returns:
    ------------
    None

    '''
    # create temp folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)
    mask_values = []

    # landsat from collection 1
    if '_01_' in scene_path:
        # list of values related with medium or high cloud confidence and cirrus values
        cloud_values = [2800, 2804, 2808, 2812, 6896, 6900, 6904, 6908]
        cloud_shadow_values = [2976, 2980, 2984, 2988, 3008, 3012,
                               3016, 3020, 7072, 7076, 7080, 7084, 7104, 7108, 7112, 7116]
        cirrus_values = [6816, 6820, 6824, 6828, 6848, 6852, 6856, 6860, 6896, 6900, 6904, 6908, 7072, 7076, 7080, 7084, 7104,
                         7108, 7112, 7116, 7840, 7844, 7848, 7852, 7872, 7876, 7880, 7884]

        if cloud_mask_level == '0':
            mask_values = [-1]
        if cloud_mask_level == '1':
            mask_values = mask_values+cloud_values
        if cloud_mask_level == '2':
            mask_values = cloud_values+cirrus_values+cloud_shadow_values

        qa_band_path = getBandPath(scene_path, 'BQA')

    if '_02_' in scene_path:
        # list of values related with medium or high cloud confidence and cirrus values
        cloud_values = [22280, 24082, 22080]
        cloud_shadow_values = [23888, 23826, 24144]
        cirrus_values = [55052, 56854]

        if cloud_mask_level == '0':
            mask_values = [-1]
        if cloud_mask_level == '1':
            mask_values = mask_values+cloud_values
        if cloud_mask_level == '2':
            mask_values = cloud_values+cirrus_values+cloud_shadow_values

        qa_band_path = getBandPath(scene_path, 'QA_PIXEL')

    qa_ds = gdal.Open(qa_band_path)
    qa_data = qa_ds.GetRasterBand(1).ReadAsArray()

    cloud_mask = np.isin(qa_data, mask_values)
    # remove objects with less than 10 connected pixels
    morphology.remove_small_objects(
        cloud_mask, min_size=10, connectivity=1, in_place=True)

    path = pathlib.PurePath(scene_path)
    # name of the binary cloud mask
    name_cloud_mask = path.name+'_cmask.tif'
    # image template to copy resolution, bounding box and coordinate reference system.
    base_path = getBandPath(scene_path, 'B6')
    outFileName = str(pathlib.Path(
        os.path.join(output_folder, name_cloud_mask)))
    saveMask(cloud_mask, outFileName, base_path)


def createCloudMaskS2(scene_path, cloud_mask_level):
    '''
    Description:
    ------------
    Creates binary cloud mask image for S2 according the image of
    cloud classification (QA60.gml or MSK_CLASSI_B00.jp2 file for S2MSI1C and SCL band for S2MSI2A).
    Saves the cloud mask to the processing folder (folder "temp"
    relative to the scene path).

    WARNING: The cloud classification for 2A product is more much accurate and
    confident than the cloud classification for the product 1C. It seems that
    this situation can change on 26th of october of 2021.
    More information:
    https://sentinels.copernicus.eu/web/sentinel/-/copernicus-sentinel-2-major-products-upgrade-upcoming

    Arguments:
    ------------
    - scene_path (string): path to scene folder
    - cloud_mask_level (string): 0, 1 or 2. Level of cloud masking

    Returns:
    ------------
    None

    '''

    # create temp folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)
    mask_values = []

    # compute cloud mask for S2MSI1C product
    if 'MSIL1C' in scene_path:
        # in this case, the .gml file (vector) has to be converted to .shp file
        # and then to a binary raster image.
        gml_path = getBandPath(scene_path, 'QA60.gml')
        msk_classi_path = getBandPath(scene_path, 'B00')
        if gml_path is None:  # cloud mask in raster format
            path = pathlib.PurePath(scene_path)
            raster_path = path.name+'_cmask.tif'
            output_path = str(pathlib.Path(
                os.path.join(output_folder, raster_path)))
            raster_template = getBandPath(scene_path, 'B11')

            if cloud_mask_level == '1' or cloud_mask_level == '2':
                mask_values = [1]
                msk_classi_path = getBandPath(scene_path, 'CPM.jp2')
                msk_ds = gdal.Open(msk_classi_path)
                msk_data = msk_ds.GetRasterBand(1).ReadAsArray()
                cloud_mask = np.isin(msk_data, mask_values)
                name_cloud_mask = path.name+'_cmask.tif'
                # image template to copy resolution, bounding box and coordinate reference system.
                base_path = getBandPath(scene_path, 'B11')
                outFileName = str(pathlib.Path(
                    os.path.join(output_folder, name_cloud_mask)))
                saveMask(cloud_mask, outFileName, base_path)
            else:
                createEmptyCloudMaskS2(output_path, raster_template)
        else:  # cloud mask in gml format
            path = pathlib.PurePath(scene_path)
            shp_path = path.name+'_cmask.shp'
            raster_path = path.name+'_cmask.tif'
            output_path = str(pathlib.Path(
                os.path.join(output_folder, raster_path)))
            raster_template = getBandPath(scene_path, 'B11')
            output_shp = str(pathlib.Path(
                os.path.join(output_folder, shp_path)))

            if cloud_mask_level == '1' or cloud_mask_level == '2':
                # prevents .gml files without clouds
                res = gmlToShp(gml_path, output_shp)
                if res:
                    rasterizeShapefile(
                        output_shp, output_path, raster_template)
                else:
                    createEmptyCloudMaskS2(output_path, raster_template)
            else:
                createEmptyCloudMaskS2(output_path, raster_template)

    # compute cloud mask for S2MSI2A product
    if 'MSIL2A' in scene_path:
        # cloud codes: 2 -> dark area; 3 -> cloud shadow; 8 -> cloud medium; 9 -> cloud high; 10 -> cirrus

        if cloud_mask_level == '0':
            mask_values = [-1]
        if cloud_mask_level == '1':
            mask_values = [9]
        if cloud_mask_level == '2':
            mask_values = [3, 8, 9, 10]

        scl_band_path = getBandPath(scene_path, 'SCL')
        scl_ds = gdal.Open(scl_band_path)
        scl_data = scl_ds.GetRasterBand(1).ReadAsArray()

        cloud_mask = np.isin(scl_data, mask_values)
        # remove objects with less than 10 connected pixels
        morphology.remove_small_objects(
            cloud_mask, min_size=10, connectivity=1, in_place=True)

        path = pathlib.PurePath(scene_path)
        # name of the binary cloud mask
        name_cloud_mask = path.name+'_cmask.tif'
        # image template to copy resolution, bounding box and coordinate reference system.
        base_path = getBandPath(scene_path, 'B11')
        outFileName = str(pathlib.Path(
            os.path.join(output_folder, name_cloud_mask)))
        saveMask(cloud_mask, outFileName, base_path)


def createEmptyCloudMaskS2(output_path, raster_template):
    '''
    Description:
    ------------
    Creates an image of 0 values with the same features of the input
    template image (resolution, bounding box and coordinate reference system).

    Arguments:
    ------------
    - output_path (string): path to the output file
    - raster_template (string): path of the image template

    Returns:
    ------------
    None

    '''
    ds_in = gdal.Open(raster_template)
    band_in = ds_in.GetRasterBand(1)
    data_in = band_in.ReadAsArray()
    trs = ds_in.GetGeoTransform()
    rows, cols = data_in.shape
    data_out = np.zeros((rows, cols)).astype(np.uint8)
    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(output_path, cols, rows, 1, gdal.GDT_Byte)
    ds_out.SetGeoTransform(trs)
    ds_out.SetProjection(ds_in.GetProjection())
    ds_out.GetRasterBand(1).WriteArray(data_out)
    ds_out.FlushCache()
    ds_in = None
    band_in = None
    ds_out = None


def aweishLandsat(scene_path):
    '''
    Description:
    ------------
    Computes L8 aweish water index for the analysis area. Involved bands: B2(blue)
    B3(green), B5(nir), B6(swir1), B7(swir2)

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B6')

    # getting bands data
    band_blue = gdal.Open(getBandPath(scene_path, 'B2'))
    data_blue = band_blue.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_green = gdal.Open(getBandPath(scene_path, 'B3'))
    data_green = band_green.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_nir = gdal.Open(getBandPath(scene_path, 'B5'))
    data_nir = band_nir.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B6'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir2 = gdal.Open(getBandPath(scene_path, 'B7'))
    data_swir2 = band_swir2.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # getting parameters to convert Dns to TOA values from MTL file
    rmb, rab, se = getTOAParameters(getBandPath(scene_path, 'MTL.txt'), '1')

    # DNs to TOA conversion
    se_factor = sin(se*(np.pi/180.0))
    db_toa = (data_blue*rmb+rab) / se_factor
    dg_toa = (data_green*rmb+rab) / se_factor
    dn_toa = (data_nir*rmb+rab) / se_factor
    ds1_toa = (data_swir1*rmb+rab) / se_factor
    ds2_toa = (data_swir2*rmb+rab) / se_factor

    # computing water index
    aweish = db_toa + (2.5 * dg_toa) - \
        (1.5 * (dn_toa + ds1_toa)) - (0.25 * ds2_toa)

    # saving water index
    saveIndex(aweish, outFileName, base_path)


def aweinshLandsat(scene_path):
    '''
    Description:
    ------------
    Computes L8 aweinsh water index for the analysis area. Involved bands: B3(green)
    B5(nir), B6(swir1), B7(swir2)

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B6')

    # getting bands data
    band_green = gdal.Open(getBandPath(scene_path, 'B3'))
    data_green = band_green.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_nir = gdal.Open(getBandPath(scene_path, 'B5'))
    data_nir = band_nir.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B6'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir2 = gdal.Open(getBandPath(scene_path, 'B7'))
    data_swir2 = band_swir2.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # getting parameters to convert Dns to TOA values from MTL file
    rmb, rab, se = getTOAParameters(getBandPath(scene_path, 'MTL.txt'), '1')

    # DNs to TOA conversion
    se_factor = sin(se*(np.pi/180.0))
    dg_toa = (data_green*rmb+rab) / se_factor
    dn_toa = (data_nir*rmb+rab) / se_factor
    ds1_toa = (data_swir1*rmb+rab) / se_factor
    ds2_toa = (data_swir2*rmb+rab) / se_factor

    # computing water index
    aweinsh = 4 * (dg_toa - ds1_toa) - (0.25 * dn_toa + 2.75 * ds2_toa)

    # saving water index
    saveIndex(aweinsh, outFileName, base_path)


def mndwiLandsat(scene_path):
    '''
    Description:
    ------------
    Computes L8 mndwi water index for the analysis area. Involved bands: B3(green)
    B6(swir1)

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''
    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B6')

    # getting bands data
    band_green = gdal.Open(getBandPath(scene_path, 'B3'))
    data_green = band_green.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B6'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # getting parameters to convert Dns to TOA values from MTL file
    rmb, rab, se = getTOAParameters(getBandPath(scene_path, 'MTL.txt'), '1')

    # DNs to TOA conversion
    se_factor = sin(se*(np.pi/180.0))
    dg_toa = (data_green*rmb+rab) / se_factor
    ds1_toa = (data_swir1*rmb+rab) / se_factor

    # computing water index
    mndwi = (dg_toa - ds1_toa) / (dg_toa + ds1_toa)

    # saving water index
    saveIndex(mndwi, outFileName, base_path)


def aweinshS2(scene_path):
    '''
    Description:
    ------------
    Computes S2 aweinsh water index for the analysis area. Involved bands: B03(green)
    B08(nir), B11(swir1), B12(swir2).
    Downscales some bands if it is needed to homogenize spatial resolutions

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B11')

    # getting bands data and downscaling if it is needed
    band_green = gdal.Open(getBandPath(scene_path, 'B03'))
    pix_size = band_green.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_green = downScaling(getBandPath(
            scene_path, 'B03')).astype(np.float32)
    else:
        data_green = band_green.GetRasterBand(
            1).ReadAsArray().astype(np.float32)

    band_nir = gdal.Open(getBandPath(scene_path, 'B08'))
    pix_size = band_nir.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_nir = downScaling(getBandPath(
            scene_path, 'B08')).astype(np.float32)
    else:
        data_nir = band_nir.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B11'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir2 = gdal.Open(getBandPath(scene_path, 'B12'))
    data_swir2 = band_swir2.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # computing water index
    aweinsh = 4 * (data_green - data_swir1) - \
        (0.25 * data_nir + 2.75 * data_swir2)

    # saving water index
    saveIndex(aweinsh, outFileName, base_path)


def aweishS2(scene_path):
    '''
    Description:
    ------------
    Computes S2 aweish water index for the analysis area. Involved bands: B02(blue),
    B03(green), B08(nir), B11(swir1), B12(swir2).
    Downscales some bands if it is needed to homogenize spatial resolutions

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B11')

    # getting bands data and downscaling if it is needed
    band_blue = gdal.Open(getBandPath(scene_path, 'B02'))
    pix_size = band_blue.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_blue = downScaling(getBandPath(
            scene_path, 'B02')).astype(np.float32)
    else:
        data_blue = band_blue.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_green = gdal.Open(getBandPath(scene_path, 'B03'))
    pix_size = band_green.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_green = downScaling(getBandPath(
            scene_path, 'B03')).astype(np.float32)
    else:
        data_green = band_green.GetRasterBand(
            1).ReadAsArray().astype(np.float32)

    band_nir = gdal.Open(getBandPath(scene_path, 'B08'))
    pix_size = band_nir.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_nir = downScaling(getBandPath(
            scene_path, 'B08')).astype(np.float32)
    else:
        data_nir = band_nir.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B11'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    band_swir2 = gdal.Open(getBandPath(scene_path, 'B12'))
    data_swir2 = band_swir2.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # computing water index
    aweish = data_blue + (2.5 * data_green) - \
        (1.5 * (data_nir + data_swir1)) - (0.25 * data_swir2)

    # saving water index
    saveIndex(aweish, outFileName, base_path)


def mndwiS2(scene_path):
    '''
    Description:
    ------------
    Computes S2 mndwi water index for the analysis area. Involved bands: B03(green),
    B11(swir1).
    Downscales some bands if it is needed to homogenize spatial resolutions

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # prevents numpy errors for invalid values or divide by zero
    np.seterr(divide='ignore', invalid='ignore')

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_wi.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B11')

    # getting bands data and downscaling if it is needed
    band_green = gdal.Open(getBandPath(scene_path, 'B03'))
    pix_size = band_green.GetGeoTransform()[1]
    if pix_size == 10.0:
        data_green = downScaling(getBandPath(
            scene_path, 'B03')).astype(np.float32)
    else:
        data_green = band_green.GetRasterBand(
            1).ReadAsArray().astype(np.float32)

    band_swir1 = gdal.Open(getBandPath(scene_path, 'B11'))
    data_swir1 = band_swir1.GetRasterBand(1).ReadAsArray().astype(np.float32)

    # computing water index
    mndwi = (data_green - data_swir1) / (data_green + data_swir1)

    # saving water index
    saveIndex(mndwi, outFileName, base_path)


def getIndexMask(index_path, thr_method, tol_area=300):
    '''
    Description:
    ------------
    Computes binary mask from water index using the standar value 0 for
    segmentation

    Arguments:
    ------------
    - index_path (string): path to water index
    - thr_method (string): method to segmentation threshold computation
      {'0': standard zero, '1': otsu bimodal, '2': otsu multimodal with 3 clases}
    - tol_area (int): tolerance to remove small holes. Default: 300

    Returns:
    ------------
    - imgmask (numpy matrix): if area removing is enabled
    - index_copy (numpy matrix): if area removing is disabled

    '''

    index_ds = gdal.Open(index_path)
    band = index_ds.GetRasterBand(1)
    index_data = band.ReadAsArray()
    index_data[index_data == float('-inf')] = 0.

    # tolerance for segmentation
    if thr_method == '0':
        tol = 0
    if thr_method == '1':
        cimg = index_data.copy()
        vec_cimg = cimg.reshape(cimg.shape[0]*cimg.shape[1])
        tol = threshold_otsu(vec_cimg)
    if thr_method == '2':
        th_otsu_multi = threshold_multiotsu(index_data, 3)
        if abs(th_otsu_multi[0]) < abs(th_otsu_multi[1]):
            tol = th_otsu_multi[0]
        else:
            tol = th_otsu_multi[1]

    # image binarization according threshold
    index_copy = index_data.copy()
    index_copy[index_data < tol] = 0.  # land
    index_copy[index_data >= tol] = 1.  # water

    if tol_area != 0:  # area removing
        img_mask = removeHolesByArea(index_copy.astype(np.byte), tol_area)
        return img_mask
    else:
        return index_copy


def createPixelLine(method, mask, cmask):
    '''
    Description:
    ------------
    Computes binary pixel mask for rough shoreline from water index mask.
    Removes clouds areas from binary pixel mask.

    Arguments:
    ------------
    - method (string): erosion or dilation
    - mask (numpy matrix): binary water index mask
    - cmask (numpy matrix): binary cloud mask

    Returns:
    ------------
    - pixel_line (numpy matrix): binary rough pixel shoreline without cloud areas

    '''

    # getting cloud mask data
    cmask_ds = gdal.Open(cmask)
    cmask_band = cmask_ds.GetRasterBand(1)
    cmask_data = cmask_band.ReadAsArray()

    # kernel for cloud mask buffering
    kernel = np.ones((9, 9), np.uint8)

    # getting pixel shoreline mask
    if method == 'erosion':
        erosion = morphology.binary_erosion(mask)
        pixel_line = mask-erosion
    if method == 'dilation':
        dilation = morphology.binary_dilation(mask)
        pixel_line = dilation-mask

    # cleaning pixel line mask using buffer of cloud areas
    cmask_dilation = morphology.binary_dilation(cmask_data, kernel)
    pixel_line[cmask_dilation == 1] = 0
    pixel_line[pixel_line == 1] = 255
    cmask_data = None

    return pixel_line


def removeHolesByArea(img, area):
    '''
    Description:
    ------------
    Removes litle holes from binary images.

    Arguments:
    ------------
    - img (numpy matrix): input image
    - area (int): area tolerance for connected pixels

    Returns:
    ------------
    - pixel_line (numpy matrix): binary rough pixel shoreline without cloud areas

    '''

    img_closed = morphology.area_closing(img, area, connectivity=1)
    img_closed = morphology.area_closing(~img_closed, area, connectivity=1)
    return ~img_closed


def maskPixelLine(pixel_line_path, mask_path):
    '''
    Description:
    ------------
    Masks binary rough pixel shoreline with beaches mask. All pixels outside of
    a beach area will be removed (converted to value 0).

    Arguments:
    ------------
    - pixel_line_path (string): path to the pixel line mask
    - mask_path (string): path to the binary beaches mask

    Returns:
    ------------
    None

    '''

    pl_ds = gdal.Open(pixel_line_path, gdal.GA_Update)
    pl_band = pl_ds.GetRasterBand(1)
    pl_data = pl_band.ReadAsArray()

    mask_ds = gdal.Open(mask_path)
    mask_band = mask_ds.GetRasterBand(1)
    mask_data = mask_band.ReadAsArray()

    pl_data[mask_data == 0] = 0
    pl_band.WriteArray(pl_data)
    pl_ds.FlushCache()
    pl_ds = None


def saveIndex(in_array, out, template_path, dType=gdal.GDT_Float32):
    '''
    Description:
    ------------
    Saves water index to tiff image

    Arguments:
    ------------
    - in_array (numpy matrix): water index data
    - out (string): output path to the tiff image
    - template_path (string): template image to copy resolution, bounding box
    and coordinate reference system.
    - dType: data type format (default: float 32 bits)

    Returns:
    ------------
    None

    '''

    if os.path.exists(out):
        os.remove(out)

    template = gdal.Open(template_path)
    driver = gdal.GetDriverByName('GTiff')
    shape = in_array.shape
    dst_ds = driver.Create(
        out, xsize=shape[1], ysize=shape[0], bands=1, eType=dType)
    proj = template.GetProjection()
    geo = template.GetGeoTransform()
    dst_ds.SetGeoTransform(geo)
    dst_ds.SetProjection(proj)
    dst_ds.GetRasterBand(1).WriteArray(in_array)
    dst_ds.FlushCache()
    dst_ds = None


def saveMask(img_out, outFilename, base_path):
    '''
    Description:
    ------------
    Saves binary mask to tiff image (byte data type)

    Arguments:
    ------------
    - img_out(numpy matrix): binary mask data
    - out_filename (string): output path to the tiff image
    - base_path (string): template image to copy resolution, bounding box
    and coordinate reference system.

    Returns:
    ------------
    None

    '''

    g = gdal.Open(base_path)
    geoTransform = g.GetGeoTransform()
    geoProjection = g.GetProjection()
    driver = gdal.GetDriverByName("GTiff")
    newDataset = driver.Create(outFilename, g.RasterXSize, g.RasterYSize,
                               1, gdal.GDT_Byte, options=["COMPRESS=Deflate"])
    newDataset.SetGeoTransform(geoTransform)
    newDataset.SetProjection(geoProjection)
    newDataset.GetRasterBand(1).WriteArray(img_out.astype(np.uint8))
    newDataset.FlushCache()
    newDataset = None


def getTOAParameters(mtl, band):
    '''
    Description:
    ------------
    Get TOA parameters from MTL metadata file

    Arguments:
    ------------
    - mtl (string): path to the MTL file
    - band (string): band name

    Returns:
    ------------
    - rmb (float): reflectance value to multiply for a single band
    - rab (float): reflectance value to add for a single band
    - se (float): sun elevation for the image scene

    '''

    with open(mtl, "r") as mtl:
        for line in mtl:
            if "REFLECTANCE_MULT_BAND_"+band in line:
                rmb = float(line.strip().split("=")[1].strip())
            if "REFLECTANCE_ADD_BAND_"+band in line:
                rab = float(line.strip().split("=")[1].strip())
            if "SUN_ELEVATION" in line:
                se = float(line.strip().split("=")[1].strip())
    return rmb, rab, se


def getColors(nc):
    '''
    makes one array os nc colors between 0 and 255
    '''
    colors = []
    salto = int(255/nc)
    for i in range(0, 255, salto+1):
        colors.append(i)
    return colors


def getWaterClass(Z, km_classes, colors):
    '''
    returns the water class based on the kmeans class with the minimum mean value
    on the swir1 band
    '''
    means = []
    for c in colors:
        means.append(np.mean(Z[km_classes == c]))
    water_class = np.argmin(means)
    return colors[water_class]


def computeKmeansS2(scene_path):
    '''
    Description:
    ------------
    Computes kmeans clasterization method. Involved bands: Swir1 (B11)

    The number of classes by default are 3

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''
    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_kmeans_mask.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B11')

    # open image
    img = gdal.Open(base_path)
    img_data = img.GetRasterBand(1).ReadAsArray()
    img_data[img_data == float('-inf')] = 0.0

    # compute percentiles for image contrast enhancing
    p1, p2 = np.percentile(img_data, (0.5, 99.5))
    img_data = exposure.rescale_intensity(img_data, in_range=(p1, p2))

    # convert image to only one row array
    w, h = img_data.shape
    Z = img_data.reshape((-1, 1))
    Z = np.float32(Z)

    # compute kmeans classification
    number_clusters = 3
    km = KMeans(number_clusters)
    km.fit(Z)
    labels = km.labels_

    # assign codes (colors) to each clase
    colors = getColors(number_clusters)
    km_classes = np.zeros(w*h, dtype='uint8')
    for ix in range(km_classes.shape[0]):
        km_classes[ix] = colors[labels[ix]]

    # get water class
    water_class = getWaterClass(Z, km_classes, colors)

    # compute and save mask
    km_classes = km_classes.reshape((w, h))
    binary_mask = np.where(km_classes == water_class, 1, 0)
    saveMask(binary_mask, outFileName, base_path)


def computeKmeansLandsat(scene_path):
    '''
    Description:
    ------------
    Computes kmeans clasterization method. Involved bands: Swir1 (B6)

    The number of classes by default are 3

    Arguments:
    ------------
    - scene_path (string): path to the scene folder

    Returns:
    ------------
    None

    '''

    # create output folder if it is needed
    output_folder = str(pathlib.Path(os.path.join(scene_path, 'temp')))
    createFolderCheck(output_folder)

    # output file name setting
    path = pathlib.PurePath(scene_path)
    name = path.name+'_kmeans_mask.tif'
    outFileName = str(pathlib.Path(os.path.join(output_folder, name)))

    # template image to copy resolution, bounding box and coordinate reference system
    base_path = getBandPath(scene_path, 'B6')

    # open image
    img = gdal.Open(base_path)
    img_data = img.GetRasterBand(1).ReadAsArray()
    img_data[img_data == float('-inf')] = 0.0

    # compute percentiles for image contrast enhancing
    p1, p2 = np.percentile(img_data, (0.5, 99.5))
    img_data = exposure.rescale_intensity(img_data, in_range=(p1, p2))

    # convert image to only one row array
    w, h = img_data.shape
    Z = img_data.reshape((-1, 1))
    Z = np.float32(Z)

    # compute kmeans classification
    number_clusters = 3
    km = KMeans(number_clusters)
    km.fit(Z)
    labels = km.labels_

    # assign codes (colors) to each clase
    colors = getColors(number_clusters)
    km_classes = np.zeros(w*h, dtype='uint8')
    for ix in range(km_classes.shape[0]):
        km_classes[ix] = colors[labels[ix]]

    # get water class
    water_class = getWaterClass(Z, km_classes, colors)

    # compute and save mask
    km_classes = km_classes.reshape((w, h))
    binary_mask = np.where(km_classes == water_class, 1, 0)
    saveMask(binary_mask, outFileName, base_path)

# *******************************************************************************
# SECTION: HELPER FUNCTIONS
# *******************************************************************************


def getBandPath(scene_path, band_name):
    '''
    Description:
    ------------
    Get absolute path from a single band or file

    Arguments:
    ------------
    - scene_path (string): path to the target folder
    - band_name (string): band name to search

    Returns:
    ------------
    - band_path (string): path to the band

    '''
    file_list = recursiveFileSearch(scene_path, '*.*')
    band_path = [i for i in file_list if (
        band_name in i) and (not 'xml' in i)]
    if len(band_path) != 0:
        return str(pathlib.Path(band_path[0]))
    else:
        return None


def getBandData(band_path):
    '''
    Returns the data matrix from a band path
    '''
    band = gdal.Open(band_path)
    band_data = band.GetRasterBand(1).ReadAsArray()
    return band_data


def getLandsatBbox(footprint):
    '''
    Description:
    ------------
    Converts roi footprint in wkb format to (xmin,ymin,xmax,ymax) tuple

    Arguments:
    ------------
    - footprint (string): roi footprint in wkb format

    Returns:
    ------------
    - bbox (tuple of strings): tuple of 4 strings

    '''

    coords = footprint.split('((')[1].split('))')[0].split(',')
    coords_3 = coords[3].split(' ')
    xmin = float(coords_3[0])
    ymin = float(coords_3[1])
    coords_1 = coords[1].split(' ')
    xmax = float(coords_1[0])
    ymax = float(coords_1[1])
    bbox = (xmin, ymin, xmax, ymax)
    return bbox


def getSourceEpsg():
    '''
    Description:
    ------------
    Gets geopgraphic WGS84 coordinate reference system in EPSG format.
    By default, this crs follows the rule (lat,long). It is needed to change
    this to the rule (long,lat).

    Arguments:
    ------------
    None

    Returns:
    ------------
    - source_epsg (object): osr spatial reference object

    '''
    source_epsg = osr.SpatialReference()
    source_epsg.ImportFromEPSG(4326)
    # be careful -> traditional gis order = (long, lat); default = (lat, long)
    source_epsg.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return source_epsg


def getTargetEpsg(scene_path, band_name):
    '''
    Description:
    ------------
    Gets coordinate reference system in EPSG format for a single band

    Arguments:
    ------------
    - scene_path (string): path to the scene folder
    - band_name (string): name of the target band

    Returns:
    ------------
    - source_epsg (object): osr spatial reference object

    '''

    file_list = recursiveFileSearch(scene_path, '*.*')
    band_path = [i for i in file_list if band_name in i][0]
    ds = gdal.Open(band_path)
    epsg_code = osr.SpatialReference(
        wkt=ds.GetProjection()).GetAttrValue('AUTHORITY', 1)
    target_epsg = osr.SpatialReference()
    target_epsg.ImportFromEPSG(int(epsg_code))
    ds = None
    return target_epsg


def reprojectShp(input_shp, output_shp, inSpatialRef, outSpatialRef):
    '''
    Description:
    ------------
    Reprojects one shapefile between two CRSs

    Arguments:
    ------------
    - input_shp (string): path to the input shapefile
    - output_shp (string): path to the output shapefile
    - inSpatialRef (osr spatial reference object): source crs
    - outSpatialRef (osr spatial reference object): target crs

    Returns:
    ------------
    None

    '''
    # getting transformation matrix between source and target crs
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # create output shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    inDataSet = driver.Open(input_shp)
    inLayer = inDataSet.GetLayer()
    inLayer_geomtype = inLayer.GetGeomType()
    inLayer.ResetReading()

    if os.path.exists(output_shp):
        driver.DeleteDataSource(output_shp)
    outDataSet = driver.CreateDataSource(output_shp)
    outLayer = outDataSet.CreateLayer(
        " ", geom_type=inLayer_geomtype)  # ogr.wkbPolygon

    # copy input shapefile database structure to the output shapefile
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    outLayerDefn = outLayer.GetLayerDefn()

    # copy attribute values and geometries from input shapefile
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        geom = inFeature.GetGeometryRef()
        if not geom is None:
            geom.Transform(coordTrans)
            outFeature = ogr.Feature(outLayerDefn)
            outFeature.SetGeometry(geom)
            for i in range(0, outLayerDefn.GetFieldCount()):
                outFeature.SetField(outLayerDefn.GetFieldDefn(
                    i).GetNameRef(), inFeature.GetField(i))
            outLayer.CreateFeature(outFeature)
            outFeature = None
        inFeature = inLayer.GetNextFeature()

    inLayer.ResetReading()
    inLayer = None
    inDataSet.Destroy()
    outLayer = None
    outDataSet.Destroy()

    # create ESRI.prj file
    outSpatialRef.MorphToESRI()
    file_name = os.path.basename(output_shp)
    dir_name = os.path.dirname(output_shp)
    prj_name = str(pathlib.Path(os.path.join(
        dir_name, file_name.split('.')[0]+'.prj')))
    with open(prj_name, 'w') as prj_file:
        prj_file.write(outSpatialRef.ExportToWkt())


def clipShapefile(input_shp, output_shp, clip_shp):
    '''
    Description:
    ------------
    Clips one single shapefile using a second shapefile

    Arguments:
    ------------
    - input_shp (string): path to the input shapefile
    - output_shp (string): path to the output shapefile
    - clip_shp (string): path to the shapefile to clip input shapefile

    Returns:
    ------------
    None

    '''

    driver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = driver.Open(input_shp, 0)
    inLayer = inDataSource.GetLayer()

    inClipSource = driver.Open(clip_shp, 0)
    inClipLayer = inClipSource.GetLayer()

    outDataSource = driver.CreateDataSource(output_shp)
    outLayer = outDataSource.CreateLayer('clip', geom_type=ogr.wkbMultiPolygon)

    ogr.Layer.Clip(inLayer, inClipLayer, outLayer)

    # create ESRI.prj file
    outSpatialRef = inLayer.GetSpatialRef()
    outSpatialRef.MorphToESRI()
    file_name = os.path.basename(output_shp)
    dir_name = os.path.dirname(output_shp)
    prj_name = str(pathlib.Path(os.path.join(
        dir_name, file_name.split('.')[0]+'.prj')))
    with open(prj_name, 'w') as prj_file:
        prj_file.write(outSpatialRef.ExportToWkt())

    inDataSource = None
    inClipSource = None
    outDataSource = None


def createShapefileFromRasterFootprint(raster_path, output_shp, target_epsg, geom_type='polygon'):
    '''
    Description:
    ------------
    Creates a shapefile from a raster footprint

    Arguments:
    ------------
    - raster_path (string): path to the input raster
    - output_shp (string): path to the output shapefile
    - target_epsg (osr spatial reference object): crs for the output shp
    - geom_type (string): type of geometry for the output shapefile

    Returns:
    ------------
    None

    '''
    footprint = getRasterFootprint(raster_path)
    dic_geom = {'polygon': ogr.wkbPolygon,
                'point': ogr.wkbPoint, 'line': ogr.wkbLineString}
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(output_shp)
    layer = data_source.CreateLayer(' ', target_epsg, dic_geom[geom_type])
    layer.CreateField(ogr.FieldDefn("Iden", ogr.OFTInteger))
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("Iden", 1)
    feature.SetGeometry(footprint)
    layer.CreateFeature(feature)
    feature = None
    data_source = None


def getRasterFootprint(raster_path):
    '''
    Description:
    ------------
    Gest raster footprint as polygon geometry in wkb format

    Arguments:
    ------------
    - raster_path (string): path to the input raster

    Returns:
    ------------
    - footprint (string): polygon geometry in wkb format

    '''

    # Get raster geometry
    raster = gdal.Open(raster_path)
    transform = raster.GetGeoTransform()
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    xLeft = transform[0]
    yTop = transform[3]
    xRight = xLeft+cols*pixelWidth
    yBottom = yTop+rows*pixelHeight

    # image footprint
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xLeft, yTop)
    ring.AddPoint(xLeft, yBottom)
    ring.AddPoint(xRight, yBottom)
    ring.AddPoint(xRight, yTop)
    ring.AddPoint(xLeft, yTop)
    footprint = ogr.Geometry(ogr.wkbPolygon)
    footprint.AddGeometry(ring)

    return footprint


def rasterizeShapefile(input_shp, output_raster, raster_template, bc):
    '''
    Description:
    ------------
    Converts input shapefile to raster TIFF file according to
    the raster template features (spatial resolution, bounding box
    and CRS).
    All geometries in the shapefile will be rasterized to create a binary
    raster with values 0 or 1.

    Arguments:
    ------------
    - input_shp (string): path to the input shapefile
    - output_raster (string): path to the output raster file
    - raster_template (string): path to the raster template file
    - bc (string): code to filter especific beach

    Returns:
    ------------
    NONE

    '''

    driver = ogr.GetDriverByName("ESRI Shapefile")
    shp_ds = driver.Open(input_shp, 0)
    template_ds = gdal.Open(raster_template)
    lyr = shp_ds.GetLayer()
    geot = template_ds.GetGeoTransform()
    prj = template_ds.GetProjection()
    driver = gdal.GetDriverByName('GTiff')
    new_raster_ds = driver.Create(
        output_raster, template_ds.RasterXSize, template_ds.RasterYSize, 1, gdal.GDT_Byte)
    new_raster_ds.SetGeoTransform(geot)
    new_raster_ds.SetProjection(prj)
    # filter by beach code if needed
    if bc is None:
        gdal.RasterizeLayer(new_raster_ds, [1], lyr)
    else:
        lyr.SetAttributeFilter("BEACH_CODE IN "+bc)
        if lyr.GetFeatureCount() == 0:
            lyr.SetAttributeFilter('')
        gdal.RasterizeLayer(new_raster_ds, [1], lyr)
    new_raster_ds.GetRasterBand(1).SetNoDataValue(2)
    new_raster_data = new_raster_ds.GetRasterBand(1).ReadAsArray()
    new_raster_data[new_raster_data == 255] = 1
    new_raster_data[new_raster_data != 1] = 0
    new_raster_ds.GetRasterBand(1).WriteArray(new_raster_data)
    new_raster_ds.FlushCache()
    new_raster_ds = None
    new_raster_data = None


def gmlToShp(gml_path, shp_path):
    '''
    Description:
    ------------
    Converts .gml file to .shp file.
    .gml files contains cloud mask for S2MSI1C product

    Arguments:
    ------------
    - gml_path (string): path to the input .gml file
    - shp_path (string): path to the output .shp file

    Returns:
    ------------
    NONE

    '''

    gml_ds = ogr.Open(gml_path)
    gml_layer = gml_ds.GetLayer()
    if gml_layer is None:
        # prevents empty .gml files (no clouds)
        gml_ds.Destroy()
        return False
    else:
        # get projection and database definition from .gml file
        gml_projection = gml_layer.GetSpatialRef()
        gml_layer_defn = gml_layer.GetLayerDefn()

        # creates shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(shp_path):
            driver.DeleteDataSource(shp_path)
        shp_ds = driver.CreateDataSource(shp_path)
        shp_layer = shp_ds.CreateLayer(
            ' ', geom_type=gml_layer_defn.GetGeomType(), srs=gml_projection)
        in_field_count = gml_layer_defn.GetFieldCount()

        # clones database definition
        for fld_index in range(in_field_count):
            src_fd = gml_layer_defn.GetFieldDefn(fld_index)
            fd = ogr.FieldDefn(src_fd.GetName(), src_fd.GetType())
            fd.SetWidth(src_fd.GetWidth())
            fd.SetPrecision(src_fd.GetPrecision())
            shp_layer.CreateField(fd)

        # copy attributte values and geometries
        in_feat = gml_layer.GetNextFeature()
        while in_feat is not None:
            geom = in_feat.GetGeometryRef().Clone()
            out_feat = ogr.Feature(feature_def=shp_layer.GetLayerDefn())

            for fld_index in range(in_field_count):
                src_fd = gml_layer_defn.GetFieldDefn(fld_index)
                name = src_fd.GetName()
                value = in_feat.GetField(fld_index)
                out_feat.SetField(name, value)

            out_feat.SetGeometry(geom)
            shp_layer.CreateFeature(out_feat)
            out_feat.Destroy()
            in_feat.Destroy()
            in_feat = gml_layer.GetNextFeature()

        gml_ds.Destroy()
        shp_ds.Destroy()
        return True


def addIdField(input_path):
    '''
    Description:
    ------------
    Adds the field "ID_FEAT" to the input shapefile

    Arguments:
    ------------
    - input_path (string): path to the input shapefile

    Returns:
    ------------
    None

    '''
    # open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds_shp = driver.Open(input_path, 1)
    layer_shp = ds_shp.GetLayer()

    layer_shp_defn = layer_shp.GetLayerDefn()
    field_names = [layer_shp_defn.GetFieldDefn(
        i).GetName() for i in range(layer_shp_defn.GetFieldCount())]

    # adds id field
    if not 'ID_FEAT' in field_names:  # to ensure that attribute "ID_FEAT" exists
        new_field = ogr.FieldDefn('ID_FEAT', ogr.OFTInteger)
        layer_shp.CreateField(new_field)

    # populates id field
    id_feat_shp = 0
    for feat_shp in layer_shp:
        feat_shp.SetField('ID_FEAT', id_feat_shp)
        layer_shp.SetFeature(feat_shp)
        id_feat_shp += 1

    ds_shp.FlushCache()
    ds_shp = None


def exportToGeojson(shp_path):
    '''
    Description:
    ------------
    Exports shapefile to GeoJson format

    Arguments:
    ------------
    - shp_path (string): path to the input shapefile

    Returns:
    ------------
    None

    '''
    base_name = os.path.basename(shp_path).split('.')[0]
    dir_name = os.path.dirname(shp_path)
    geojson_path = str(pathlib.Path(
        os.path.join(dir_name, base_name+'.json')))
    with shapefile.Reader(shp_path) as shp:
        geojson_data = shp.__geo_interface__
        with open(geojson_path, 'w') as geojson_file:
            geojson_file.write(json.dumps(geojson_data))


def copyShpToFolder(input_path, output_path, source_epsg):
    '''
    Description:
    ------------
    Copies final point and line shapefiles of the extracted shoreline to the
    SDS folder. The copied file is previously reprojected to the WGS84 lat-long
    spatial reference system

    Arguments:
    ------------
    - source_epsg(object): spatial reference system for input shapefile
    - input_path (string): path to the input shapefiles
    - output_path (string): path to the output shapefiles

    Returns:
    ------------
    None

    '''
    cp_files = recursiveFileSearch(input_path, '*_cp.shp')
    cl_files = recursiveFileSearch(input_path, '*_cl.shp')

    root_name = os.path.basename(cp_files[0]).split('_cp.')[0]
    output_folder = str(pathlib.Path(os.path.join(output_path, root_name)))
    createFolderCheck(output_folder)
    target_epsg = getSourceEpsg()

    for cp_file in cp_files:
        addIdField(cp_file)
        reprojectShp(cp_file, str(pathlib.Path(
            os.path.join(output_folder, root_name+'_points.shp'))), source_epsg, target_epsg)
        exportToGeojson(str(pathlib.Path(
            os.path.join(output_folder, root_name+'_points.shp'))))

    for cl_file in cl_files:
        addIdField(cl_file)
        reprojectShp(cl_file, str(pathlib.Path(
            os.path.join(output_folder, root_name+'_lines.shp'))), source_epsg, target_epsg)
        exportToGeojson(str(pathlib.Path(
            os.path.join(output_folder, root_name+'_lines.shp'))))


def writeHtml(image_url_list, titles, output_search_folder, product_name):
    '''
    Description:
    ------------
    HTML file with the quicklook images found in the request are created. This file
    will be launched by the default web browser.

    Arguments:
    ------------
    - image_url_list(list of strings): list of the quicklooks
    - titles (list of strings): titles for each quicklook
    - output_search_folder (string): path to the output html file
    - product_name (string): name of the product for Landsat (to avoid name 
      differences on quicklooks).

    Returns:
    ------------
    html_file (object): html file

    '''

    html = "<html><head>"
    html += "<style> table, th, td {border: 1px solid black; border-collapse: collapse;}</style>"
    html += "<title>Results</title></head><body>"
    html += '<center><table>'
    columns = 2
    rows = int(len(image_url_list)/columns)+1
    for r in range(1, rows+1):
        html += '<tr>'
        for c in range(1, columns+1):
            i = (r-1)*columns+c
            if i <= len(image_url_list):
                image_url = image_url_list[i-1]
                if product_name == 'landsat_ot_c2_l2':
                    image_url = image_url.replace('L2SP', 'L1TP')
                html += f'<td><center><img src="{image_url}" align="center" width="100%" height="100%"></center></td>'
            else:
                html += f'<td></td>'
        html += '</tr>'
        html += '<tr>'
        for c in range(1, columns+1):
            i = (r-1)*columns+c
            if i <= len(image_url_list):
                title = titles[i-1]
                html += f'<td><center><p><font color="blue">[{i-1}] {title}</font></p></center></td>'
            else:
                html += f'<td></td>'
        html += '</tr>'
    html += '</table></center>'
    html += "</body></html>"
    html_file = str(pathlib.Path(os.path.join(
                    os.getcwd(), output_search_folder, 'search_result.html')))
    with open(html_file, 'w') as outputfile:
        outputfile.write(html)
    return html_file


def getS2QuicklookUrlExternalServer(titles):
    # title
    # [1] S2B_MSIL1C_20211203T104319_N0301_R008_T30SYJ_20211203T113738 CC:0.07 Days: -12 Offline
    image_url_list = []
    for title in titles:
        parts = title.split('_')
        level = parts[1]
        date = parts[2][0:8]
        server = 'https://roda.sentinel-hub.com/sentinel-s2-l1c/tiles/'
        tile = parts[5][1:]
        tile = tile[0:2]+'/'+tile[2:3]+'/'+tile[3:]
        date = date[0:4]+'/'+date[4:6]+'/'+date[6:]
        server = server+tile+'/'+date+'/0/preview.jpg'
        image_url_list.append(server)
    return image_url_list


def writeHtmlS2(image_url_list, titles, output_search_folder, product_name, s2_quicklook_server):
    '''
    Description:
    ------------
    HTML file with the quicklook images found in the request are created (for S2 images). This file
    will be launched by the default web browser.

    Arguments:
    ------------
    - image_url_list(list of strings): list of the quicklooks
    - titles (list of strings): titles for each quicklook
    - output_search_folder (string): path to the output html file
    - product_name (string): name of the product for Landsat (to avoid name 
      differences on quicklooks).
    - s2_quicklook_ser (string): quicklook source (Copernicus or external)


    Returns:
    ------------
    html_file (object): html file

    '''

    html = "<html><head>"
    html += "<style> table, th, td {border: 1px solid black; border-collapse: collapse;}</style>"
    html += "<title>Results</title></head><body>"
    html += '<center><table>'
    columns = 2
    rows = int(len(image_url_list)/columns)+1
    for r in range(1, rows+1):
        html += '<tr>'
        for c in range(1, columns+1):
            i = (r-1)*columns+c
            if s2_quicklook_server == '2':
                image_url_list = getS2QuicklookUrlExternalServer(titles)
            if i <= len(image_url_list):
                image_url = image_url_list[i-1]
                if s2_quicklook_server != '2':
                    if product_name == 'landsat_ot_c2_l2':
                        image_url = image_url.replace('L2SP', 'L1TP')
                html += f'<td><center><img src="{image_url}" align="center" width="100%" height="100%"></center></td>'
            else:
                html += f'<td></td>'
        html += '</tr>'
        html += '<tr>'
        for c in range(1, columns+1):
            i = (r-1)*columns+c
            if i <= len(image_url_list):
                title = titles[i-1]
                html += f'<td><center><p><font color="blue">[{i-1}] {title}</font></p></center></td>'
            else:
                html += f'<td></td>'
        html += '</tr>'
    html += '</table></center>'
    html += "</body></html>"
    html_file = str(pathlib.Path(os.path.join(
                    os.getcwd(), output_search_folder, 'search_result.html')))
    with open(html_file, 'w') as outputfile:
        outputfile.write(html)
    return html_file

# *******************************************************************************
# SECTION: FOLDER MANAGEMENT FUNCTIONS
# *******************************************************************************


def createFolderTree(output_data_folder):
    '''
    Description:
    ------------
    Creates folder tree structure for storaging and processing.
    The basic tree folder structure is:
        output_data
         |_ data -> contains downloaded and processed L8 and S2 scenes
         |_ sds -> contains SDS as processing result

    Arguments:
    ------------
    None

    Returns:
    ------------
    - data_path (string):
    - sds_path (string):

    '''

    data_path = str(pathlib.Path(os.path.join(output_data_folder, 'data')))
    createFolderCheck(data_path)
    sds_path = str(pathlib.Path(os.path.join(output_data_folder, 'sds')))
    createFolderCheck(sds_path)
    return data_path, sds_path,


def createFolder(folder_path):
    '''
    Description:
    ------------
    Creates a new folder. If the path already exists, then it will be
    firstly removed.

    Arguments:
    ------------
    - folder_path (string): folder to be created

    Returns:
    ------------
    None

    '''

    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)
    os.makedirs(folder_path)


def createFolderCheck(folder_path):
    '''
    Description:
    ------------
    Creates a new folder only in case this folder does not exist

    Arguments:
    ------------
    - folder_path (string): folder to be created

    Returns:
    ------------
    None

    '''

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def recursiveFileSearch(rootdir='.', pattern='*'):
    '''
    Description:
    ------------
    search for files recursively based in pattern strings

    Arguments:
    ------------
    - rootdir (string): path to the base folder
    - pattern (string): pattern to search files

    Returns:
    ------------
    - matches (list of strings): list of absolute paths to each found file

    '''

    matches = []
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(str(pathlib.Path(os.path.join(root, filename))))
    return matches


def recursiveFolderSearch(rootdir='.', pattern='*'):
    '''
    Description:
    ------------
    search for folders recursively based in pattern strings

    Arguments:
    ------------
    - rootdir (string): path to the base folder
    - pattern (string): pattern to search files

    Returns:
    ------------
    - matches (list of strings): list of absolute paths to each found folder

    '''

    matches = []
    for root, dirnames, filenames in os.walk(rootdir):
        for dirname in fnmatch.filter(dirnames, pattern):
            matches.append(str(pathlib.Path(os.path.join(root, dirname))))
    if len(matches) > 0:
        return matches[0]
    else:
        return ''


def exportLandsatResultsToJson(txt_string):
    '''
    Description:
    ------------
    Exports results for Landsat Searching to JSON file

    Arguments:
    ------------
    - txt_string (string): string of results

    Returns:
    ------------
    - results (dictionary): results in json format

    '''
    txt_list = txt_string.split('\n')[:-1]
    scenes = []
    for lin in txt_list:
        if not 'Central' in lin:
            id = lin.split('[')[1].split(']')[0]
            scene = lin.split('Scene:')[1].split('Cloud')[0].strip()
            ccp = lin.split('cover:')[1].split('%')[0].strip()
            days_off = lin.split('%')[1].split('days')[0].strip()
            scene = {'id': id, 'scene': scene, 'cloud_cover_percentage': ccp,
                     'days_shift_from_central_date': days_off}
            scenes.append(scene)
        else:
            cd = lin.split('date:')[1].strip()
    results = {'saet_results': {'central_date': cd, 'scenes': scenes}}
    return results


def exportSentinelResultsToJson(txt_string):
    '''
    Description:
    ------------
    Exports results for Sentinel Searching to JSON file

    Arguments:
    ------------
    - txt_string (string): string of results

    Returns:
    ------------
    - results (dictionary): results in json format

    '''
    cd = ''
    txt_list = txt_string.split('\n')[:-1]
    availability = []
    scenes = []
    flag = 0
    for lin in txt_list:
        if len(lin) > 0:
            if lin[0] == 'S':
                id = lin.split('Scene:')[1].split('Cloud')[0].strip()
                ccp = lin.split('coverage:')[1].split(
                    'availability:')[0].strip()
                available = lin.split('availability:')[1].strip()
                availability.append(
                    {'scene': id, 'cloud_cover_percentage': ccp, 'available': available})
            if lin[0] == '[':
                if not 'Central' in lin:
                    id = lin.split('[')[1].split(']')[0]
                    scene_id = lin.split('Scene:')[1].split('Cloud')[0].strip()
                    ccp = lin.split('coverage:')[1].split(
                        'availability:')[0].strip()
                    days_off = lin.split('%')[1].split('days')[0].strip()
                    scene = {'id': id, 'scene': scene_id, 'cloud_cover_percentage': ccp,
                             'days_shift_from_central_date': days_off}
                    scenes.append(scene)
                else:
                    cd = lin.split('date:')[1].strip()

    results = {'saet_results': {'central_date': cd,
                                'availability': availability, 'scenes': scenes}}
    return results

# *******************************************************************************
# SECTION: TIME FUNCTIONS
# *******************************************************************************


def startClock():
    '''
    Description:
    ------------
    Inits clock time. It works along with endClock() function

    Arguments:
    ------------
    None

    Returns:
    ------------
    None

    '''

    global _start_time
    _start_time = time.time()


def endClock():
    '''
    Description:
    ------------
    Ends clock time. It works along with startClock() function

    Arguments:
    ------------
    None

    Returns:
    ------------
    None

    '''

    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    message = 'Time passed: {}hour:{}min:{}sec'.format(t_hour, t_min, t_sec)
    # print(message)
    logging.info('Time passed: {}hour:{}min:{}sec'.format(
        t_hour, t_min, t_sec)+'\n')


def getStartEndDate(central_date, offset):
    '''
    Description:
    ------------
    Returns start date and end date from a central date acoording an offset in
    days

    Arguments:
    ------------
    - central_date (string): date in 'YYYYmmdd' format
    - offset (int): days

    Returns:
    ------------
    - start_date (string): start date in 'YYYYmmdd' format
    - end_date (string): end date in 'YYYYmmdd' format

    '''
    central_date = datetime.datetime.strptime(central_date, '%Y%m%d')
    start_date = datetime.datetime.strftime(
        central_date-datetime.timedelta(days=offset), '%Y%m%d')
    end_date = datetime.datetime.strftime(
        central_date+datetime.timedelta(days=offset), '%Y%m%d')
    return start_date, end_date


def valid_date(s):
    '''
    Description:
    ------------
    Validate date argument for central date (--cd) to ensure
    correct format ('YYYYmmdd')

    Arguments:
    ------------
    - central_date (string): date supposed to be validated

    Returns:
    ------------
    - If test passes: valid date in Int type.
    - If test fails: prints error message and exit

    '''
    try:
        if len(s) != 8:
            print(f'Not a valid date: {s}')
            sys.exit()
        d = datetime.datetime.strptime(s, "%Y%m%d")
        return int(datetime.datetime.strftime(d, "%Y%m%d"))
    except ValueError:
        print(f'Not a valid date: {s}')
        sys.exit()


# *******************************************************************************
# SECTION: SHORELINE EXTRACTION FUNCTIONS
# *******************************************************************************

def extractPoints(source_path, pl_path, processing_path, kernel, ppp, degree):
    '''
    Description:
    ------------
    Extract subpixel shoreline points based on kernel analysis over swir1 band, taking as
    template the binary mask of rough shoreline pixel line.
    Values standar used in most of the previous studies with good results are:
    3 (kernel), 4 (ppp), 3 (degree).

    for more information about this algorithm and some results:

    - "Automatic extraction of shorelines from Landsat TM and ETM+ multi-temporal images with subpixel precision". 2012.
      Remote Sensing of Environment. Josep E.Pardo-Pascual, Jaime Almonacid-Caballer, Luis A.Ruiz, Jesus Palomar-Vazquez.

    - "Assessing the Accuracy of Automatically Extracted Shorelines on Microtidal Beaches from Landsat 7,
    Landsat 8 and Sentinel-2 Imagery". 2018. Remote Sensing. Josep E. Pardo-Pascual, Elena Sanchez-Garcia, Jaime Almonacid-Caballer
    Jesus Palomar-Vazquez, Enrique Priego de los Santos, Alfonso Fernández-Sarría, Angel Balaguer-Beser.


    Arguments:
    ------------
    - source_path (string): path to the swir1 band
    - pl_path (string): path to the binary mask of rough shoreline pixel line.
    - processing_path (string): path to the folder to storage results (for each scene,
      this folder is named "temp".
    - kernel (int): kernel size in pixels. Must be an odd number
    - ppp (int): points per pixel. Number of points per pixel extracted. 4 points
      in a 20 m size resolution image means 1 point every 5 meters.
    - degree (int): degree for the mathematical fitting function. Standard values
      are 3 or 5.


    Returns:
    ------------
    - True or False (extraction was success or not)

    '''

    # opens swir1 image
    source_ds = gdal.Open(source_path)
    source_band = source_ds.GetRasterBand(1)
    source_data = source_band.ReadAsArray()

    # opens pixel line mask image
    pl_ds = gdal.Open(pl_path)
    pl_band = pl_ds.GetRasterBand(1)
    pl_data = pl_band.ReadAsArray()

    # creates output text file for coordinate points
    base_name = os.path.basename(source_path).split('.')[0]
    source_data[source_data == float('-inf')] = 0
    if os.path.isfile(str(pathlib.Path(os.path.join(processing_path, base_name+'.d')))):
        os.remove(str(pathlib.Path(os.path.join(processing_path, base_name+'.d'))))
    file_coord = open(str(pathlib.Path(os.path.join(
        processing_path, base_name+'.d'))), 'a')

    # gets swir1 features
    geoTrans = source_ds.GetGeoTransform()
    minXimage = geoTrans[0]
    maxYimage = geoTrans[3]
    dim = source_data.shape
    rows = dim[0]
    columns = dim[1]

    offset = 10  # number of rows and columns preserved to avoid overlapping in adjacent scenes
    c1 = f1 = offset
    c2 = columns - offset
    f2 = rows - offset
    resol_orig = geoTrans[1]  # pixel size
    resol = float(geoTrans[1])/ppp  # point distance
    gap = int(kernel/2)
    points_x = []
    points_y = []
    wm = computeWeights(kernel, ppp)  # weights matrix
    white_pixel = False
    for f in tqdm(range(f1, f2)):
        for c in range(c1, c2):
            valor = pl_data[f, c]
            if valor == 255:  # pixel belongs to the rough pixel line
                white_pixel = True
                nf = f
                nc = c
                # sub-matrix based on kernel size
                sub = source_data[nf-gap:nf+kernel-gap, nc-gap:nc+kernel-gap]
                # sub-matrix resampling based on ppp value
                sub_res = rescale(sub, scale=ppp, order=3, mode='edge')
                cx, cy, cz = createData(sub_res, resol)  # resampled data
                m = polyfit2d(cx, cy, cz, deg=degree)  # fitting data

                # computes laplacian function
                dx = deriva(m, 'x')
                d2x = deriva(dx, 'x')
                dy = deriva(m, 'y')
                d2y = deriva(dy, 'y')
                laplaciano = d2x+d2y

                # get contour points for laplacian = 0
                v = verticeslaplaciano(cx, cy, laplaciano, kernel, ppp)

                if v != None:
                    if len(v) != 0:
                        # if there are more than one contour, we select the contour with highest slope and more centered
                        indice = mejor_curva_pendiente3(v, dx, dy, wm, resol)
                        if indice != -1:
                            linea = v[indice]
                            for i in range(0, len(linea)):
                                par = linea[i]
                                points_x.append(par[0])
                                points_y.append(par[1])

                            # writes the contour points to the text file
                            escribeCoords(points_x, points_y, minXimage, maxYimage,
                                          resol_orig, resol, wm, nf, nc, kernel, file_coord)
                            points_x = []
                            points_y = []
    file_coord.close()

    if white_pixel:
        # variable release
        del source_ds
        del source_band
        del source_data
        del pl_data
        del sub
        del sub_res
        del m
        del cx
        del cy
        del cz
        del wm
        gc.collect()
        return True

    return False


def computeWeights(kernel, ppp):
    '''
    Description:
    ------------
    Computes a matrix with values that follows a normal distribution.
    It is used to ponderate the extracted points based on the distance
    of each point to the center of the image kernel

    Arguments:
    ------------
    - kernel (int): kernel size in pixels. Must be an odd number
    - ppp (int): points per pixel. Number of points per pixel extracted.

    Returns:
    ------------
    - p (numpy matrix): weights matrix

    '''

    p = np.zeros((kernel*ppp, kernel*ppp))
    f, c = p.shape
    cont_i = cont_j = 1.0
    for i in range(0, f):
        for j in range(0, c):
            d = np.sqrt((cont_i-(float(f)+1.0)/2.0)**2 +
                        (cont_j-(float(c)+1.0)/2.0)**2)
            p[i, j] = normcdf(-d, 0, 3)*2
            cont_j += 1
        cont_i += 1
        cont_j = 1
    return p


def normcdf(x, mu, sigma):
    '''
    Description:
    ------------
    Computes the normal distribution value

    Arguments:
    ------------
    - x (float): distance from the center of kernel
    - mu: mean of the normal distribution
    - sigma: standar deviation of the normal distribution

    Returns:
    ------------
    - y (float): normal distribution value

    '''

    t = x-mu
    y = 0.5*erfcc(-t/(sigma*sqrt(2.0)))
    if y > 1.0:
        y = 1.0
    return y


def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196 +
                                                 t*(.09678418+t*(-.18628806+t*(.27886807 +
                                                                               t*(-1.13520398+t*(1.48851587+t*(-.82215223 +
                                                                                                               t*.17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r


def createData(image, resol):
    '''
    Description:
    ------------
    Creates x, y, z arrays of the resampled kernel

    Arguments:
    ------------
    - image (numpy matrix): resampled kernel
    - resol (float): swir1 spatial resolution

    Returns:
    ------------
    - z, y, z (float arrays)

    '''
    inicio = resol-(resol/2.0)  # pixel center
    z = (np.ravel(image)).astype(float)
    tamdata = int(np.sqrt(len(z)))
    x, y = np.meshgrid(np.arange(inicio, tamdata*resol, resol),
                       np.arange(inicio, tamdata*resol, resol))
    x = (np.ravel(x)).astype(float)
    y = (np.ravel(y)).astype(float)
    return x, y, z


def deriva(m, axis):
    '''
    Description:
    ------------
    Computes derivative function of a matrix in a particular axis

    Arguments:
    ------------
    - m (numpy matrix): input matrix
    - axis (string): axis to compute the derivative function

    Returns:
    ------------
    - nm (numpy array): derivative function

    '''

    f, c = m.shape
    if axis == 'x':
        factores = range(1, c)
        nm = m[:, range(1, c)]
        nm = nm*factores
        ceros = np.zeros((f,), dtype=np.float)
        nm = np.vstack((nm.T, ceros.T)).T
        return nm

    if axis == 'y':
        factores = range(1, f)
        nm = m[range(1, f), :]
        nm = (nm.T*factores).T
        ceros = np.zeros((c,), dtype=np.float)
        nm = np.vstack((nm, ceros))
        return nm


def verticeslaplaciano(x, y, m, kernel, ppp):
    '''
    Description:
    ------------
    Computes contour points from laplacian function = 0.
    Uses matplotlib contour function

    Arguments:
    ------------
    - x, y, m (numpy 1D arrays): x, y , z coordinates for laplacian function
    - axis (string): axis to compute the derivative function
    - kernel (int): kernel size in pixels. Must be an odd number
    - ppp (int): points per pixel. Number of points per pixel extracted.

    Returns:
    ------------
    - v (list): list of contour vertices

    '''
    clf()
    v = []
    zz = polyval2d(x, y, m)
    x = np.reshape(x, (kernel*ppp, kernel*ppp))
    y = np.reshape(y, (kernel*ppp, kernel*ppp))
    zz = np.reshape(zz, (kernel*ppp, kernel*ppp))
    try:  # Prevents errors in contour computing
        CS = contour(x, y, zz, 0, colors='y')
        curvas = get_contour_verts(CS)
        for curva in curvas:
            for parte in curva:
                v.append(parte)
        return v
    except:
        return None


def get_contour_verts(cn):
    '''
    Description:
    ------------
    Extract vertices from a contour

    Arguments:
    ------------
    - cn (object): matplotlib contour object

    Returns:
    ------------
    - contours (list): list of contours vertices

    '''

    contours = []
    # for each contour line
    for cc in cn.collections:
        paths = []
        # for each separate section of the contour line
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            paths.append(np.vstack(xy))
        contours.append(paths)
    return contours


def mejor_curva_pendiente3(v, dx, dy, mp, resol):
    '''
    Description:
    ------------
    Select best contour based on the highest mean slope and centrality criteria

    Arguments:
    ------------
    - v (list): list of contours
    - dx (numpy matrix): first derivative of the fitting function in X axis (slope)
    - dy (numpy matrix): first derivative of the fitting function in Y axis (slope)
    - mp (numpy matrix): weight matrix (centrality criteria)

    Returns:
    ------------
    - candidate (int): index of the selected contour

    '''
    pendientes = []
    p_max = 0
    candidate = -1
    for i, curva in enumerate(v):
        for par in curva:
            x = par[0]
            y = par[1]
            px = abs(polyval2d([x], [y], dx))
            py = abs(polyval2d([x], [y], dy))
            p = np.sqrt(px**2+py**2)
            peso = mp[int(x/resol), int(y/resol)]
            pendientes.append(p*peso)
        p_med = np.average(pendientes)
        if p_med >= p_max:
            p_max = p_med
            candidate = i
        pendientes = []
    return candidate


def escribeCoords(x, y, xmin, ymax, resol_orig, resol, wm, fil, col, kernel, output_file):
    '''
    Description:
    ------------
    Write extracted contours vertices coordinates to the .txt file.
    The original points are in subpixel image coordinates. They have to be
    converted to world coordinates

    Arguments:
    ------------
    - x (list): list of X coordinates
    - y (list): list of Y coordinates
    - xmin (float): minimum X coordinate of the swir1 image
    - ymin (float): minimum Y coordinate of the swir1 image
    - resol_orig (float): spatial resolution of the swir1 image
    - resol (float): map distance among each extracted point (pixel size / ppp)
    - wm (numpy matrix): weight matrix
    - fil (int): row coordinate for the center pixel in the current kernel
    - col (int): column coordinate for the center pixel in the current kernel
    - kernel (int): kernel size in pixels
    - output_file (string): path to the output file

    Returns:
    ------------
    None

    '''

    for i in range(0, len(x)):
        # coordenadas punto sobre imagen global
        rx = xmin+(col-int(kernel/2.0))*resol_orig+x[i]
        ry = ymax-(fil-int(kernel/2.0))*resol_orig-y[i]
        peso = wm[int(x[i]/resol), int(y[i]/resol)]
        output_file.write(str(rx)+","+str(ry)+","+str(peso)+'\n')


# *******************************************************************************
# SECTION: AVERAGE POINTS FUNCTIONS
# *******************************************************************************

def averagePoints(source_path, processing_path, cluster_distance, min_cluster_size):
    '''
    Description:
    ------------
    Takes the subpixel rough extracted points and computes the average points.
    The algorithm scan in X and Y direction points with the same coordinates and
    makes groups of points by using clustering criteria as maximum distance among
    points and minimum number of points in a cluster.
    Creates a .txt file with the results

    Arguments:
    ------------
    - source_path (string): path to the scene
    - processing_path (string): path to the processing folder
    - cluster_distance (int): maximum distance between two points to be consider a cluster
    - min_cluster_size (int): minimum number of points in a cluster

    Returns:
    ------------
    None

    '''

    # reading coordinates from extracted subpixel points in the kernel analysis
    base_name = os.path.basename(source_path).split('.')[0]
    file_name = str(pathlib.Path(
        os.path.join(processing_path, base_name+'.d')))
    with open(file_name, 'r') as fichero:
        iter1 = csv.reader(fichero, delimiter=',')
        datos = np.asarray([[dato[0], dato[1]]
                           for dato in iter1]).astype(float)
        fichero.seek(0)
        iter2 = csv.reader(fichero, delimiter=',')
        pesos = np.asarray([dato[2] for dato in iter2]).astype(float)

    ejex = np.unique(datos[:, 0])  # unique values on the x-axis
    ejey = np.unique(datos[:, 1])  # unique values on the y-axis

    # computing clusters
    medias = creaCluster(datos, pesos, ejex, ejey,
                         cluster_distance, min_cluster_size)

    # writes results to the output file (average x and average y of every cluster)
    with open(str(pathlib.Path(os.path.join(processing_path, base_name+'.m'))), 'w') as fichero:
        for media in medias:
            fichero.write(str(media[0])+","+str(media[1])+"\n")


def creaCluster(d, p, ex, ey, cluster_distance, min_cluster_size):
    '''
    Description:
    ------------
    Makes groups of points acording clustering criteria (maximum distance among
    points and minimum number of points in a cluster). From each group, the algorithm
    computes a ponderate average value for x and y coordinates based on a weight matrix.

    Arguments:
    ------------
    - d (numpy array): list of X-Y coordinates
    - p (numpy matrix): weight matrix
    - ex (numpy array): unique values on the x-axis
    - ey (numpy array): unique values on the y-axis
    - cluster_distance (int): maximum distance between two points to be consider a cluster
    - min_cluster_size (int): minimum number of points in a cluster

    Returns:
    ------------
    - average_points (list): list of average points for each cluster

    '''

    tol = cluster_distance
    average_points = []

    # clustering in x-axis
    for x in ex:
        id_x = np.nonzero(d[:, 0] == x)
        cy = d[:, 1][id_x]
        pey = p[id_x]
        if len(cy) >= 2:
            orig_coord, pos = getClusters(cy, tol)
            for cp in pos:
                cluster = orig_coord[cp]
                if len(cluster) >= min_cluster_size:
                    p_cluster = pey[cp]
                    media_y = np.average(cluster, weights=p_cluster)
                    average_points.append([x, media_y])

    # clustering in y-axis
    for y in ey:
        id_y = np.nonzero(d[:, 1] == y)
        cx = d[:, 0][id_y]
        pex = p[id_y]
        if len(cx) >= 2:
            orig_coord, pos = getClusters(cx, tol)
            for cp in pos:
                cluster = orig_coord[cp]
                if len(cluster) >= min_cluster_size:
                    p_cluster = pex[cp]
                    media_x = np.average(cluster, weights=p_cluster)
                    average_points.append([media_x, y])
    return average_points


def getClusters(coord, tol):
    '''
    Description:
    ------------
    Makes groups of points based on a maximum distance.

    Arguments:
    ------------
    - coord (list): list of point coordinates with the same x or y value
    - tol (int): cluster distance (maximum distance between two points to
      be consider a cluster)

    Returns:
    ------------
    - orig_coord (list): list of point coordinates with the same x or y value
    - pos (list): index of points that belong tho the same cluster

    '''

    clusters = []
    cluster = []
    orig_coord = coord.copy()
    coord.sort()
    cluster.append(0)
    for i in range(0, len(coord)-1):
        current = coord[i]
        siguiente = coord[i+1]
        dist = siguiente-current
        if dist <= tol:
            cluster.append(i+1)
        else:
            clusters.append(cluster)
            cluster = []
            cluster.append(i+1)
    clusters.append(cluster)
    parcial = []
    pos = []
    for c in clusters:
        for iden in c:
            a, = np.where(orig_coord == coord[iden])
            parcial.append(a[0])
        pos.append(parcial)
        parcial = []
    return orig_coord, pos


# *******************************************************************************
# SECTION: SHAPEFILE CREATION FUNCTIONS
# *******************************************************************************

def createShpFromAverageFile(source_path, processing_path):
    '''
    Description:
    ------------
    Converts average point coordinates from .txt file to .shp file.
    The name of the shapefile is based on the name of the .xtx file.
    In order to copy the beach code to each average point, an attribute
    field call "BEACH_CODE" is added.

    Arguments:
    ------------
    - source_path (string): path to the template swir1 image
    - processing_path (string): path to the processing output path.

    Returns:
    ------------
    - shp_path (string): path to the .shp file

    '''

    # gets projection from template image
    source_ds = gdal.Open(source_path)
    prj = osr.SpatialReference()
    prj.ImportFromWkt(source_ds.GetProjectionRef())

    # reads average point coordinates from .txt file
    base_name = os.path.basename(source_path).split('.')[0]
    file_name = str(pathlib.Path(
        os.path.join(processing_path, base_name+'.m')))
    with open(file_name, 'r') as f:
        data = csv.reader(f, delimiter=',')
        coords = np.asarray([[dato[0], dato[1]]
                            for dato in data]).astype(float)

    # creates a new shapefile and adds a database structure
    driver = ogr.GetDriverByName("ESRI Shapefile")
    shp_path = str(pathlib.Path(os.path.join(
        processing_path, base_name+'.shp')))
    if os.path.exists(shp_path):
        driver.DeleteDataSource(shp_path)
    data_source = driver.CreateDataSource(shp_path)
    layer = data_source.CreateLayer(' ', geom_type=ogr.wkbPoint, srs=prj)
    id_field = ogr.FieldDefn("Id_pnt", ogr.OFTInteger)
    id_field.SetWidth(10)
    layer.CreateField(id_field)
    id_field2 = ogr.FieldDefn("BEACH_CODE", ogr.OFTInteger)
    id_field2.SetWidth(10)
    layer.CreateField(id_field2)
    layer_defn = layer.GetLayerDefn()

    # creates shapefile features
    for i in range(0, len(coords)):
        values = coords[i]
        feat = ogr.Feature(layer_defn)
        pnt = ogr.Geometry(ogr.wkbPoint)
        pnt.AddPoint(values[0], values[1])
        feat.SetGeometry(pnt)
        feat.SetField('Id_pnt', i)
        feat.SetField('BEACH_CODE', 0)
        layer.CreateFeature(feat)
    data_source.FlushCache()
    source_ds = None
    return shp_path


def copyShpIdentifiers(shp_polygons, shp_points):
    '''
    Description:
    ------------
    Copies the beach code from the beach shapefile to each average point using
    two geometry intersection test.

    Arguments:
    ------------
    - shp_polygons (string): path to beaches shapefile
    - shp_points (string): path to the points shapefile.

    Returns:
    ------------
    None

    '''

    # opens both polygons and points shapefiles
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds_pol = driver.Open(shp_polygons, 0)
    ds_point = driver.Open(shp_points, 1)

    layer_pol = ds_pol.GetLayer()
    layer_point = ds_point.GetLayer()

    layer_point_defn = layer_point.GetLayerDefn()
    field_names = [layer_point_defn.GetFieldDefn(
        i).GetName() for i in range(layer_point_defn.GetFieldCount())]

    # adds beach code field
    if not 'BEACH_CODE' in field_names:  # to ensure that attribute "BEACH_CODE" exists
        new_field = ogr.FieldDefn('BEACH_CODE', ogr.OFTInteger)
        layer_point.CreateField(new_field)

    # populates beach code field
    for feat_pol in layer_pol:
        id_feat_pol = feat_pol.GetField('BEACH_CODE')
        geom_pol = feat_pol.GetGeometryRef()
        geom_envelope = envelopeToGeom(geom_pol)
        for feat_point in layer_point:
            geom_point = feat_point.GetGeometryRef()
            if geom_point.Intersect(geom_envelope):  # first intersection test
                if geom_point.Intersect(geom_pol):  # second intersection test
                    feat_point.SetField('BEACH_CODE', id_feat_pol)
                    layer_point.SetFeature(feat_point)

    ds_point.FlushCache()
    ds_point = None
    ds_pol = None


def envelopeToGeom(geom):
    '''
    Description:
    ------------
    Returns the bounding box of a polygon geometry.

    Arguments:
    ------------
    - geom (objetc): ogr geometry object

    Returns:
    ------------
    poly_envelope (object): ogr geometry object

    '''

    # gets bounding box from geometry
    (minX, maxX, minY, maxY) = geom.GetEnvelope()

    # creates ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(minX, minY)
    ring.AddPoint(maxX, minY)
    ring.AddPoint(maxX, maxY)
    ring.AddPoint(minX, maxY)
    ring.AddPoint(minX, minY)

    # creates polygon
    poly_envelope = ogr.Geometry(ogr.wkbPolygon)
    poly_envelope.AddGeometry(ring)
    return poly_envelope


# *******************************************************************************
# SECTION: POINT CLEANING FUNCTIONS
# *******************************************************************************

def cleanPoints2(shp_path, tol_rba, level):
    '''
    Description:
    ------------
    Remove outliers points based on two criteria:
    - longest spanning tree algorithm (LST).
    - angle tolerance.

    To improve the performance, the algorithm uses an initial Delaunay triangulation
    to create a direct graph.

    Two versions of cleaned points shapefile is created: point and line versions

    More information:
    "An efficient protocol for accurate and massive shoreline definition from
    mid-resolution satellite imagery". 2020. Coastal Engineering. E. Sanchez-García,
    J.M. Palomar-Vazquez, J.E. Pardo-Pascual, J. Almonacid-Caballer, C. Cabezas-Rabadan,
    L. Gomez-Pujol.

    Arguments:
    ------------
    - shp_path (string): path to the points shapefile
    - tol_rba (int): angle tolerance
    - level (int): takes one point every n (level) points. Speeds the process

    Returns:
    ------------
    None

    '''
    # opens the shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    source_ds = driver.Open(shp_path, 0)
    prj = source_ds.GetLayer().GetSpatialRef()
    base_name = os.path.basename(shp_path).split('.')[0]
    dir_name = os.path.dirname(shp_path)
    layer = source_ds.GetLayer()

    # gest list of unique BEACH_CODE values
    ids = []
    for feature in layer:
        id_feat = feature.GetField('BEACH_CODE')
        if not id_feat is None:
            ids.append(id_feat)
    ids = list(set(ids))
    layer.ResetReading()

    ids.sort()
    # prevents from points with BEACH_CODE 0 (outside of any beach area).
    if ids[0] == 0:
        ids.remove(0)  # removes points

    # creates groups of points with the same BEACH_CODE value
    groups = []
    identifiers = []
    for id_feat in ids:
        geometries = []
        layer.SetAttributeFilter("BEACH_CODE = "+str(id_feat))
        for feature in layer:
            geom = feature.GetGeometryRef()
            if not geom is None:
                geometries.append(geom.Clone())
        groups.append(geometries)
        identifiers.append(id_feat)

    # process each group separately
    clean_geometries = []
    level = 1
    for i in range(0, len(groups)):
        group = groups[i]
        identifier = identifiers[i]
        coords = []
        ng = float(len(group))
        # prevents from too much long numer of points in a group
        # level = ceil(ng/group_size)
        for i in range(0, len(group), level):
            geom = group[i].Clone()
            coords.append([geom.GetX(), geom.GetY()])
        points = np.array(coords)
        if len(points >= 4):  # delaunay triangulation needs 4 or more points
            try:
                tri = Delaunay(points)
                # list of triangles
                lista_tri = tri.simplices
                # list of ids of the connected points wiht LST
                lst = computeLST(lista_tri, points)
                # remove LST point by angle tolerance
                clean_points = cleanPointsByAngle(lst, points, tol_rba)
                # list of cleaned points with its identifier
                clean_geometries.append([clean_points, identifier])
            except:
                pass

    # crates point and line versions of the cleaned points
    makePointShp(str(pathlib.Path(os.path.join(dir_name, base_name +
                 '_cp.shp'))), clean_geometries, prj)
    makeLineShp(str(pathlib.Path(os.path.join(dir_name, base_name +
                '_cl.shp'))), clean_geometries, prj)
    source_ds = None


def computeLST(lista_tri, puntos):
    '''
    Description:
    ------------
    Computes Longest Spanning Tree (LST) of a directed graph.

    Arguments:
    ------------
    - lista_tri (list): list of triangles from Delaunay triangulation
    - puntos (list): list of original points (unclean)

    Returns:
    ------------
    final (list): list of index of points belonging to the LST

    '''

    # creates a graph from triangles
    G = nx.Graph()
    for tri in lista_tri:
        G.add_edge(tri[0], tri[1], weight=dist(puntos[tri[0]], puntos[tri[1]]))
        G.add_edge(tri[1], tri[2], weight=dist(puntos[tri[1]], puntos[tri[2]]))
        G.add_edge(tri[2], tri[0], weight=dist(puntos[tri[2]], puntos[tri[0]]))
    # computes Minimum Spanning tree (MST)
    MST = nx.minimum_spanning_tree(G)
    MST_directed_graph = nx.Graph(MST).to_directed()
    # select the first end point of the grahp from and arbitrary point
    lpath1 = nx.single_source_bellman_ford_path_length(MST_directed_graph, 1)
    indice1 = list(lpath1.keys())[-1]
    # select the second end point of the graph from the fisrt end point
    lpath2 = nx.single_source_bellman_ford_path_length(
        MST_directed_graph, indice1)
    indice2 = list(lpath2.keys())[-1]
    # computes the MST between the first and second end points (LST)
    final = nx.dijkstra_path(MST_directed_graph, indice1, indice2)
    return final


def dist(p1, p2):
    '''
    Description:
    ------------
    Returns distance between two couple of coordinates.

    Arguments:
    ------------
    - p1 (list): first couple of coordinates
    - p2 (list): second couple of coordinates

    Returns:
    ------------
    - (float): distance

    '''
    return sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)


def cleanPointsByAngle(lst, puntos, tol_rba):
    '''
    Description:
    ------------
    Removes a point if the angle formed between the anterior and posterior
    point is greater than a tolerance.

    Arguments:
    ------------
    - lst (list): list of index in LST graph
    - puntos (list): list of LST points
    - tol_rba (int): angle tolerance

    Returns:
    ------------
    - clean_points (list): list of point ogr geometries

    '''

    clean_points = []
    lineas = []
    for iden in lst:
        lineas.append(puntos[iden])
    if tol_rba != None:
        lineas = removeByAngle(lineas, tol_rba)

    for i in range(0, len(lineas)):
        new_geom = ogr.Geometry(ogr.wkbPoint)
        pt = lineas[i]
        new_geom.AddPoint(pt[0], pt[1])
        clean_points.append(new_geom.Clone())
    return clean_points


def removeByAngle(puntos, tol):
    '''
    Description:
    ------------
    Removes a point if the angle formed between the anterior and posterior
    point is greater than a tolerance. Function included in cleanPointsByAngle()

    Arguments:
    ------------
    - puntos (list): list of LST points
    - tol (int): angle tolerance

    Returns:
    ------------
    - final_points (list): list of cleaned points

    '''

    final_points = []
    for i in range(0, len(puntos)-2):
        a = puntos[i]
        b = puntos[i+1]
        c = puntos[i+2]
        ang = degrees(atan2(c[1]-b[1], c[0]-b[0]) -
                      atan2(a[1]-b[1], a[0]-b[0]))
        if ang < 0:
            ang = ang + 360
        if ang >= tol:
            final_points.append(puntos[i+1])
    final_points.insert(0, puntos[0])
    final_points.append(puntos[len(puntos)-1])
    return final_points


def makePointShp(shp_path, source_list, prj):
    '''
    Description:
    ------------
    Creates a point shapefile from a list of points.

    Arguments:
    ------------
    - shp_path (string): path to the output shapefile
    - source_list (list of list): list of points
    - prj (object): ogr coordinate spatial reference object

    Returns:
    ------------
    None

    '''
    # creates empty point shapefile
    data_source, layer = createEmptyShp(shp_path, 'point', prj)
    layer_defn = layer.GetLayerDefn()

    # adds points to new shapefile along with the BEACH_CODE attribute
    for i in range(0, len(source_list)):
        points = source_list[i][0]
        identifier = source_list[i][1]
        for pt in points:
            feat = ogr.Feature(layer_defn)
            feat.SetGeometry(pt)
            feat.SetField('BEACH_CODE', identifier)
            layer.CreateFeature(feat)
    data_source.FlushCache()
    data_source = None


def makeLineShp(shp_path, source_list, prj):
    '''
    Description:
    ------------
    Creates a polyline shapefile from a list of points.

    Arguments:
    ------------
    - shp_path (string): path to the output shapefile
    - source_list (list of list): list of points
    - prj (object): ogr coordinate spatial reference object

    Returns:
    ------------
    None

    '''
    # creates empty point shapefile
    data_source, layer = createEmptyShp(shp_path, 'line', prj)
    layer_defn = layer.GetLayerDefn()

    # crates lines from  each point sublist to new shapefile along with
    # the BEACH_CODE attribute
    for i in range(0, len(source_list)):
        points = source_list[i][0]
        identifier = source_list[i][1]
        feat = ogr.Feature(layer_defn)
        new_geom = ogr.Geometry(ogr.wkbLineString)
        for point_geom in points:
            new_geom.AddPoint(point_geom.GetX(), point_geom.GetY())
        feat.SetGeometry(new_geom)
        feat.SetField('BEACH_CODE', identifier)
        layer.CreateFeature(feat)
    data_source.FlushCache()
    data_source = None


def createEmptyShp(shp_path, geom_type, prj):
    '''
    Description:
    ------------
    Creates an empty shapefile of a specific geometry and projection.
    Adds the field 'BEACH_CODE'.

    Arguments:
    ------------
    - shp_path (string): path to the output shapefile
    - geom_type (string): type of geometry ('point', 'line' or 'polygon')
    - prj (object): ogr coordinate spatial reference object

    Returns:
    ------------
    - data_source (object): ogr data source object
    - layer (object): ogr layer object

    '''

    dict_geom = {'point': ogr.wkbPoint,
                 'line': ogr.wkbLineString, 'polygon': ogr.wkbPolygon}
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(shp_path):
        driver.DeleteDataSource(shp_path)
    data_source = driver.CreateDataSource(shp_path)
    layer = data_source.CreateLayer(
        ' ', geom_type=dict_geom[geom_type], srs=prj)
    id_field = ogr.FieldDefn("BEACH_CODE", ogr.OFTInteger)
    id_field.SetWidth(10)
    layer.CreateField(id_field)

    return data_source, layer
