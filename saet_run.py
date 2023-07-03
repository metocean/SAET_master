'''
SAET (Shoreline Analysis and Extraction Tool). V2.0
Author: Geo-Environmental Cartography and Remote Sensing Group (CGAT).
Universitat Politecnica de Valencia (Spain). October 2022.
https://cgat.webs.upv.es/

This software is the intellectual property of the UPV and has been developed
within the framework of the European ECFAS project by the following authors:
Jesus Palomar Vazquez, Jaime Almonacid Caballer, Josep E. Pardo Pascual and
Carlos Cabezas Rabadan.

This is the main module. This module must be runned in command line style.
The parameters, optional or mandatory, and its formats are explained below.
To obtain help about the usage, type: "python saet_run.py -h"

# *******************************************************************************
# PARAMETERS
# *******************************************************************************


  -h, --help    show this help message and exit

  Arguments:
  --rm RM       Run mode (search [os] / download and process [dp] / only download [od] / reprocess downloaded images [rp] / offline S2 retrieval [or]). --rm=os / --rm=dp / --rm=od / --rm=rp / --rm=or. Default: os
  --fp FP       Footprint for searching scenes.
                - path of the roi file for searching scenes (--fp=c:\data\roi.geojson)
                - coordinates long/lat in this format: fp=long_min,lat_min,long_max,lat_max
                - NONE
  --sd SD       Start date for searching scenes (YYYYMMDD). --cd=20210101.
  --cd CD       Central date for storm event (YYYYMMDD). --cd=20210101.
  --ed SD       End date for searching scenes (YYYYMMDD). --cd=20210101.
  --mc [0-100]  maximum cloud coverture for the whole scene [0-100]. --mc=10
  --lp LP       Product type for Landsat 8. (landsat_ot_c2_l1 / 'landsat_ot_c2_l2' / NONE). --lp=landsat_ot_c2_l1. Default: landsat_ot_c2_l1
  --ll LL       Scene list for Landsat 8 (NNNNNN). --ll=198032,198033. Default: NONE
  --sp SP       Product type for Sentinel 2 (S2MSI1C / S2MSI2A / NONE). --sp=S2MSI1C. Default: S2MSI1C
  --bc BC       List of numbers to filter the extraction process by beach codes. --bc=1680,1758. Default: NONE 
  --sl SL       Scene list for Sentinel 2 (NNSSS). --ll=30TYJ,30TYK. Default: NONE
  --of OF       Output folder path. --of=c:\data (windows) | --of=/data (linux)
  --wi WI       Water index type (aweish, aweinsh,mndwi,kmeans). --wi=aweinsh. Default: aweinsh
  --th TH       Thresholding method (0: standard 0 value, 1: Otsu bimodal, 2: Otsu multimodal 3 classes). --th=0. Default: 0
  --mm MM       Morphological method (erosion, dilation). --mm=dilation, Default: dilation
  --cl CL       Cloud mask level (0: no masking, 1: only opaque clouds, 2: opaque clouds + cirrus + cloud shadows). Default: 0
  --ks KS       kernel size for points extraction (3,5). --ks=3. Default: 3
  --np NP       List of number of products for download (only if --rm=od or --rm=dp). [0,2,5,3] / [*] / [5-10]. Default: NONE
  --oa OA       Offline S2 activation (check / activate). Only if --rm=or. Default: check

Examples:
file "examples_of_use.txt"

'''
# imports
import os
import sys
import argparse
from pathlib import Path
import logging
import logging.handlers
import saet_config
from saet_tools import (recursiveFolderSearch,
                        startClock,
                        endClock,
                        startSearchForLandsat8,
                        startSearchForSentinel2,
                        createFolderTree,
                        createFolderCheck,
                        processLandsat8Scenes,
                        processSentinel2Scenes,
                        checkCommandArgumments,
                        valid_date,
                        parseNumberProducts,
                        downloadLandsat8_2,
                        downloadS2_2,
                        filterScenesInfolder,
                        processSentinel2SceneFromPath,
                        processLandsatSceneFromPath,
                        offlineS2Activation
                        )


def path_venv():
    '''
    Determine if SAET is running into a virtual environment (recommended)
    '''
    try:
        return os.environ['VIRTUAL_ENV']
    except:
        return ''


# to avoid problems with proj.db
if sys.platform == 'win32':
    pvenv = path_venv()
    if pvenv == '':
        python_path = Path(os.path.dirname(str(sys.executable)))
        path_proj = python_path / 'Lib/site-packages/osgeo/data/proj'
        os.environ['PROJ_LIB'] = str(path_proj)
    else:
        pvenv = Path(pvenv)
        path_proj = pvenv / 'Lib/site-packages/osgeo/data/proj'
        os.environ['PROJ_LIB'] = str(path_proj)
# if (sys.platform == 'linux') or (sys.platform == 'linux2'):
#    path_proj = '/usr/lib/python3/dist-packages'
#    os.environ['PROJ_LIB'] = str(path_proj)

# init logging
logger = logging.getLogger("")


def init_logger(level):
    '''
    Init the logging process. The level of logging can be changed from configuration file "saet_config.py"
    '''

    if level == '10':
        log_level = logging.DEBUG
    if level == '20':
        log_level = logging.INFO
    if level == '30':
        log_level = logging.WARNING
    if level == '40':
        log_level = logging.ERROR
    if level == '50':
        log_level = logging.CRITICAL

    logger.setLevel(log_level)

    # console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    if (log_level == logging.INFO) or (log_level == logging.WARNING):
        console_formatting = logging.Formatter(
            '%(asctime)s %(levelname)s %(message)s')
    else:
        console_formatting = logging.Formatter(
            '%(asctime)s %(filename)s-%(funcName)s %(levelname)s %(message)s')

    console_handler.setFormatter(console_formatting)
    logger.addHandler(console_handler)


def parse_args():
    '''
    Function to validate the input parameters. These parameters are obtained from the command-line
    sentence. 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--rm',
                        type=str,
                        choices=['os', 'dp', 'od', 'op', 'or'],
                        help='Run mode (only search [s] / download and process [dp] / only donwload [od] / only process [op] / offline S2 retrieval [or]). --rm=os / --rm=dp / --rm=od / --rm=op / --rm=or. Default: os',
                        default='os',
                        required=True)
    parser.add_argument('--fp',
                        type=str,
                        help='path of the roi file for searching scenes (fp=c:\data\roi.geojson), coordinates long/lat in this format: fp=long_min,lat_min,long_max,lat_max. Default: NONE',
                        default='NONE',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--sd',
                        help='Start date for searching scenes (YYYYMMDD). --sd=20210101. Default:20200101',
                        type=valid_date,
                        default='20200101',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--cd',
                        help='Central date for storm (YYYYMMDD). --sd=20210101. Default:20200102',
                        type=valid_date,
                        default='20200102',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--ed',
                        help='End date for searching scenes (YYYYMMDD). --sd=20210101. Default:20200103',
                        type=valid_date,
                        default='20200103',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--mc',
                        type=int,
                        choices=range(0, 101),
                        metavar='[0-100]',
                        help='maximum cloud coverture for the whole scene [0-100]. --mc=10',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--lp',
                        type=str,
                        choices=['landsat_ot_c2_l1',
                                 'landsat_ot_c2_l2', 'NONE'],
                        help='Landsat 8 product type. landsat_ot_c2_l1 or landsat_ot_c2_l2 or NONE. Default: landsat_ot_c2_l1',
                        default='landsat_ot_c2_l1',
                        required=not (('--rm=op' in sys.argv) or ('--rm=or' in sys.argv)))
    parser.add_argument('--ll',
                        type=str,
                        help='List of scenes for Landsat 8 (number of 6 digits). --ll=198032,199031. Default: NONE',
                        default='NONE',
                        required=not (('--rm=op' in sys.argv) or ('--rm=or' in sys.argv)))
    parser.add_argument('--sp',
                        type=str,
                        choices=['S2MSI1C', 'S2MSI2A', 'NONE'],
                        help='Sentinel 2 product type (S2MSI1C / S2MSI2A). --s2=S2MSI1C / --s2=S2MSI2A / NONE. Default: S2MSI1C',
                        default='S2MSI1C',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--sl',
                        type=str,
                        help='List of scenes for Sentinel 2 (string of 5 characters). --sl=31TCF,30TYK. Default: NONE',
                        default='NONE',
                        required=not '--rm=op' in sys.argv)
    parser.add_argument('--bc',
                        type=str,
                        help='beach code filter list. --bc=520,548 Default: NONE',
                        default='NONE',
                        required=False)  # ('--rm=dp' in sys.argv) or ('--rm=op') in sys.argv)
    parser.add_argument('--of',
                        type=str,
                        help='output data folder. --of=c:\data (windows) --of=/data. Default: SAET_HOME_PATH',
                        default='SAET_HOME_PATH',
                        required=False)
    parser.add_argument('--wi',
                        type=str,
                        choices=['aweish', 'aweinsh', 'mndwi', 'kmeans'],
                        help='Water index type (aweish, aweinsh,mndwi,kmeans). --wi=aweinsh. Default: aweinsh',
                        default='aweinsh',
                        required=False)
    parser.add_argument('--th',
                        type=str,
                        choices=['0', '1', '2'],
                        help='Thresholding method (0: standard 0 value, 1: Otsu bimodal, 2: Otsu multimodal 3 classes). --th=0. Default: 0',
                        default='0',
                        required=False)
    parser.add_argument('--mm',
                        type=str,
                        choices=['erosion', 'dilation'],
                        help='Morphological method (erosion, dilation). --mm=dilation, Default: dilation',
                        default='dilation',
                        required=False)
    parser.add_argument('--cl',
                        type=str,
                        choices=['0', '1', '2'],
                        help='Cloud mask level (0: no masking, 1: only opaque clouds, 2: opaque clouds + cirrus + cloud shadows). Default: 0',
                        default='0',
                        required=False)
    parser.add_argument('--ks',
                        type=str,
                        choices=['3', '5'],
                        help='Kernel size for points extraction. Default: 3',
                        default='3',
                        required=False)
    parser.add_argument('--np',
                        type=str,
                        help='List of number of products for download (only if --rm=d and --rm=p). [0,2,5,3] / [*] / [5-10]. Default: NONE',
                        default='NONE',
                        required=('--rm=dp' in sys.argv) or ('--rm=od' in sys.argv))
    parser.add_argument('--oa',
                        type=str,
                        choices=['check', 'activate'],
                        help='Offline S2 activation (only if --rm=or). "check" / "activate". Default: "check"',
                        default='check',
                        required=('--rm=or' in sys.argv))

    return parser.parse_args()


# main function
def run_algo(args):

    # check arguments for some parameters
    run_parameters = checkCommandArgumments(args)

    # get some configuration parameters
    # these parameters are configured using saet_config.py file
    user_esa = os.getenv('USER_ESA')
    pass_esa = os.getenv('PASS_ESA')
    user_usgs = os.getenv('USER_USGS')
    pass_usgs = os.getenv('PASS_USGS')
    beaches_path = os.getenv('SHP_BEACHES_PATH')
    l8grid_path = os.getenv('SHP_LANDSAT_GRID_PATH')
    s2grid_path = os.getenv('SHP_SENTINEL2_GRID_PATH')
    quicklook = os.getenv('QUICKLOOK')
    s2_quicklook_server = os.getenv('S2_QUICKLOOK_SERVER')
    output_res = os.getenv('OUT_RES')
    ridap = os.getenv('RIDAP')

    # managing output data folder and search data folder from input parameters
    if args.of == 'SAET_HOME_PATH':
        saet_home_path = str(
            Path(os.path.dirname(os.path.realpath(__file__))))
        output_data_folder = Path(os.path.join(saet_home_path, 'output_data'))
    else:
        custom_data_folder = str(Path(args.of))
        output_data_folder = Path(os.path.join(
            custom_data_folder, 'output_data'))
    output_search_folder = str(
        Path(os.path.join(output_data_folder, 'search_data')))

    # these processing parameters are exposed by command line options
    water_index = args.wi
    thresholding_method = args.th
    morphology_method = args.mm
    cloud_mask_level = args.cl
    kernel_size = args.ks
    offline_activation = args.oa

    # adding credentials to run_parameters
    run_parameters['user_esa'] = user_esa
    run_parameters['pass_esa'] = pass_esa
    run_parameters['user_usgs'] = user_usgs
    run_parameters['pass_usgs'] = pass_usgs

    # auxiliar data paths
    if not os.path.isfile(beaches_path):
        logger.warning(
            f'The file {beaches_path} does not exist'+'\n')
        sys.exit(1)
    if not os.path.isfile(l8grid_path):
        logger.warning(
            f'The file {l8grid_path} does not exist'+'\n')
        sys.exit(1)
    if not os.path.isfile(s2grid_path):
        logger.warning(
            f'The file {s2grid_path} does not exist'+'\n')
        sys.exit(1)

    # adding aux_path shp data to run parameters
    run_parameters['l8grid_path'] = l8grid_path
    run_parameters['s2grid_path'] = s2grid_path

    # adding output results mode to run parameters
    run_parameters['output_res'] = output_res
    run_parameters['output_data_folder'] = output_data_folder
    run_parameters['output_search_folder'] = output_search_folder
    run_parameters['quicklook'] = quicklook
    run_parameters['s2_quicklook_server'] = s2_quicklook_server

    # init clock
    startClock()

    l8_scenes = []
    s2_scenes = []

    # ONLY SEARCH FOR IMAGES _______________________________________________________________
    if run_parameters['run_mode'] == 'os':  # only search mode
        createFolderCheck(output_search_folder)
        if run_parameters['l8_product'] != 'NONE':
            l8_scenes = startSearchForLandsat8(run_parameters)
        if run_parameters['s2_product'] != 'NONE':
            s2_scenes, s2_titles = startSearchForSentinel2(
                run_parameters, len(l8_scenes))

    # OFFLINE S2 RETRIEVAL _________________________________________________________________
    if run_parameters['run_mode'] == 'or':  # only search mode
        if run_parameters['s2_product'] != 'NONE':
            s2_scenes, s2_titles = startSearchForSentinel2(
                run_parameters, len(l8_scenes))
            if offline_activation == 'activate':
                offlineS2Activation(s2_scenes, run_parameters)

    # ONLY DOWNLOAD IMAGES _________________________________________________________________
    if run_parameters['run_mode'] == 'od':  # only donwload mode
        createFolderCheck(output_search_folder)
        if run_parameters['l8_product'] != 'NONE':
            l8_scenes = startSearchForLandsat8(run_parameters)
        if run_parameters['s2_product'] != 'NONE':
            s2_scenes, s2_titles = startSearchForSentinel2(
                run_parameters, len(l8_scenes))

        if len(l8_scenes) != 0 or len(s2_scenes) != 0:
            # making folder structure for data processing
            data_path, sds_path = createFolderTree(output_data_folder)
        if len(l8_scenes) != 0:
            filtered_l8_scenes = parseNumberProducts(args.np, l8_scenes)
            if len(filtered_l8_scenes) != 0:
                print(str(len(filtered_l8_scenes)) +
                      ' Landsat scenes to download'+'\n')
                for scene_l8 in filtered_l8_scenes:
                    print(scene_l8['display_id'])
                    scene_path, scene_id = downloadLandsat8_2(
                        scene_l8, data_path, user_usgs, pass_usgs)
        if len(s2_scenes) != 0:

            # selection of suitable bands to be downloaded
            if run_parameters['s2_product'] == 'S2MSI2A':
                bands = ['B02', 'B03', 'B08', 'B11', 'B12', 'SCL']

            if run_parameters['s2_product'] == 'S2MSI1C':
                bands = ['B02', 'B03', 'B08', 'B11', 'B12', 'QA60', 'CPM']

            filtered_s2_scenes = parseNumberProducts(args.np, s2_scenes)
            if len(filtered_s2_scenes) != 0:
                print(str(len(filtered_s2_scenes)) +
                      ' S2 scenes to download'+'\n')
                for filtered_s2_scene in filtered_s2_scenes:
                    # print(s2_titles[s2_scenes.index(filtered_s2_scene)])
                    print(
                        s2_titles[list(s2_scenes.values()).index(filtered_s2_scene)])
                    scene_path, title = downloadS2_2(
                        run_parameters['s2_product'], filtered_s2_scene, bands, data_path, user_esa, pass_esa)

    # DOWNLOAD AND PROCESS SCENES _________________________________________________________________
    if run_parameters['run_mode'] == 'dp':  # donwload mode
        # search for scenes
        createFolderCheck(output_search_folder)
        if run_parameters['l8_product'] != 'NONE':
            l8_scenes = startSearchForLandsat8(run_parameters)
        if run_parameters['s2_product'] != 'NONE':
            s2_scenes, s2_titles = startSearchForSentinel2(
                run_parameters, len(l8_scenes))

        if len(l8_scenes) != 0 or len(s2_scenes) != 0:
            # making folder structure for data processing
            data_path, sds_path = createFolderTree(output_data_folder)
        if len(l8_scenes) != 0:
            filtered_l8_scenes = parseNumberProducts(args.np, l8_scenes)
            if len(filtered_l8_scenes) != 0:
                print(str(len(filtered_l8_scenes)) +
                      ' Landsat scenes to download'+'\n')
                for scene_l8 in filtered_l8_scenes:
                    print(scene_l8['display_id'])
                processLandsat8Scenes(filtered_l8_scenes, data_path, sds_path, beaches_path, user_usgs, pass_usgs, run_parameters['run_mode'],
                                      water_index, thresholding_method, cloud_mask_level, morphology_method, kernel_size, run_parameters['beach_code_filter'])

        if len(s2_scenes) != 0:
            filtered_s2_scenes = parseNumberProducts(args.np, s2_scenes)
            if len(filtered_s2_scenes) != 0:
                print(str(len(filtered_s2_scenes)) +
                      ' S2 scenes to download'+'\n')
                for filtered_s2_scene in filtered_s2_scenes:
                    # print(s2_titles[s2_scenes.index(filtered_s2_scene)])
                    print(
                        s2_titles[list(s2_scenes.values()).index(filtered_s2_scene)])
                processSentinel2Scenes(run_parameters['s2_product'], filtered_s2_scenes, data_path, sds_path, beaches_path, user_esa, pass_esa,
                                       run_parameters['run_mode'], s2_titles, water_index, thresholding_method, cloud_mask_level, morphology_method, kernel_size, run_parameters['beach_code_filter'])

    # ONLY PROCESS (process previous dowloaded scenes) __________________________________
    if run_parameters['run_mode'] == 'op':
        # search for scenes in the data folder
        data_path_s2 = recursiveFolderSearch(output_data_folder, Path('s2'))
        data_path_l8 = recursiveFolderSearch(output_data_folder, Path('l8'))
        data_path_l9 = recursiveFolderSearch(output_data_folder, Path('l9'))

        list_of_total_scenes = []

        if data_path_s2 != '':
            list_of_paths_s2 = [Path(f.path)
                                for f in os.scandir(data_path_s2) if f.is_dir()]
            if list_of_paths_s2 != []:
                list_of_total_scenes += list_of_paths_s2

        if data_path_l8 != '':
            list_of_paths_l8 = [Path(f.path)
                                for f in os.scandir(data_path_l8) if f.is_dir()]
            if list_of_paths_l8 != []:
                list_of_total_scenes += list_of_paths_l8

        if data_path_l9 != '':
            list_of_paths_l9 = [Path(f.path)
                                for f in os.scandir(data_path_l9) if f.is_dir()]
            if list_of_paths_l9 != []:
                list_of_total_scenes += list_of_paths_l9

        if list_of_total_scenes != []:
            filtered_scenes = filterScenesInfolder(list_of_total_scenes)
            if filtered_scenes != []:
                print('')
                print('Scenes to be processed: ')
                print('')
                for filtered_scene in filtered_scenes:
                    # print(filtered_scene)
                    data_path = recursiveFolderSearch(
                        output_data_folder, Path(filtered_scene))
                    sds_path = recursiveFolderSearch(
                        output_data_folder, Path('sds'))
                    # print(data_path)
                    if 's2' in data_path:
                        processSentinel2SceneFromPath(data_path, filtered_scene, sds_path, beaches_path, water_index,
                                                      thresholding_method, cloud_mask_level, morphology_method, kernel_size, run_parameters['beach_code_filter'])
                    else:
                        processLandsatSceneFromPath(data_path, filtered_scene, sds_path, beaches_path, water_index,
                                                    thresholding_method, cloud_mask_level, morphology_method, kernel_size, run_parameters['beach_code_filter'])

            # removing intermediate data
    if ridap == '1':
        createFolderTree(output_data_folder)

    # end clock
    endClock()


if __name__ == '__main__':
    """ Entrance to the algorithm workflow:
        - gets the user arguments
        - inits the logger
        - runs the algorithm
    """
    # getting user parameters
    args = parse_args()
    # inits logger with a required level
    log_level = os.getenv('LOG_LEVEL')
    init_logger(log_level)
    logger.info('Starting SAET algorithm...'+'\n')
    # runs the main algorithm
    run_algo(args)
    logger.info('SAET algorithm have finished successfully.'+'\n')
    sys.exit(0)
