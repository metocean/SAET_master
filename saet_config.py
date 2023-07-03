############################ CONFIGURATION FILE #######################################

import os
from pathlib import Path

# variables declared with os.getenv(key,default=value) are exposed to change its value
# from other scripts without modify the configuration file

# CREDENTIALS FOR IMAGERY SERVERS (You must set your own credentials.) ################
os.environ['USER_ESA'] = os.getenv('USER_ESA', '********')
os.environ['PASS_ESA'] = os.getenv('PASS_ESA', '********')
os.environ['USER_USGS'] = os.getenv('USER_USGS', '********')
os.environ['PASS_USGS'] = os.getenv('PASS_USGS', '********')

# HOME FOLDER #########################################################################
os.environ["SAET_HOME_PATH"] = str(
    Path(os.path.dirname(os.path.realpath(__file__))))

# REMOVE INTERMEDIATE DATA AFTER PROCESSING (RIDAP) ###################################
# 1 -> remove data, 0 -> ramain data
# This parameters doesn't affect SAET for local users.
os.environ['RIDAP'] = os.getenv('RIDAP', '0')

# AUXILIARY DATA ######################################################################
# Folder for auxiliary shapefiles needed for SAET (relative paths)
os.environ['AUX_DATA_FOLDER_PATH'] = str(
    Path(os.path.join(os.getenv("SAET_HOME_PATH"), 'aux_data')))
# Names of every auxiliary shapefile
os.environ['SHP_BEACHES_PATH'] = str(
    Path(os.path.join(os.getenv("AUX_DATA_FOLDER_PATH"), 'beaches.shp')))
os.environ['SHP_LANDSAT_GRID_PATH'] = str(
    Path(os.path.join(os.getenv("AUX_DATA_FOLDER_PATH"), 'landsat_grid.shp')))
os.environ['SHP_SENTINEL2_GRID_PATH'] = str(
    Path(os.path.join(os.getenv("AUX_DATA_FOLDER_PATH"), 'sentinel2_grid.shp')))

# LOGGING #############################################################################
# Level for logging
# {NOTSET:0, DEBUG:10, INFO:20, WARNING:30, ERROR:40, CRITICAL:50}
os.environ['LOG_LEVEL'] = os.getenv('LOG_LEVEL', '20')

# RESULTS #############################################################################
# Quicklook activation. {0: deactivated, 1: activated}
os.environ['QUICKLOOK'] = os.getenv('QUICKLOOK', '1')
# S2 Quicklook server. {1: copernicus, 2: external}
os.environ['S2_QUICKLOOK_SERVER'] = os.getenv('S2_QUICKLOOK_SERVER', '1')

# Way to see the results in search mode
# {CONSOLE:0, TXT FILE:1, JSON FILE:2}
os.environ['OUT_RES'] = os.getenv('OUT_RES', '0')
