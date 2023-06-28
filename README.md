# SAET_master
SHORELINE ANALYSIS AND EXTRACTION TOOL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7488654.svg)](https://doi.org/10.5281/zenodo.7488654)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## 1. INTRODUCTION
This software has been developed as part of the ECFAS (European Coastal Flood Awareness System) project by the Geo-Environmental Cartography and Remote Sensing Group (CGAT) at the Universitat Politècnica de València, Spain. It contains the core algorithm for shoreline extraction at a sub-pixel level. For detailed information on the algorithm, please refer to the following papers:

- Palomar-Vázquez, J.; Pardo-Pascual, J.E.; Almonacid-Caballer, J.; Cabezas-Rabadán, C. Shoreline Analysis and Extraction Tool (SAET): A New Tool for the Automatic Extraction of Satellite-Derived Shorelines with Subpixel Accuracy. *Remote Sens.* 2023, 15, 3198. https://doi.org/10.3390/rs15123198
- Pardo-Pascual, J.E., Almonacid-Caballer, J., Ruiz, L.A., Palomar-Vázquez, J. Automatic extraction of shorelines from Landsat TM and ETM multi-temporal images with subpixel precision. *Remote Sensing of Environment*, 123, pp. 1-11, 2012. https://doi.org/10.1016/j.rse.2012.02.024
- Pardo-Pascual, J.E., Sánchez-García, E., Almonacid-Caballer, J., Palomar-Vázquez, J.M., Priego de los Santos, E., Fernández-Sarría, A., Balaguer-Beser, Á. Assessing the Accuracy of Automatically Extracted Shorelines on Microtidal Beaches from Landsat 7, Landsat 8 and Sentinel-2 Imagery. *Remote Sensing*, 10(2), 326, 2018. https://doi.org/10.3390/rs10020326
- Sánchez-García, E., Palomar-Vázquez, J. M., Pardo-Pascual, J. E., Almonacid-Caballer, J., Cabezas-Rabadán, C., & Gómez-Pujol, L. An efficient protocol for accurate and massive shoreline definition from mid-resolution satellite imagery. *Coastal Engineering*, 103732, 2020. https://doi.org/10.1016/j.coastaleng.2020.103732

- CGAT: https://cgat.webs.upv.es
- ECFAS: https://www.ecfas.eu

## Copyright and License
This open-source software is distributed under the GNU license and is the copyright of the UPV. It has been developed within the framework of the European ECFAS project by the following authors: Jesús Palomar Vázquez, Jaime Almonacid Caballer, Josep E. Pardo Pascual, and Carlos Cabezas Rabadán.

Please note that this software is designed specifically for the automatic extraction of shorelines in pre-storm and post-storm events and is not intended for massive extraction purposes.

## How to Cite SAET
To cite SAET in your research, please use the following reference:

J. Palomar-Vázquez, J. Almonacid-Caballer, J.E. Pardo-Pascual, and C. Cabezas-Rabadán (2021).
SAET (V 1.0). Open-source code. Universitat Politècnica de València. http://www.upv.es


## 2. WORKFLOW

In this image, we can see the general workflow of the algorithm. The parameters can change depending on the expected results (see section 5).

![Alt text](https://github.com/jpalomav/SAET_master/blob/main/doc/images/workflow.jpg)


## 3. REQUIREMENTS

The tool uses Sentinel 2, Landsat 8 and Landsat 9 images as input. In this way, the first thing we need is to have and user and a password from Copernicus Scihub and USGS Landsat Explorer servers:

- **Access USGS Landsat Explorer service:** In this case, you need to do two things: register on the Landsat Explorer website and make a request to access the service “machine to machine” (m2m). For the first requirement, you must register on the website https://ers.cr.usgs.gov/register. Once you have your credentials, access the website https://earthexplorer.usgs.gov, and go to your profile settings. Click on the button “Go to Profile” and finally, on the option “Access Request”. There you can make a new request to the m2m service by filling out a form.
Once you have your credentials for both data sources providers, you can edit the file “saet_config.py” (see structure section) changing the asterisks with your own credentials:

```
os.getenv('USER_ESA', '******')
os.getenv('PASS_ESA', '******')
os.getenv('USER_USGS', '******')
os.getenv('PASS_USGS, '******')
```

## 4. FOLDER STRUCTURE

The folder SAET contains the following files and subfolders:
-	saet_run.py. Main script. This script must be executed in command line mode.
-	saet_tools.py. Module that contains all functions needed to run the algorithm.
-	polynomial.py. Module with extra functions for surface interpolations
-	saet_config. py file that contains several configuration variables (credentials for the satellite imagery servers, folder structure for results, etc.).
-	examples_of_use.txt. Text file that contains several examples of use of the tool.
-	aux_data. Folder containing:
      * beaches.shp. Shapefile with all European areas classified as beaches. Based on “Coastal uses 2018” dataset (https://land.copernicus.eu/local/coastal-zones/coastal-zones-2018). This shapefile contains a field named “BEACH_CODE”, that will be copied to the final shoreline shapefile. This file is used to focus the shoreline extraction just only in these areas.
      * landsat_grid.shp. Shapefile with all Landsat-8 footprints (https://www.usgs.gov/media/files/landsat-wrs-2-descending-path-row-shapefile). This shapefile has two fields named “PATH” and “ROW” that will be used to select a scene by its identifier.
      * sentinel2_grid.shp. Shapefile with all Sentinel-2 footprints (based on the .kml file https://sentinels.copernicus.eu/web/sentinel/missions/sentinel-2/data-products).
      * map_director.qgz. Project file for QGIS with the three previous .shp files. It is useful for planning the process of downloading scenes.
      * SAET.pdf. Document explaining the tool
      * roi.geojson. File in geojson format containing an example of ROI (region of interest) to be used for searching scenes. To make easier the creation of this file (if it is needed), you can visit this website: https://geojson.io
-	search_data. Folder containing the file for the result of the algorithm in “searching” mode. This file will have the format indicated in the configuration file (saet_config.py). In this way, the possible formats are html, txt or json.
-	landsatxplore2. This folder contains the needed files for a modified version of the landsatxplore API to access to the Landsat imagery from USGS server. This have been necessary to update this library to the new Collection 2 of Landsat, which includes Landsat-9 product.

