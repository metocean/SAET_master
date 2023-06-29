# SAET
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
os.environ['USER_ESA'] = os.getenv('USER_ESA', '********')
os.environ['PASS_ESA'] = os.getenv('PASS_ESA', '********')
os.environ['USER_USGS'] = os.getenv('USER_USGS', '********')
os.environ['PASS_USGS'] = os.getenv('PASS_USGS', '********')
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

## Configuration file

The configuration file (saet_config.py) contains some sections that allows controlling several aspects to access the imagery servers and modify the algorithm workflow. Normally it will not be needed to change this file apart from the credential values, but if you want to do it, take in account this explanation about each section:
-	Section “credentials”. **Change the asterisks characters by your own credentials to run SAET properly (see section 3).**
-	Section “home folder”. It represents the path where SAET will be installed. All other subfolders will depend on it with relative paths.
-	Section “auxiliary data”. Relative path to the auxiliar data needed to SAET. The name of each shapefile can be changed if it is required.
-	Section “logging”. This section should be changed only by expert users. It controls the level of messages (severity) that SAET can return. For testing and debugging purposes, set this level to 10.
-	Section results. It controls how the user can see the searching results.
The rest of parameters that controls SAET are exposed as command line parameters (see section 5)

## 5. RUNNING SAET

One way to do this is by opening a PowerShell window or a command window (cmd). In Windows, go to the search bar and type "powershell" or "cmd". Run the script saet_run.py with parameters:
```
python saet_run.py --parameter=value
```
**Note:** whether using one or the other method (powershel or cmd), it would be a good idea to open these tools as administrator.

## Parameters

|     Parameter    |     Description    |     Required    |     Usage / examples    |     Default value    |
|---|---|---|---|---|
|     --help    |     Shows   the tool help message.    |          |          |          |
|     --rm    |     Run   mode. “os” means only searching; “dp” downloading and processing; “od” only   downloading; “op” only processing (for previous downloaded images); “or”   offline S2 retrieval. If you want just only search for images, use the “os”   mode.     **Note:** “dp” mode must be   used along with the parameter --np (number of products). “or” mode must be   used along with the parameter –oa (offline S2 activation)    |     Yes    |     --rm=os     --rm=dp     --rm=od     --rm=op     --rm=or    |     os    |
|     --fp    |     Footprint   for scene searching. This parameter has three ways to be used:           - path   to the ROI file (see structure section).      - bounding   box in latitude and longitude coordinates separated by “comma” with the   format min_long, min_lat, max_long, max_lat.      -   using NONE value to avoid searching by ROI. In this case, we can activate the   searching by scene or tile identifiers (parameters --ll and --sl).          |     Yes   except “op” mode    |     fp=c:\data\roi.geojson.     fp= 0.28,39.23,0.39,39.33     fp =   NONE    |     NONE    |
|     --sd    |     Start   date for searching scenes in format (YYYYMMDD).    |     Yes   except “op” mode    |     --sd=20211001    |          |
|     --cd    |     Central   date for searching scenes in format (YYYYMMDD). It is assumed to be the   central date of the storm.    |     Yes   except “op” mode    |     --cd=20211001    |          |
|     --ed    |     End   date for searching scenes in format (YYYYMMDD).    |     Yes   except “op” mode    |     --ed=20211001.    |          |
|     --mc    |     Maximum   percentage of cloud coverage in the scene. It must be a number between 0 and   100.    |     Yes   except “op” mode    |     --mc=10    |          |
|     --lp    |     Product   type for Landsat scenes. By default, the tool uses the Collection 2 (level 1)   to search Landsat-8 images (landsat_ot_c2_l1), but it also can search Landsat   8 and Landsat 9 images from Collection 2 at level 2 (landsat_ot_c2_l2). This   parameter can be set up to “NONE”, to avoid Landsat scenes processing.     **Note:** -	L8 and L9 products are only available from Collection 2. Collection 1 are not available anymore as of December 30, 2022.    |     Yes   except “op” mode    |     --lp=   landsat_ot_c2_l1     --lp=   landsat_ot_c2_l2     --lp=NONE    |          |
|     --ll    |     Scene   list identifiers for Landsat images. It must be a list of numbers of 6   digits. If there is more than one identifier, they must be separated by the   “comma” character. The default value is NONE (which means that the search by   ROI will have priority).    |     Yes   except “op” mode    |     --ll=198032     --ll=198032,199031     --ll=NONE    |     NONE    |
|     --sp    |     product   type for Sentinel-2 scenes (S2). The tool uses 1C (S2MSI1C) and 2A (S2MSI2A)   products. The product by default is 1C. This parameter can be set up to   “NONE”, to avoid S2 scenes processing.    |     Yes   except “op” mode    |     --sp=   S2MSI1C     --sp=   S2MSI2A     --sp=NONE    |     NONE    |
|     --sl    |     Scene   list identifiers for S2. It must be a list of alphanumeric characters (named   “tile” identifier) composed of two numbers and three capital letters. If   there is more than one identifier, they must be separated by the “comma”   character. The default value is NONE (which means that the search by ROI will   have priority).    |     Yes   except “op” mode    |     --sl=30TYJ     --sl=30TYJ,30TYK     --sl=NONE    |     NONE    |
|     --bc    |     Beach   code list to filter the extraction process for a group of beaches. This code   is related to the “BEACH_CODE” field in the shapefile “Beaches.shp”. The   default value is NONE, which means that all beaches in the image will be   processed. In case you want process a group of beaches, you must indicate a   list of codes, separated by “comma”. If some code in the list is not correct,   it will not be considered. If all codes are incorrect, all beaches will be   processed.    |     False    |     --bc=1683     --bc=1683,2485,758     --bc=NONE    |     NONE    |
|     --of    |     Output   data folder. It indicates the folder where the images and results will be   stored    |     No    |     --of=c:\data   (windows)     --of=/data(linux)    |     SAET installation folder    |
|     --wi    |     Water   index type. SAET supports these indices: aweinsh, aweish, mndwi, kmeans   (K-means it is not a water index, but also leads to a classification mask. In   this case it is not needed a threshold value).    |     No    |     --wi=aweinsh     --wi=aweish     --wi=mndwi     --wi=kmeans    |     aweinsh    |
|     --th    |     Threshold   method to obtain the water-land mask from the water index. SAET supports   three methods: standard 0 value, Otsu bimodal and Otsu multimodal with three   classes. These methods are applied for all type of index except kmeans.    |     No    |     --th=0   (standard 0 value)     --th=1   (Otsu bimodal)     --th=2   (Otsu multimodal)    |     0    |
|     --mm    |     Morphological   method. To generate the shoreline at pixel level (SPL) from the water-land   mask. SAET can apply two methods: erosion and dilation.    |     No    |     --mm=erosion     --mm=dilation    |     dilation    |
|     --cl    |     Cloud   masking severity. This parameter controls what kind of clouds will be used to   mask the SPL. SAET supports three levels of severity: low (SAET don’t use   cloud mask), medium (only opaque clouds) and high (opaque clouds, cirrus, and   cloud shadows).     **Note:**   Landsat   8-9 Collection 2 and Sentinel-2 use algorithms to classify clouds. SAET uses   these classification layers. You must assume that sometimes this   classification can fail. This will directly affect the result.    |     No    |     --cl=0   (low)     --cl=1   (medium)     --cl=2   (high)    |     0    |
|     --ks    |     Kernel   size. The main algorithm for shoreline extraction uses a kernel analysis over   each pixel in the SPL. Users can control this size, choosing between 3 or 5   pixels.    |     No    |     --ks=3     --ks=5    |     3    |
|     --np    |     Number   of products. List of products to be processed in both “d” and “p” modes. The   list must consist of the identifiers of the products to be processed. These   identifiers can be obtained in the searching mode. This parameter supports   three different formats: list of identifiers, range of identifiers or all   identifiers.    |     No    |     --np=0,2,5,3   (list)     --np=5-10   (range)     --np=*   (all)    |     NONE    |
|     --oa    |     Offline   Sentinel-2 activation. Some S2 products could be in offline mode. So, before   downloading them, they must be activated (from offline to online mode). See “considerations”   section.    |     No   except “or” mode    |     --oa=check     --oa=activate    |     check    |

This is the text of help that appears when you run SAET with the --h parameter (python saet_run.py --h):
```
usage: saet_run.py [-h] --rm {os,dp,od,op,or} --fp FP --sd SD --cd CD --ed ED --mc [0-100] --lp {landsat_ot_c2_l1,landsat_ot_c2_l2,NONE} --ll LL --sp {S2MSI1C,S2MSI2A,NONE} --sl SL [--bc BC] [--of OF]
                   [--wi {aweish,aweinsh,mndwi,kmeans}] [--th {0,1,2}] [--mm {erosion,dilation}] [--cl {0,1,2}] [--ks {3,5}] [--np NP] [--oa {check,activate}]

optional arguments:
  -h, --help            show this help message and exit
  --rm {os,dp,od,op,or}
                        Run mode (only search [s] / download and process [dp] / only donwload [od] / only process [op] / offline S2 retrieval [or]). --rm=os / --rm=dp / --rm=od / --rm=op / --rm=or.
                        Default: os
  --fp FP               path of the roi file for searching scenes (fp=c:\data oi.geojson), coordinates long/lat in this format: fp=long_min,lat_min,long_max,lat_max. Default: NONE
  --sd SD               Start date for searching scenes (YYYYMMDD). --sd=20210101. Default:20200101
  --cd CD               Central date for storm (YYYYMMDD). --sd=20210101. Default:20200102
  --ed ED               End date for searching scenes (YYYYMMDD). --sd=20210101. Default:20200103
  --mc [0-100]          maximum cloud coverture for the whole scene [0-100]. --mc=10
  --lp {landsat_ot_c2_l1,landsat_ot_c2_l2,NONE}
                        Landsat 8 product type. landsat_ot_c2_l1 or landsat_ot_c2_l2 or NONE. Default: landsat_ot_c2_l1
  --ll LL               List of scenes for Landsat 8 (number of 6 digits). --ll=198032,199031. Default: NONE
  --sp {S2MSI1C,S2MSI2A,NONE}
                        Sentinel 2 product type (S2MSI1C / S2MSI2A). --s2=S2MSI1C / --s2=S2MSI2A / NONE. Default: S2MSI1C
  --sl SL               List of scenes for Sentinel 2 (string of 5 characters). --sl=31TCF,30TYK. Default: NONE
  --bc BC               beach code filter list. --bc=520,548 Default: NONE
  --of OF               output data folder. --of=c:\data (windows) --of=/data. Default: SAET_HOME_PATH
  --wi {aweish,aweinsh,mndwi,kmeans}
                        Water index type (aweish, aweinsh,mndwi,kmeans). --wi=aweinsh. Default: aweinsh
  --th {0,1,2}          Thresholding method (0: standard 0 value, 1: Otsu bimodal, 2: Otsu multimodal 3 classes). --th=0. Default: 0
  --mm {erosion,dilation}
                        Morphological method (erosion, dilation). --mm=dilation, Default: dilation
  --cl {0,1,2}          Cloud mask level (0: no masking, 1: only opaque clouds, 2: opaque clouds + cirrus + cloud shadows). Default: 0
  --ks {3,5}            Kernel size for points extraction. Default: 3
  --np NP               List of number of products for download (only if --rm=d and --rm=p). [0,2,5,3] / [*] / [5-10]. Default: NONE
  --oa {check,activate}
                        Offline S2 activation (only if --rm=or). "check" / "activate". Default: "check"
```

## Examples of use

* Searching for all Sentinel-2 (level 1C) scenes inside an area of interest (tile 30SYJ) and with less than 15% of cloud coverage:
```
python saet_run.py --rm=os --fp=NONE --sd=20230401 --cd=20230415 --ed=20230430 --mc=15 --lp=NONE --ll=NONE --sp=S2MSI1C --sl=30SYJ

2023-06-28 16:29:02,542 INFO Starting SAET algorithm...

2023-06-28 16:29:02,880 INFO Searching for S2 images...

2023-06-28 16:29:04,066 INFO Found 6 products
Scene: S2B_MSIL1C_20230420T104619_N0509_R051_T30SYJ_20230420T125145 Cloud coverage: 9.2 availability: online
Scene: S2B_MSIL1C_20230417T103629_N0509_R008_T30SYJ_20230417T123856 Cloud coverage: 0.05 availability: online
Scene: S2A_MSIL1C_20230412T103621_N0509_R008_T30SYJ_20230412T155912 Cloud coverage: 0.0 availability: online
Scene: S2B_MSIL1C_20230407T103629_N0509_R008_T30SYJ_20230407T124138 Cloud coverage: 0.0 availability: online
Scene: S2A_MSIL1C_20230405T105031_N0509_R051_T30SYJ_20230405T160934 Cloud coverage: 0.01 availability: online
Scene: S2A_MSIL1C_20230402T103621_N0509_R008_T30SYJ_20230402T160048 Cloud coverage: 0.0 availability: online


[0] Scene: S2B_MSIL1C_20230420T104619_N0509_R051_T30SYJ_20230420T125145 Cloud coverage: 9.2% 5 days Online
[1] Scene: S2B_MSIL1C_20230417T103629_N0509_R008_T30SYJ_20230417T123856 Cloud coverage: 0.05% 2 days Online
[*******] Central date:20230415
[2] Scene: S2A_MSIL1C_20230412T103621_N0509_R008_T30SYJ_20230412T155912 Cloud coverage: 0.0% -3 days Online
[3] Scene: S2B_MSIL1C_20230407T103629_N0509_R008_T30SYJ_20230407T124138 Cloud coverage: 0.0% -8 days Online
[4] Scene: S2A_MSIL1C_20230405T105031_N0509_R051_T30SYJ_20230405T160934 Cloud coverage: 0.01% -10 days Online
[5] Scene: S2A_MSIL1C_20230402T103621_N0509_R008_T30SYJ_20230402T160048 Cloud coverage: 0.0% -13 days Online


2023-06-28 16:29:06,192 INFO Time passed: 0hour:0min:3sec

2023-06-28 16:29:06,192 INFO SAET algorithm have finished successfully.
```
Results show a first list of products with its availability (online or offline). For every image, this last list contains the identifier (number of order in the list), the name, the cloud coverage percentage, and the difference in days between the central date and de image date. Besides, depending on the values of the variable “QUICKLOOK” (configuration file), an .html file with the quicklook images is opened automatically in the default browser. 

**Note:** sometimes firefox may experiment problems showing quicklooks. If this is the case, try chrome as default browser.

Once the searching results have been obtained, if we want to download and process just only the closest products to the central date (products 1 and 2), the suitable command line will be the following:
```
python saet_run.py --rm=dp --fp=NONE --sd=20230401 --cd=20230415 --ed=20230430 --mc=15 --lp=NONE --ll=NONE --sp=S2MSI1C --sl=30SYJ --np=1,2
```

In case we want to reprocess any previous downloaded image, for example using other water index or other thresholding method, we only have to change these parameters and repeat the same command line but changing the run mode (--rm) to “op”. In this case, only the images previously downloaded will be reprocessed.

```
python saet_run.py --rm=op --wi=mndwi
2023-06-28 17:11:26,536 INFO Starting SAET algorithm...

List of scenes in the data folder:

[0] S2A_MSIL1C_20230405T105031_N0509_R051_T30SYJ_20230405T160934
[1] S2B_MSIL1C_20230420T104619_N0509_R051_T30SYJ_20230420T125145

Number of images to be reprocessed?: 0,1
```
SAET will display a list of images already stored in the output data folder and will ask for the images to be processed. You must use the same syntax as the --np parameter.

**Note:** more examples can be found in the file “examples_of_use.txt”.

## run mode election

Next picture shows the workflow to run SAET in the most convenient way. The recommendation is:
* Select your area of analysis and product of interest. The file map_director.qgz (QGIS) will be very useful to decide what scene (Landsat) or tile (Sentinel-2) will be 
* Always start with the "only searching mode".
* If you are going to analyse just a few images (one or two images before and after the storm peak date), you can follow with the "downloading and processing" run mode.
* In case you want to analyse a time series and in order to prevent conexion problems form the serves, it is recommendable to use "only download" run mode instead "downloading and processing".
* Anyway, if you want to reprocess any image previously downloaded, you can use the "only processing" run mode. This run mode will allow you to reprocess the images with other parameters (water index, threshol method, etc.) without having to download them again.
* Only in case you want to analyse Sentinel-2 images and some of them are in "offline" mode, you can use the "Offline S2 activation" run mode, along with the parameter --oa, first time with the value "--oa=activate" to activate the product, and from time to time, with the value "--oa=check", to check if the products are online. **Note:** Only S2 online products can be downloaded. Once a product has been activated, remains in online mode during just few days. Finally, if you prefer, you also can do the activation process from the Landsat Earth Explorer platform.

![Alt text](https://github.com/jpalomav/SAET_master/blob/main/doc/images/run_modes.jpg)
