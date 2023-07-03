Python------------------------
python 397 (64 bits)

Package            min Version
------------------ ---------
python-dateutil    2.8.2
numpy		   1.21.2
matplotlib         3.4.3
GDAL               3.3.1 (*)
sentinelsat        1.1.0
shapely		   1.7.1
pyshp		   2.1.3
scikit-image       0.18.3
scikit-learn	   1.0.2
scipy              1.7.1
networkx           2.6.2

(*) IMPORTANT NOTE:
------------------------------------------------------
GDAL installation with pip command can be problematic.
If errors occur during the installation of GDAL, try 
installing the corresponding wheel file, according to
your operative system (Windows or Linux).

These files are in the folder "gdal_python_wheels":

Windows: GDAL-3.3.3-cp39-cp39-win_amd64.whl
Linux: GDAL-3.4.1-cp39-cp39-manylinux_2_5_x86_64.manylinux1_x86_64.whl

The installation can be done using the pip command.



