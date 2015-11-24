# mapnik-tiler
Prepare map data for web.

This script is a front end for mapnik's generate_tiles_multiprocess.py. Groups of images can be merged using a GDAL VRT and passed at the command line.

For eg: gdalbuildvrt -vrtnodata "0 0 0" tifs.vrt *.tif
