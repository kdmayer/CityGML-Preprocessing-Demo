# Geospatial Tools for the Processing of 3D-Geometries and related files

This repository deals with the processing of 3D-geometries provided by the state of North Rhine-Westphalia ([here](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/3d_gebaeudemodelle/index.html))

## Dependencies

The repository requires the packages GeoPandas, scipy, geopy, pandas, numpy, and related typical packages used. 


## File Descriptions

```bash
address_geocoding.py
```
This script provides a script to directly geocode addresses using OSM API. 
```bash
geo_utils.py
```
This file is a collection of relevant helper functions for the execution of the 3D geoprocessing and other processes.
```bash
mastr2landkreis.py
```
This script leverages the polygons of the administrative districts in Nord Rhine-Westphalia to create subsets of the public solar registry within those polygons. 
```bash
pvs2landkreis.py
```
This script leverages the polygons of the administrative districts in Nord Rhine-Westphalia to create subsets of the segmented PV polygons within those polygons. 
```bash
rooftopPolygonConverter.py
```
Converts 3D CityGML files to two kinds of CSV-Files: First, all relevant features are extracted to create a building dictionary containing wall geometries, layout geometries, and rooftop geometries. This is saved as a CSV-File. Second, these files are filtered for rooftop geometries, which are then processed in the next step.
```bash
shapefileconversion.py
```
Converts rooftop-CSV files to shapefiles for geoprocessing purposes. 


## License
[MIT](https://choosealicense.com/licenses/mit/)
