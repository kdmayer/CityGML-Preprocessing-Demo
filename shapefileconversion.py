import geopandas as gpd
import pandas as pd
from shapely import wkt
import glob
import os
from geo_utils import add_missing_addresses_to_rooftopdata

#create shapefile for rooftops from csv file
#-- VARS
CRS = {'init': 'epsg:4326'} #epsg for middle Europe
DIRECTORY = '/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/Building3DData/rooftopData/' #directory with previously created rooftop files
DIRECTORY_OUTPUT = 'Shapefile Output' #set output directory
#-- CODE
def converter(df, col_name = 'RooftopPolygon_2d', prefix = 'nan', city = 'nan', output_dir = DIRECTORY_OUTPUT, crs = {'init': 'epsg:4326'}):
    #convert string geometry into valid geometry
    geometry = gpd.GeoSeries(df[col_name]).apply(wkt.loads)
    #create geodataframe
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs = crs)
    #only use relevant columns
    gdf = gdf[['Area', 'Azimuth', 'Building_ID', 'City',
       'PostalCode', 'RoofTopID', 'RooftopType',
        'Street', 'StreetNumber', 'Tilt', 'geometry']]
    #for missing addresses find the nearest address in dataframe and use it as a proxy
    gdf = add_missing_addresses_to_rooftopdata(gdf)
    #create shapefile
    gdf.to_file(driver='ESRI Shapefile',
                filename= output_dir + "{prefix}_{city}.shp".format(prefix = prefix, city = city))


if __name__ == "__main__":
    os.chdir(DIRECTORY)
    #iterates over DIRECTORY to create
    for file in glob.glob('*.csv'):
        rooftopPolygonDF = pd.read_csv(DIRECTORY + file)
        city = file.split('_')[1][:-4]
        #convert CSV rooftopfile to shapefile, uses 2D rooftopPolygon as geometry; do not change this in the initialization of the function
        converter(rooftopPolygonDF, prefix = 'polygon', city = city, crs = CRS)
        print('Processed {}'.format(file))
