import geopandas as gpd
import pandas as pd
from shapely import wkt
import glob
import os
from geo_utils import add_missing_addresses_to_rooftopdata

#-- VARS
CRS = {'init': 'epsg:4326'} #epsg for middle Europe
DIRECTORY = '/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/Building3DData/rooftopData/'
DIRECTORY_OUTPUT = 'Shapefile Output'
#-- CODE
def converter(df, col_name = 'RooftopPolygon_2d', prefix = 'nan', city = 'nan', output_dir = DIRECTORY_OUTPUT, crs = {'init': 'epsg:4326'}):
    geometry = gpd.GeoSeries(df[col_name]).apply(wkt.loads)
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs = crs)
    gdf = gdf[['Area', 'Azimuth', 'Building_ID', 'City',
       'PostalCode', 'RoofTopID', 'RooftopType',
        'Street', 'StreetNumber', 'Tilt', 'geometry']]
    gdf = add_missing_addresses_to_rooftopdata(gdf)
    gdf.to_file(driver='ESRI Shapefile',
                filename= output_dir + "{prefix}_{city}.shp".format(prefix = prefix, city = city))


if __name__ == "__main__":
    os.chdir(DIRECTORY)
    for file in glob.glob('*.csv'):
        rooftopPolygonDF = pd.read_csv(DIRECTORY + file)
        city = file.split('_')[1][:-4]
        converter(rooftopPolygonDF, prefix = 'polygon', city = city, crs = CRS)
        print('Processed {}'.format(file))
