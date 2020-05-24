import geopandas as gpd
import glob
import pandas as pd
from shapely import wkt
import time

#SOURCE: https://www.opengeodata.nrw.de/produkte/geobasis/vkg/dvg/dvg1/ (Shapefile)
#DESCRIPTION https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/verwaltungskarten/verwaltungsgrenzen/index.html
start = time.time()
landkreise = gpd.read_file('landkreise/reprojected_landkreise.shp', encoding = 'utf8', crs = {'init': 'epsg:4326'}) #read in administrative zone polygons
landkreise = landkreise[['GN', 'geometry']]


pv_all = pd.DataFrame()
path = 'AWS_OpenNRW/'
for file in glob.glob(path + '*.csv'): #read in of all segmented PV polygon files in the defined path
    '''
    this loop takes very long and requires large computing capacities
    only execute this code if necessary
    '''
    print(file)
    pv_helper = pd.read_csv(file, sep = '\t', usecols = [1,2], header = None) #only use relevant cols to reduce memory required
    pv_all = pv_all.append(pv_helper)

interim_1 = time.time()

print('Loaded PV polygon files. This took {} seconds'.format(interim_1 - start))

pv_all = pv_all.rename(columns = {1: 'PV_coords', 2: 'polygon_coords'}) #name columns for indexing
pv_all = pv_all[pv_all['polygon_coords'] != 'POLYGON ((6.341371801241112 51.82340640781572, 6.341(8.070480264504688, 51.11216465078038, 8.07391834953523, 51.114323022634586)'] #exclude corrupted geometry

pv_all = pv_all.reset_index(drop = True)
geometry = gpd.GeoSeries(pv_all['PV_coords'].apply(wkt.loads)) #conversion of string coordinates to shapely Polygon
pv_all = gpd.GeoDataFrame(pv_all['polygon_coords'], geometry = geometry, crs = {'init': 'epsg:4326'}) #creating center image point geometry
return_df = gpd.sjoin(pv_all, landkreise, how = 'inner', op = 'within') #spatial join image center with "Landkreis" in North-Rhine Westphalia
interim_3 = time.time()
print('Spatial join to administrative zones. This took {} seconds'.format(interim_3 - start))
liste_key_kreis = landkreise['GN'].values.tolist()
del return_df['index_right']
for i in range(len(liste_key_kreis)): #creating PV polygon file per administrative zone
    try:
        print(liste_key_kreis[i])
        helper_df_sjoin = return_df[return_df['GN'] == liste_key_kreis[i]] #filtering of joined DF for a specific Landkreis
        pd.DataFrame(helper_df_sjoin).to_csv('/Users/benni/PycharmProjects/geospatial_tools/AWS_OpenNRW/Landkreise_filtered/PV_Landkreis_{}.csv'.format(liste_key_kreis[i])) #save as csv
        print('Shape:', helper_df_sjoin.shape[0])
        helper_df_sjoin.geometry = gpd.GeoSeries(helper_df_sjoin['polygon_coords'].apply(wkt.loads))#convert string to polygon
        del helper_df_sjoin['polygon_coords']
        duplicate_geometry = helper_df_sjoin.geometry #helper GeoSeries for duplicate dropping
        duplicate_geometry = duplicate_geometry.centroid #using centroid of respective polygon for duplicate search
        duplicate_geometry = duplicate_geometry.apply(lambda geom: geom.wkb) #create a hash of centroid point
        helper_df_sjoin = helper_df_sjoin.loc[duplicate_geometry.drop_duplicates().index] #find duplicate hashes and filter dataframe
        print('Shape after deleting duplicates:', helper_df_sjoin.shape[0])
        helper_df_sjoin.to_file(driver='ESRI Shapefile', filename="/Users/benni/PycharmProjects/geospatial_tools/AWS_OpenNRW/Landkreise_filtered_shapefiles/PV_Landkreis_{}.shp".format(liste_key_kreis[i]))
    except ValueError as e:
        print(e)

interim_4 = time.time()

print('The whole conversion took {} seconds'.format(interim_4 - start))