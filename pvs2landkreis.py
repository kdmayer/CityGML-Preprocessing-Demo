import geopandas as gpd
import glob
import pandas as pd
from shapely import wkt
from shapely.errors import WKTReadingError

#SOURCE: https://www.opengeodata.nrw.de/produkte/geobasis/vkg/dvg/dvg1/ (Shapefile)
#DESCRIPTION https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/verwaltungskarten/verwaltungsgrenzen/index.html
landkreise = gpd.read_file('landkreise/reprojected_landkreise.shp', encoding = 'utf8', crs = {'init': 'epsg:4326'})
landkreise = landkreise[['GN', 'geometry']]

pv_all = pd.DataFrame()
path = 'AWS_OpenNRW/'
for file in glob.glob(path + '*.csv'):
    print(file)
    pv_helper = pd.read_csv(file, sep = '\t', names = ['identifier','PV_coords', 'polygon_coords'])
    pv_all = pv_all.append(pv_helper)

pv_all = pv_all.reset_index(drop = True)

geometry = gpd.GeoSeries(pv_all['PV_coords'].apply(wkt.loads)) #string parsing of coordinates
pv_all_points = gpd.GeoDataFrame(geometry = geometry, crs = {'init': 'epsg:4326'}) #creating center image point geometry

return_df = gpd.sjoin(pv_all_points, landkreise, how = 'inner', op = 'within') #spatial join to map image centroids with "Landkreis" in North-Rhine Westphalia

liste_key_kreis = landkreise['GN'].values.tolist()

for i in range(len(liste_key_kreis)):
    try:
        print(liste_key_kreis[i])
        helper_df_sjoin = return_df[return_df['GN'] == liste_key_kreis[i]] #filtering of joined DF for a specific Landkreis
        helper_df_pv_filtered = pv_all.iloc[helper_df_sjoin.index, :] #filtering of whole segmented PVs dataframe for the lines which lie in the Landkreis
        helper_df_pv_filtered['Landkreis'] = liste_key_kreis[i] #add Landkreis information
        helper_df_pv_filtered.to_csv('/Users/benni/Desktop/AWS_OpenNRW/Landkreise_filtered/PV_Landkreis_{}.csv'.format(liste_key_kreis[i])) #save as csv
        geometry_list = helper_df_pv_filtered['polygon_coords'].values.tolist()
        geometrys = []
        for geometry in geometry_list:
            try:
                geometrys.append(wkt.loads(geometry))
            except WKTReadingError:
                geometrys.append('corrupted geometry')
        helper_df_pv_filtered['geometry'] = pd.Series(geometrys)
        helper_df_pv_filtered = helper_df_pv_filtered[helper_df_pv_filtered['geometry'] != 'corrupted geometry']
        geometry = gpd.GeoSeries(helper_df_pv_filtered['geometry'])
        helper_gpd_pv_filtered = gpd.GeoDataFrame(helper_df_pv_filtered[['identifier', 'PV_coords', 'Landkreis']], crs={'init': 'epsg:4326'})
        helper_gpd_pv_filtered.to_file(driver='ESRI Shapefile', filename="/Users/benni/Desktop/AWS_OpenNRW/Landkreise_filtered_shapefiles/PV_Landkreis_{}.shp".format(liste_key_kreis[i]))
    except ValueError as e:
        print(e)