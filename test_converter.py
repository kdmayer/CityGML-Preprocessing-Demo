from lxml import etree
from tqdm import tqdm
from shapely.geometry import Polygon
from shapely import wkt
from geo_utils import Building, ns_bldg, ns_citygml, ns_gml, convert_3D_2D, add_missing_addresses_to_rooftopdata

import geopandas as gpd
import pandas as pd
import glob
import os

# Configuration

DATA_DIR = "/Users/Kevin/Projects/CityGML_Processing/data/"
TOWN = "TEST_GML"

rooftop_data = {}
building_data = {}

check_sum = 0

gml_file_paths = []

GML_DATA_DIR = DATA_DIR + "GML/"

for file_name in os.listdir(GML_DATA_DIR):
    if file_name.endswith(".gml"):
        gml_file_paths.append(GML_DATA_DIR + file_name)

for file_path in gml_file_paths:

    try:
        # Read in file and build tree
        CITYGML = etree.parse(file_path)
        root = CITYGML.getroot()

    except Exception as e:
        # Some files have corrupted XML codes; these files are skipped
        print(f'File {file_path} cannot be parsed due to XMLSyntaxError: {e}')

    city_objects = []
    buildings = []

    # Create a list of all the objects in the gml file
    for obj in root.getiterator('{%s}cityObjectMember' % ns_citygml):
        city_objects.append(obj)

    # Create a list of all the buildings in the gml file
    for city_object in city_objects:
        for child in city_object.getchildren():
            if child.tag == '{%s}Building' % ns_bldg:
                buildings.append(child)

    print(f'There are {len(buildings)} Building(s) in this CityGML file.')

    # Iterate over building classes (see geo_utils BuildingClass for more information)
    building_classes = []
    check_sum += len(buildings)

    for building in buildings:
        identifier = building.attrib['{%s}id' % ns_gml]
        building_classes.append(Building(building, identifier))

    print("Create building dictionary")
    for building_class in tqdm(building_classes):
        building_data[building_class.id] = {
            'Building_ID': building_class.id, 'City': building_class.city,
            'Street': building_class.streetName, 'StreetNumber': building_class.streetNumber,
            'Gemeindeschluessel': building_class.gemeindeschluessel,
            'RoofData': building_class.roofdata, 'WallData': building_class.walldata,
            'GroundData': building_class.grounddata, 'Datenquelle_Dachhoehe': building_class.datenquelle_dachhoehe,
            'DatenquelleBodenhoehe': building_class.datenquelle_bodenhoehe,
            'DatenquelleLage': building_class.datenquelle_lage,
            'BuildingFunction': building_class.bldg_function, 'RooftopType': building_class.bldg_roofType,
            'MeasuredHeight': building_class.bldg_measuredHeight, 'SourceFile': file_path.split("/")[-1]
        }

    print("Create rooftop dictionary")
    # Create rooftopDictionary with relevant rooftop information
    for building_key in tqdm(building_data):

        for roof_key in building_data[building_key]["RoofData"]:
            roof = building_data[building_key]["RoofData"][roof_key]

            rooftop_data[roof_key] = {
                'Building_ID': building_data[building_key]['Building_ID'],
                'City': building_data[building_key]['City'],
                'Street': building_data[building_key]['Street'],
                'StreetNumber': building_data[building_key]['StreetNumber'],
                'Gemeindeschluessel': building_data[building_key]['Gemeindeschluessel'],
                'RooftopType': building_data[building_key]['RooftopType'],
                'RoofTopID': roof_key, 'Area': roof['area'],
                'Azimuth': roof['azimuth'],
                'Tilt': roof['tilt'],
                'RooftopPolygon': roof['polygon'],
                'Source_file': file_path.split("/")[-1]
            }

# Convert dictionaries to dataframes
building_data_df = pd.DataFrame(building_data).transpose()
rooftop_data_df = pd.DataFrame(rooftop_data).transpose()
rooftop_data_df['RooftopPolygon'] = rooftop_data_df['RooftopPolygon'].apply(Polygon)

# Create 2D rooftoppolygon
rooftop_data_df['RooftopPolygon_2d'] = convert_3D_2D(rooftop_data_df['RooftopPolygon'])

# Check sum to validate quality of extraction
print(f'There should be {check_sum} buildings available in the dataframe. There are {len(rooftop_data_df.Building_ID.unique())} buildings available.')


# Create GeoJSON for rooftops from CSV
def convert_df_to_gdf(rooftop_df, crs={'init': 'epsg:4326'}):
    rooftop_df = rooftop_df.rename(columns={'RooftopPolygon_2d': 'geometry', 'Gemeindeschluessel': 'PostalCode'})
    gdf = gpd.GeoDataFrame(rooftop_df, geometry='geometry', crs=crs)

    # Only use relevant columns
    gdf = gdf[['Area', 'Azimuth', 'Building_ID', 'City',
               'PostalCode', 'RoofTopID', 'RooftopType',
               'Street', 'StreetNumber', 'Tilt', 'geometry']]

    # For missing addresses find the nearest address in dataframe and use it as a proxy
    gdf = add_missing_addresses_to_rooftopdata(gdf)

    return gdf


gdf = convert_df_to_gdf(rooftop_data_df)

# Create GeoJSON
filename = DATA_DIR + "GeoJSON/" + f"_{TOWN}.geojson"
gdf.to_file(filename, driver="GeoJSON")

