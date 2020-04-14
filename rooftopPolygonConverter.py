from shapely.geometry import Polygon
import pandas as pd
from geo_utils import Building, ns_bldg, ns_citygml, ns_gml, convert_3D_2D
from lxml import etree

DIRECTORY = 'your directory'
FOLDER = 'Folder of 3D Data ZIP containing n gml xml files LoD 2 standard'
TOWN = FOLDER.split('_')[3]
import os
import glob
os.chdir(DIRECTORY + FOLDER)
rooftopData_key = {}
buildingData = {}
check_sum = 0
for f in glob.glob('*.gml'):
    FILENAME = f[:f.rfind('.')]
    FULLPATH = DIRECTORY + FOLDER + '/' + f
    print('Parsing GML file ', FULLPATH)
    try:
        CITYGML = etree.parse(FULLPATH)
        root = CITYGML.getroot()
    except Exception as e:
        print('File {} cannot be parsed due to XMLSyntaxError:'.format(FILENAME), e)
        continue

    cityObjects = []
    buildings = []

    for obj in root.getiterator('{%s}cityObjectMember'% ns_citygml):
        cityObjects.append(obj)


    for cityObject in cityObjects:
        for child in cityObject.getchildren():
            if child.tag == '{%s}Building' %ns_bldg:
                buildings.append(child)

    print('There are', len(buildings),'Building(s) in this CityGML file')
    buildingclasses = []
    check_sum += len(buildings)
    for b in buildings:
        id = b.attrib['{%s}id' %ns_gml]
        buildingclasses.append(Building(b, id))

    for bu in buildingclasses:
        buildingData[bu.id] = {'Building_ID': bu.id, 'City': bu.city, 'PostalCode': bu.postalCode,
                           'Street': bu.streetName, 'StreetNumber': bu.streetNumber,
                           'Gemeindeschluessel': bu.gemeindeschluessel,
                            'RoofData': bu.roofdata,
                            'WallData': bu.walldata, 'GroundData': bu.grounddata,
                           'Datenquelle_Dachhoehe': bu.datenquelle_dachhoehe,
                           'DatenquelleBodenhoehe': bu.datenquelle_bodenhoehe, 'DatenquelleLage': bu.datenquelle_lage,
                           'BuildingFunction': bu.bldg_function, 'RooftopType': bu.bldg_roofType,
                           'MeasuredHeight': bu.bldg_measuredHeight, 'SourceFile': FILENAME
                           }
    for i in buildingData:
        for j in buildingData[i]["RoofData"]:
            rooftopData = buildingData[i]["RoofData"][j]
            rooftopData_key[j] = {'Building_ID': buildingData[i]['Building_ID'],
                                  'City': buildingData[i]['City'],
                                  'PostalCode': buildingData[i]['PostalCode'],
                                  'Street': buildingData[i]['Street'],
                                  'StreetNumber': buildingData[i]['StreetNumber'],
                                  'Gemeindeschluessel': buildingData[i]['Gemeindeschluessel'],
                                  'RooftopType': buildingData[i]['RooftopType'],
                                   'RoofTopID': j, 'Area': rooftopData['area'],
                                   'Azimuth': rooftopData['azimuth'],
                                   'Tilt': rooftopData['tilt'],
                                  'RooftopPolygon': rooftopData['polygon'],
                                   'Source_file': FILENAME}

buildingData_DF = pd.DataFrame(buildingData).transpose()
rooftopData_DF = pd.DataFrame(rooftopData_key).transpose()
rooftopData_DF['RooftopPolygon'] = rooftopData_DF['RooftopPolygon'].apply(Polygon)
rooftopData_DF['RooftopPolygon_2d'] = convert_3D_2D(rooftopData_DF['RooftopPolygon'])

print('There should be ', check_sum, 'buildings available in the dataframe. There are ',
      len(rooftopData_DF.Building_ID.unique()), 'buildings available.')

os.chdir('/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/Building3DData')
rooftopData_DF.to_csv('rooftopData_{}.csv'.format(TOWN), index = False)
buildingData_DF.to_csv('buildingData_{}.csv'.format(TOWN), index = False, sep = ';')