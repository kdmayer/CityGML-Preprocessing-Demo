from shapely.geometry import Polygon
import pandas as pd
from geo_utils import Building, ns_bldg, ns_citygml, ns_gml, convert_3D_2D
from lxml import etree

DIRECTORY = 'your directory'
FOLDER = 'Folder of 3D Data directory containing multiple gml xml files LoD 2 standard'
TOWN = FOLDER.split('_')[3] #creates suffix for administrative zone
import os
import glob
os.chdir(DIRECTORY + FOLDER)
rooftopData_key = {} #dictionary initialization
buildingData = {}
check_sum = 0
for f in glob.glob('*.gml'):#iterate over all GML files
    FILENAME = f[:f.rfind('.')]#extract filename
    FULLPATH = DIRECTORY + FOLDER + '/' + f
    print('Parsing GML file ', FULLPATH)
    try:
        CITYGML = etree.parse(FULLPATH) #read in file and build tree
        root = CITYGML.getroot()
    except Exception as e:
        print('File {} cannot be parsed due to XMLSyntaxError:'.format(FILENAME), e) #some files have corrupted XML codes; these files are skipped
        continue

    cityObjects = []
    buildings = []
    #create list of objects in gml file
    for obj in root.getiterator('{%s}cityObjectMember'% ns_citygml):
        cityObjects.append(obj)

    #create list of buildings in gml file (usually the same number as 3D data in NRW only contains buildings, in other regions these files may contain other objects
    for cityObject in cityObjects:
        for child in cityObject.getchildren():
            if child.tag == '{%s}Building' %ns_bldg:
                buildings.append(child)

    print('There are', len(buildings),'Building(s) in this CityGML file')
    #iterate over building classes (see geo_utils BuildingClass for more information
    buildingclasses = []
    check_sum += len(buildings)
    for b in buildings:
        id = b.attrib['{%s}id' %ns_gml]
        buildingclasses.append(Building(b, id))
    #extract relevant information from Buildingclasses and create another dictionary
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
    #create rooftopDictionary with relevant rooftop information
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
#convert dictionaries to dataframes
buildingData_DF = pd.DataFrame(buildingData).transpose()
rooftopData_DF = pd.DataFrame(rooftopData_key).transpose()
rooftopData_DF['RooftopPolygon'] = rooftopData_DF['RooftopPolygon'].apply(Polygon)
#create 2D rooftoppolygon
rooftopData_DF['RooftopPolygon_2d'] = convert_3D_2D(rooftopData_DF['RooftopPolygon'])

#check sum to validate quality of extraction
print('There should be ', check_sum, 'buildings available in the dataframe. There are ',
      len(rooftopData_DF.Building_ID.unique()), 'buildings available.')

os.chdir('/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/Building3DData')
rooftopData_DF.to_csv('rooftopData_{}.csv'.format(TOWN), index = False)
buildingData_DF.to_csv('buildingData_{}.csv'.format(TOWN), index = False, sep = ';')