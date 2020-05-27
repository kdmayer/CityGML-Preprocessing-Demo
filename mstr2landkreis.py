import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
filtering = False

#read in mastr with geocoded columns
mstr_geocoded = pd.read_csv('mstr_bing.csv', encoding = 'utf8')

#create geometry list
geometry = [(xy) for xy in zip(mstr_geocoded['1'].astype(float), mstr_geocoded['0'].astype(float))]

#convert each coordinate tuple to shapely.Point
geometry = [Point(x) for x in geometry]

#convert pandas DataFrame to geopandas DataFrame
mstr_geocoded = gpd.GeoDataFrame(mstr_geocoded[['EinheitMastrNummer', 'Standort', 'EegInbetriebnahmedatum',
                                               'InstallierteLeistung', 'VerknuepfteEinheit',
                                               ]], geometry = geometry, crs = {'init': 'epsg:4326'})

#read in administrative zone shape file
landkreise = gpd.read_file('landkreise/reprojected_landkreise.shp', encoding = 'UTF-8', crs = {'init': 'epsg:4326'})
landkreise = landkreise[['GN', 'geometry']]

#create spatial join between geocoded MaStR and administrative zone polygons --> determination in which AZ MaStR points lie
mstr_geocoded = gpd.sjoin(mstr_geocoded, landkreise, how = 'inner', op = 'within')

admin_zones = mstr_geocoded.GN.unique().tolist()

#iterate over administrative zones and create distinct mastr file per administrative zone
for zone in sorted(admin_zones):
    print(zone)
    mstr_geocoded_helper = mstr_geocoded[mstr_geocoded.GN == zone]
    mstr_geocoded_helper['hash'] = mstr_geocoded_helper.geometry.apply(lambda geom: geom.wkb)
    hashes = mstr_geocoded_helper[mstr_geocoded_helper['hash'].duplicated(keep=False)]['hash'].unique().tolist()
    mstr_helper_corrected = pd.DataFrame()
    print(mstr_geocoded_helper.shape)
    #fitler most likely corrupted entries, set variable in the beginning; filtering is slow
    #FILTER: filters entries with the same installation capacity installed at the same address in the same year and month
    if filtering:
        duplicated = pd.Series()
        for hasher in hashes:
            # filter mstr with Geohash
            adress_helper = mstr_geocoded_helper[mstr_geocoded_helper['hash'] == hasher]
            # create rounded installation capacity column
            adress_helper['rounded'] = round(adress_helper['Installier'])
            # prepare iteration over years and months
            years = adress_helper['EegInbetri'].dt.year.unique().tolist()
            months = adress_helper['EegInbetri'].dt.month.unique().tolist()
            # iterate over years
            for year in years:
                datehelper = adress_helper[adress_helper['EegInbetri'].dt.year == year]
                for month in months:
                    # iterate over months in year
                    datehelper = datehelper[datehelper['EegInbetri'].dt.month == month]
                    # only consider multiple entries at an address in the same year and month
                    if len(datehelper.index) > 1:
                        # extract duplicates
                        duplicated = duplicated.append(datehelper['rounded'].duplicated(keep='first'))
            del adress_helper['rounded']
        indices = duplicated[duplicated == True].index #extract duplicate indices, keep first
        mstr_geocoded_helper = mstr_geocoded_helper[mstr_geocoded_helper.index[indices]] #filter out duplicate indices
        mstr_helper_corrected.to_file('mstr_split_landkreise/mstr_{}.shp'.format(zone), encoding = 'ISO-8859-1')
    else:
        mstr_geocoded_helper.to_file('mstr_split_landkreise/mstr_{}.shp'.format(zone), encoding = 'ISO-8859-1')







