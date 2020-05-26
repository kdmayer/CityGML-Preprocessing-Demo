import pandas as pd
from geo_utils import get_coordinates
import csv
from geopy.exc import GeocoderTimedOut, GeocoderQuotaExceeded
import time
import os

#load dataframe for geocoding
bnetzmstr_nrw = pd.read_csv('/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/bnetzmstr_nrw_filtered.csv', sep = ';')
#extract addresses to be geocoded from dataframe
adressen = bnetzmstr_nrw['Standort'].values.tolist()

#to secure results a csv file, which saves the results temporarily is created
if os.path.exists("coords.csv"):
    coords = pd.read_csv('coords.csv')
    adressen = adressen[coords.shape[0]:]
else:
    print('No coords.csv in directory - geocoding starts')


#iterate over addresses and retrieve geocodes using OSM API
for i in range(len(adressen)):
    try:
        line = get_coordinates(adressen[i]) #save coordinates
        with open('./coords.csv', "a") as csv_file:
            writer = csv.writer(csv_file, lineterminator="\n")
            writer.writerow(line) #write coordinates into csv file
            print(line)
            time.sleep(2)
    except (Exception, GeocoderQuotaExceeded, GeocoderTimedOut) as e:
        print(e)
        with open('./coords.csv', "a") as csv_file:
            writer = csv.writer(csv_file, lineterminator="\n")
            writer.writerow(',') #mark results, which could not be retrieved out of different exceptions (see except)
        time.sleep(10)
        continue
