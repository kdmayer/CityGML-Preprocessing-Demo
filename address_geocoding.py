import pandas as pd
from geo_utils import get_coordinates
import csv
from geopy.exc import GeocoderTimedOut, GeocoderQuotaExceeded
import time
import os

bnetzmstr_nrw = pd.read_csv('/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/bnetzmstr_nrw_filtered.csv', sep = ';')
adressen = bnetzmstr_nrw['Standort'].values.tolist()

if os.path.exists("coords.csv"):
    coords = pd.read_csv('coords.csv')
    adressen = adressen[coords.shape[0]:]
else:
    print('No coords.csv in directory - geocoding starts')



for i in range(len(adressen)):
    try:
        line = get_coordinates(adressen[i])
        with open('./coords.csv', "a") as csv_file:
            writer = csv.writer(csv_file, lineterminator="\n")
            writer.writerow(line)
            print(line)
            time.sleep(2)
    except (Exception, GeocoderQuotaExceeded, GeocoderTimedOut) as e:
        print(e)
        with open('./coords.csv', "a") as csv_file:
            writer = csv.writer(csv_file, lineterminator="\n")
            writer.writerow(',')
        time.sleep(10)
        continue
