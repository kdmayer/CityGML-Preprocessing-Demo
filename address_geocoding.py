import pandas as pd
from geo_utils import get_coordinates
import csv
from geopy.exc import GeocoderTimedOut, GeocoderQuotaExceeded
import time
edit = False

bnetzmstr_nrw = pd.read_csv('/Users/benni/PycharmProjects/DeepSolar_3D/3d_data/bnetzmstr_nrw_filtered.csv', sep = ';')

#bnetzmstr_nrw['Standort'].str.extract(r'(?P<Street>\D+)\s?(?P<Number>\d+\S*)?\s(?P<Postal>\d{5})\s(?P<City>\D+)$')
adressen = bnetzmstr_nrw['Standort'].values.tolist()

for i in range(len(adressen)):
    try:
        line = get_coordinates(adressen[i])
        with open('./coords.csv', "a") as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerow(line)
            print(line)
            time.sleep(2)
    except (Exception, GeocoderQuotaExceeded, GeocoderTimedOut) as e:
        print(e)
        time.sleep(10)
