import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon

path = '/Users/benni/Desktop/AWS_OpenNRW/'
file = 'PVs_NRW_Sherlock.csv'

pvs_1 = pd.read_csv(path + file)
print(pvs_1.head(5))

