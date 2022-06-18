'''
Analyses data from the agrarmeteo.ch station at Strickhof, Lindau-Eschikon
and calculates the mean air temperatur and annual precipitation sum between
2004 and 2022
'''

import pandas as pd

# URL to query agrarmeteo.ch for the Strickhof weather station
# url = 'https://www.agrometeo.ch/de/meteorologie/data?stations=175&sensors=1%3Aavg,6%3Asum&from=2000-01-01&to=2022-04-04&scale=day&groupBy=station&measured=0'

fpath_weather_data = '../../shp//weather_data.csv'
df = pd.read_csv(fpath_weather_data, sep=';')

df['Datum'] = pd.to_datetime(df['Datum'])
df['year'] = df['Datum'].apply(lambda x: x.year)

# exclude current year (too few data)
del_idx = df[df.year == 2022].index
df.drop(index=del_idx, inplace=True)

mean_temp = df['LINDAU - Temperatur Durchschnitt +2 m (°C)'].mean()
# df['LINDAU - Temperatur Durchschnitt +2 m (°C)'].plot()

annual_precip = df[['year', 'LINDAU - Niederschlag (mm oder Liter/m2)']].groupby(by='year').agg('sum')
