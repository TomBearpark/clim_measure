#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import numpy as np
import pandas as pd
from geopy.distance import geodesic

from local_paths import *

## =============================== ##

YEARS = list(range(2011, 2020+1, 1))
NDAYS_MIN = 180

## =============================== ##

def compute_geodesic_distance(pair, latitudes, longitudes):
	idx1, idx2 = pair
	coord1 = (latitudes[idx1], longitudes[idx1])
	coord2 = (latitudes[idx2], longitudes[idx2])
	return geodesic(coord1, coord2).kilometers

## =============================== ##

stationlists = [\
	'GHCNh_stationlist_daymean_2011-2020_{0:03d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc00_2011-2020_{0:03d}plus.csv'.format(NDAYS_MIN),
	#'GHCNh_stationlist_utc00utc12_2011-2020_{0:03d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc12_2011-2020_{0:03d}plus.csv'.format(NDAYS_MIN)\
	]

for stationlist in stationlists:

	label = stationlist.split('_')[2]

	df = pd.read_csv(os.path.join(DATAPATH_STATIONLISTS, stationlist))

	n_stations = df.shape[0]
	latitudes = df['Latitude'].values
	longitudes = df['Longitude'].values
	distance_matrix = np.zeros((n_stations, n_stations), dtype=np.float32)

	for i in range(n_stations):
		print(i)
		coord_i = (latitudes[i], longitudes[i])
		for j in range(i + 1, n_stations):
			coord_j = (latitudes[j], longitudes[j])
			distance = geodesic(coord_i, coord_j).kilometers
			distance_matrix[i, j] = distance
			distance_matrix[j, i] = distance  # Symmetric

	# Convert the NumPy array to a DataFrame
	df_dist = pd.DataFrame(distance_matrix, index=df['Station_ID'], columns=df['Station_ID'])

	# Optionally, save to a CSV file
	ofile = 'GHCNh_pairwisedistances_{0:s}_{1:d}-{2:d}_{3:03d}plus.csv'.format(label, 2011, 2020, NDAYS_MIN)
	df_dist.to_csv(os.path.join(DATAPATH_DISTANCES, ofile))
