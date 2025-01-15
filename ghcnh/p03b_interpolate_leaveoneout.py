#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import numpy as np
import pandas as pd

from local_paths import *

## =============================== ##

YEARS = list(range(2011, 2020+1, 1))
RADIUS_KM = 500.
NDAYS_MIN = 180

## =============================== ##

def expand_df(df, dims, values):
	df = df.set_index(dims)
	multi_index = (pd.MultiIndex.from_product(
			iterables=values,
			names=dims))
	df = df.reindex(multi_index, fill_value=np.nan).reset_index()
	df = df.sort_values(by=dims, ascending=True).reset_index(drop=True)
	return df

## =============================== ##

stationlists = [\
	'GHCNh_stationlist_daymean_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc00_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc12_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN)\
	]

for stationlist in stationlists:

	label = stationlist.split('_')[2]
	print(label)

	## =============================== ##

	df_stations = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, stationlist))
	stations = df_stations['Station_ID'].unique()

	## =============================== ##

	ifile = 'GHCNh_pairwisedistances_{0:s}_{1:d}-{2:d}_180plus.csv'.format(label, YEARS[0], YEARS[-1])
	df_dist = pd.read_csv(os.path.join(DATAPATH_DISTANCES, ifile))
	df_dist = df_dist.loc[df_dist['Station_ID'].isin(stations), stations]

	## =============================== ##

	# Preprocess the distance matrix
	df_dist = df_dist.set_index('Station_ID')
	df_dist_filtered = df_dist.where(df_dist <= RADIUS_KM, 0.)
	np.fill_diagonal(df_dist_filtered.values, 0)

	# Convert to NumPy array
	distance_matrix = df_dist_filtered.values  # Shape: (n_stations, n_stations)

	# Calculate inverse distances, avoiding division by zero
	inverse_distances = np.where(distance_matrix > 0., 1. / distance_matrix, 0.)

	# Normalize the weights so that each row sums to one
	weights_sum = inverse_distances.sum(axis=1, keepdims=True)  # Shape: (n_stations, 1)
	normalized_weights = np.divide(inverse_distances, weights_sum, where=weights_sum > 0.)
	weights_sum = normalized_weights.sum(axis=1, keepdims=True)  # Shape: (n_stations, 1)
	normalized_weights = np.divide(normalized_weights, weights_sum, where=weights_sum > 0.)

	# Convert the normalized weights back to a DataFrame for alignment
	df_weights = pd.DataFrame(normalized_weights, index=df_dist_filtered.index, columns=df_dist_filtered.columns)

	## =============================== ##

	stations = list(df_dist.columns)

	for year in YEARS:

		print(f"\nProcessing Year: {year}")

		ifile = 'GHCNh_{0:d}_t2m_{1:s}.csv'.format(year, label)
		df = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
		df = df.loc[df['Station_ID'].isin(stations), :]

		dates = df['date'].unique()
		dates.sort()
		stations = df['station'].unique()
		stations.sort()
		df = expand_df(df, ['station', 'date'], [stations, dates])

		## =============================== ##

		df_pivot = df.pivot(index='date', columns='station', values='temperature')
		df_pivot = df_pivot.sort_index()

		# Ensure the order of stations matches between df_weights and df_pivot
		df_pivot = df_pivot.loc[:, df_weights.columns]

		# Convert temperature data to a NumPy array
		temperature_matrix = df_pivot.values  # Shape: (n_dates, n_stations)

		# Create masks for stations used in interpolation
		W_bool = (normalized_weights > 0.)  # (n_stations x n_stations)
		T_nan = np.isnan(temperature_matrix)  # (n_dates x n_stations)

		# Compute mask: set to NaN if any neighbor used for interpolation has NaN
		mask = (T_nan @ W_bool.T) > 0  # (n_dates x n_stations)

		# Perform matrix multiplication to get interpolated temperatures
		temperature_matrix = np.nan_to_num(temperature_matrix, nan=0.0)
		interpolated_temps = temperature_matrix @ normalized_weights.T  # Shape: (n_dates, n_stations)

		# Set interpolated_temps to NaN where mask is True
		interpolated_temps[mask] = np.nan

		# Convert interpolated_temps back to a DataFrame
		df_interpolated = pd.DataFrame(interpolated_temps, index=df_pivot.index, columns=df_pivot.columns)

		# Optional: Validate Interpolation by Comparing with Actual Temperatures
		# Melt the DataFrame to long format for comparison
		df_actual = df_pivot.reset_index().melt(id_vars='date', var_name='Station_ID', value_name='t2m')
		df_interpolated_long = df_interpolated.reset_index().melt(id_vars='date', var_name='Station_ID', value_name='t2m_interpolated')

		ofile = 'GHCNh_t2m_{0:s}_{1:d}.csv'.format(label, year)
		df_actual.to_csv(os.path.join(DATAPATH_INTERPOLATED, ofile), index=False)

		ofile = 'GHCNh_t2m_{0:s}_{1:d}_interpolated.csv'.format(label, year)
		df_interpolated_long.to_csv(os.path.join(DATAPATH_INTERPOLATED, ofile), index=False)
