#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import numpy as np
import pandas as pd

import geopandas as gpd

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from local_paths import *

## =============================== ##

YEARS = list(range(2011, 2020+1, 1))
NDAYS_MIN = 30

## =============================== ##

datapath = "/home/manuel/Research/projects/project_41/replication/data/"
datafile = "countries.shp"
gdf_countries = gpd.read_file(os.path.join(datapath, datafile))
gdf_countries = gdf_countries.loc[:, ['ISO_A3_EH', 'geometry']]

## =============================== ##

stationlists = [\
	'GHCNh_stationlist_daymean_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc00_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN),
	#'GHCNh_stationlist_utc00utc12_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN),
	'GHCNh_stationlist_utc12_2011-2020_{0:d}plus.csv'.format(NDAYS_MIN)\
	]

for stationlist in stationlists:

	label = stationlist.split('_')[2]

	df = pd.read_csv(os.path.join(DATAPATH_STATIONLISTS, stationlist))

	gdf_stations = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude)).set_crs('EPSG:4326')

	fig, ax = plt.subplots(figsize=(6,4))
	ax.set_title(None)
	gdf_countries.plot(ax=ax, alpha=1., facecolor='grey', lw=0.5, edgecolor='k')
	#gdf.plot(ax=ax, alpha=1., column=VARIABLE, lw=0.2, cmap=cmap, norm=norm, edgecolor='none', markersize=1.)
	gdf_countries.plot(ax=ax, alpha=1., facecolor='none', lw=0.5, edgecolor='k')
	gdf_stations.plot(ax=ax, marker='o', markersize=0.2, color='r')
	ax.set_xlim(-130., 180.)
	ax.set_ylim(-60., 75.)
	ax.plot(ax.get_xlim(), [-23.5, -23.5], 'k--', linewidth=0.5)
	ax.plot(ax.get_xlim(), [23.5, 23.5], 'k--', linewidth=0.5)
	ax.set_xlabel(None)
	ax.set_ylabel(None)
	#ax.annotate(text='GHCNh stations with at least 180 observations\n for 12 UTC and 0 UTC every year 2011-2020', xy=(0.9, 0.08), xycoords='axes fraction', ha='right', va='bottom', fontsize='small')
	#ax.annotate(text='ECMWF, 1985-2020', xy=(0.9, 0.02), xycoords='axes fraction', ha='right', va='bottom', fontsize='small')
	fig.savefig('../figures/map_stations_{0:s}.png'.format(label), bbox_inches='tight', dpi=400)
