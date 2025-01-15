#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import numpy as np
import pandas as pd

from local_paths import *

## ============================================================== ##

YEARS = list(range(2011, 2020+1, 1))

## ============================================================== ##

list_utc00 = []
list_utc12 = []
list_daily = []

for year in YEARS:
	print(year)

	ifile = "GHCNh_{0:04d}_meta.csv".format(year)
	df_meta = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))

	ifile = "GHCNh_{0:04d}_t2m_utc00.csv".format(year)
	df_nobs = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
	df_nobs = df_nobs.groupby('Station_ID')['temperature'].count().reset_index().rename(columns={'temperature': 'nobs_temperature'})
	df = df_meta.merge(df_nobs, on='Station_ID', how='left')
	df.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_stationlist_utc00.csv'.format(year)), index=False)
	df['year'] = year
	list_utc00.append(df)

	ifile = "GHCNh_{0:04d}_t2m_utc12.csv".format(year)
	df_nobs = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
	df_nobs = df_nobs.groupby('Station_ID')['temperature'].count().reset_index().rename(columns={'temperature': 'nobs_temperature'})
	df = df_meta.merge(df_nobs, on='Station_ID', how='left')
	df.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_stationlist_utc12.csv'.format(year)), index=False)
	df['year'] = year
	list_utc12.append(df)

	ifile = "GHCNh_{0:04d}_t2m_daily.csv".format(year)
	df_nobs = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
	df_nobs = df_nobs.groupby('Station_ID')['temperature'].count().reset_index().rename(columns={'temperature': 'nobs_temperature'})
	df = df_meta.merge(df_nobs, on='Station_ID', how='left')
	df.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_stationlist_daily.csv'.format(year)), index=False)
	df['year'] = year
	list_daily.append(df)

df_all = pd.concat(list_utc00, axis=0, ignore_index=True)
df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}-{1:d}_stationlist_utc00.csv'.format(YEARS[0], YEARS[-1])), index=False)

df_all = pd.concat(list_utc12, axis=0, ignore_index=True)
df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}-{1:d}_stationlist_utc12.csv'.format(YEARS[0], YEARS[-1])), index=False)

df_all = pd.concat(list_daily, axis=0, ignore_index=True)
df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}-{1:d}_stationlist_daily.csv'.format(YEARS[0], YEARS[-1])), index=False)
