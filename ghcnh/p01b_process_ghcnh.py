#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import getpass

import numpy as np
import pandas as pd

import datetime

from local_paths import *

## ============================================================== ##

YEARS = range(2011, 2020+1, 1)

## ============================================================== ##

for year in YEARS:
	print(year)
	ipath = os.path.join(DATAPATH_GHCNH_WITH_META, str(year))
	all_files = os.listdir(ipath)

	# merge files with station meta data
	ifiles = [f for f in all_files if '_meta.csv' in f]
	files_list = []
	for ifile in ifiles:
		df = pd.read_csv(os.path.join(ipath, ifile))
		files_list.append(df)
	df_all = pd.concat(files_list, axis=0, ignore_index=True)
	df_all = df_all.sort_values(by='Station_ID', ascending=True).reset_index(drop=True)
	df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_meta.csv'.format(year)), index=False)

	# merge files with 00 utc
	ifiles = [f for f in all_files if '_t2m_utc00.csv' in f]
	cols = ['date', 'temperature']
	files_list = []
	for ifile in ifiles:
		df = pd.read_csv(os.path.join(ipath, ifile))
		df = df.loc[:, cols]
		df['date'] = pd.to_datetime(df['date']).apply(lambda x: datetime.datetime.strftime(x, "%Y-%m-%d"))
		df['Station_ID'] = ifile.split('_')[1]
		files_list.append(df)
	df_all = pd.concat(files_list, axis=0, ignore_index=True)
	df_all = df_all.sort_values(by=['Station_ID', 'date'], ascending=True).reset_index(drop=True)
	df_all = df_all.loc[:, ['Station_ID', 'date', 'temperature']]
	df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_t2m_utc00.csv'.format(year)), index=False)
 
	# merge files with 12 utc
	ifiles = [f for f in all_files if '_t2m_utc12.csv' in f]
	cols = ['date', 'temperature']
	files_list = []
	for ifile in ifiles:
		df = pd.read_csv(os.path.join(ipath, ifile))
		df = df.loc[:, cols]
		df['date'] = pd.to_datetime(df['date']).apply(lambda x: datetime.datetime.strftime(x, "%Y-%m-%d"))
		df['Station_ID'] = ifile.split('_')[1]
		files_list.append(df)
	df_all = pd.concat(files_list, axis=0, ignore_index=True)
	df_all = df_all.sort_values(by=['Station_ID', 'date'], ascending=True).reset_index(drop=True)
	df_all = df_all.loc[:, ['Station_ID', 'date', 'temperature']]
	df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_t2m_utc12.csv'.format(year)), index=False)
 
	# merge files with daily means (first, average the 24 hourly values)
	ifiles = [f for f in all_files if '_t2m_daily.csv' in f]
	cols = ['date', 'temperature']
	files_list = []
	for ifile in ifiles:
		df = pd.read_csv(os.path.join(ipath, ifile))
		df = df.loc[:, cols]
		df['date'] = pd.to_datetime(df['date']).apply(lambda x: datetime.datetime.strftime(x, "%Y-%m-%d"))
		df = df.groupby('date')['temperature'].mean().reset_index()
		df['Station_ID'] = ifile.split('_')[1]
		files_list.append(df)
	df_all = pd.concat(files_list, axis=0, ignore_index=True)
	df_all = df_all.sort_values(by=['Station_ID', 'date'], ascending=True).reset_index(drop=True)
	df_all = df_all.loc[:, ['Station_ID', 'date', 'temperature']]
	df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_t2m_daily.csv'.format(year)), index=False)

	"""

	# merge files with min and max
	ifiles = [f for f in all_files if '_t2m_daily.csv' in f]
	cols = ['date', 'temperature']
	files_list = []
	for ifile in ifiles:
		df = pd.read_csv(os.path.join(ipath, ifile))
		df = df.loc[:, cols]
		df['date'] = pd.to_datetime(df['date']).apply(lambda x: datetime.datetime.strftime(x, "%Y-%m-%d"))
		df = (
			df.groupby('date')['temperature']
				.agg(['min', 'max'])
				.reset_index()
				.rename(columns={'min': 'temperature_min', 'max': 'temperature_max'})
		)
		df['Station_ID'] = ifile.split('_')[1]
		files_list.append(df)
	df_all = pd.concat(files_list, axis=0, ignore_index=True)
	df_all = df_all.sort_values(by=['Station_ID', 'date'], ascending=True).reset_index(drop=True)
	df_all = df_all.loc[:, ['Station_ID', 'date', 'temperature_min', 'temperature_max']]
	df_all.to_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, 'GHCNh_{0:d}_t2m_minmax.csv'.format(year)), index=False)

	"""