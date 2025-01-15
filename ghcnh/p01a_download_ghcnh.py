#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy

import numpy as np
import pandas as pd

import tarfile
import requests

from local_paths import *

## ============================================================== ##

YEARS = range(2011, 2020+1, 1)

## ============================================================== ##

# Define a function for processing
def process_archive(url, opath):

	# Step 1: Download the .tar.gz file
	filename = url.split('/')[-1]
	response = requests.get(url)
	with open(filename, 'wb') as file:
		file.write(response.content)
	
	# Step 2: Extract the file
	with tarfile.open(filename, 'r:gz') as tar:
		tar.extractall()
		extracted_files = tar.getnames()  # List of files inside the tar.gz

	# Step 3: Process each file (assuming CSV for this example)
	for extracted_file in extracted_files:
		print(extracted_file)
		if extracted_file.endswith('.psv'):  # Modify if other formats
			df = pd.read_csv(extracted_file, sep='|')
			# Example manipulation: simple summary of the data

			df_meta = df.loc[[0], ['Station_ID', 'Station_name', 'Latitude', 'Longitude', 'Elevation']]
			df_meta['Elevation'] = df_meta['Elevation'].replace(['******'], np.nan)
			df_meta['Elevation'] = df_meta['Elevation'].astype(float)
			df_meta.loc[df_meta['Elevation'] < -900, 'Elevation'] = np.nan
			df_meta.to_csv(os.path.join(opath, extracted_file.replace('.psv', '_meta.csv')), index=False)

			df['date'] = pd.to_datetime(df['DATE'])
			cols = [\
				'date',
				'temperature',
				'temperature_Quality_Code',
				'temperature_Report_Type',
				'temperature_Source_Code',
				#'temperature_Measurement_Code',
				#'temperature_Source_Station_ID',
				#'remarks',
				#'remarks_Measurement_Code',
				#'remarks_Quality_Code',
				#'remarks_Report_Type',
				#'remarks_Source_Code',
				#'remarks_Source_Station_ID'
				]
			df['temperature'] = df['temperature'].astype(float)
			df.loc[df['temperature'] < -900, 'temperature'] = np.nan
			df = df.loc[df['temperature'].notnull(), cols]
			df = df.loc[df['date'].dt.minute == 0, :]
			df = df.groupby('date').first().reset_index()
			df_t2m_utc12 = df.loc[df['date'].dt.hour == 12, :]
			df_t2m_utc00 = df.loc[df['date'].dt.hour == 00, :]
			df_t2m_utc12.to_csv(os.path.join(opath, extracted_file.replace('.psv', '_t2m_utc12.csv')), index=False)
			df_t2m_utc00.to_csv(os.path.join(opath, extracted_file.replace('.psv', '_t2m_utc00.csv')), index=False)

			df['date_str'] = df['date'].dt.date
			df['nhours'] = df.groupby('date_str')['temperature'].transform(lambda x: x.count())
			df = df.loc[df['nhours'] == 24, cols]
			df.to_csv(os.path.join(opath, extracted_file.replace('.psv', '_t2m_daily.csv')), index=False)

	# Step 4: Clean up - delete extracted files and the archive
	for extracted_file in extracted_files:
		os.remove(extracted_file)
	os.remove(filename)

## ============================================================== ##

base_url = 'https://www.ncei.noaa.gov/oa/global-historical-climatology-network/hourly/archive/'

for year in YEARS:
	print(year)
	filename = 'ghcn-hourly_v1.0.0_d{0:04d}_c20240709.tar.gz'.format(year)
	opath = os.path.join(DATAPATH_GHCNH_WITH_META, str(year))
	if not os.path.exists(opath):
		os.makedirs(opath, exist_ok=True)
	process_archive(base_url + filename, opath)
