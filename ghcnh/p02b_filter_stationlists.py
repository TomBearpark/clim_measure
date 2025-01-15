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
NDAYS_MIN = 30

## =============================== ##

# Load the data
ifile = 'GHCNh_{0:d}-{1:d}_stationlist_utc00.csv'.format(YEARS[0], YEARS[-1])
df = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
df = df.loc[df['year'].between(YEARS[0], YEARS[-1]), :]
df = df.loc[df['nobs_temperature'] >= NDAYS_MIN, :]
stations = df.groupby('Station_ID')['year'].count()
stations = stations[stations == (YEARS[-1] - YEARS[0] + 1)].index.values
df = df.loc[df['Station_ID'].isin(stations), :].groupby('Station_ID').first().reset_index().drop(columns=['year'])
df = df.sort_values(by='Station_ID', ascending=True).reset_index(drop=True)
ofile = 'GHCNh_stationlist_utc00_{0:d}-{1:d}_{2:03d}plus.csv'.format(YEARS[0], YEARS[-1], NDAYS_MIN)
df.to_csv(os.path.join(DATAPATH_STATIONLISTS, ofile), index=False)

# Load the data
ifile = 'GHCNh_{0:d}-{1:d}_stationlist_utc12.csv'.format(YEARS[0], YEARS[-1])
df = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
df = df.loc[df['year'].between(YEARS[0], YEARS[-1]), :]
df = df.loc[df['nobs_temperature'] >= NDAYS_MIN, :]
stations = df.groupby('Station_ID')['year'].count()
stations = stations[stations == (YEARS[-1] - YEARS[0] + 1)].index.values
df = df.loc[df['Station_ID'].isin(stations), :].groupby('Station_ID').first().reset_index().drop(columns=['year'])
df = df.sort_values(by='Station_ID', ascending=True).reset_index(drop=True)
ofile = 'GHCNh_stationlist_utc12_{0:d}-{1:d}_{2:03d}plus.csv'.format(YEARS[0], YEARS[-1], NDAYS_MIN)
df.to_csv(os.path.join(DATAPATH_STATIONLISTS, ofile), index=False)

# Load the data
ifile = 'GHCNh_{0:d}-{1:d}_stationlist_daily.csv'.format(YEARS[0], YEARS[-1])
df = pd.read_csv(os.path.join(DATAPATH_GHCNH_PROCESSED, ifile))
df = df.loc[df['year'].between(YEARS[0], YEARS[-1]), :]
df = df.loc[df['nobs_temperature'] >= NDAYS_MIN, :]
stations = df.groupby('Station_ID')['year'].count()
stations = stations[stations == (YEARS[-1] - YEARS[0] + 1)].index.values
df = df.loc[df['Station_ID'].isin(stations), :].groupby('Station_ID').first().reset_index().drop(columns=['year'])
df = df.sort_values(by='Station_ID', ascending=True).reset_index(drop=True)
ofile = 'GHCNh_stationlist_daymean_{0:d}-{1:d}_{2:03d}plus.csv'.format(YEARS[0], YEARS[-1], NDAYS_MIN)
df.to_csv(os.path.join(DATAPATH_STATIONLISTS, ofile), index=False)

## =============================== ##
