#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getpass

username = getpass.getuser()
if username == 'manuel':
	DATAPATH_GHCNH_WITH_META = '/media/manuel/Extreme SSD/GHCNh/' # about 15 GB
	DATAPATH_GHCNH_PROCESSED = '/home/manuel/Research/projects/project_62/data/ghcnh' # about 2.5 GB
	DATAPATH_STATIONLISTS = '../data/' # < 1 GB
	DATAPATH_DISTANCES = '../data/' # < 1 GB
	DATAPATH_INTERPOLATED = '/home/manuel/Research/projects/project_62/data/ghcnh_interpolated/' # about 2.5 GB
elif username == 'yourname':
	DATAPATH_GHCNH_WITH_META = 'your/directory/'
	DATAPATH_GHCNH_PROCESSED = 'your/directory/'
	DATAPATH_STATIONLISTS = 'your/directory/'
	DATAPATH_DISTANCES = 'your/directory/'
	DATAPATH_INTERPOLATED = 'your/directory/'
else:
	exit('Please edit the file local_paths and specify the local directories.')
