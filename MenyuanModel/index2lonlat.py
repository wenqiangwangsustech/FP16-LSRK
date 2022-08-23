#!/usr/bin/env python

"""
Created on 2022-02-24
23:59
@author: Wenqiang Wang 
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import struct
import sys
import os
from pyproj import Proj
import json
from mpl_toolkits.mplot3d import Axes3D
from pyscripts.GRID import GRID




Keys = [ '26', '29', '32', '62', '65', '68', '98', '101', '104' ]

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

DH = params["DH"]
centerX = params["centerX"]
centerY = params["centerY"]


latC = params["centerLatitude"]
lonC = params["centerLongitude"]
proj1 = Proj( proj='aeqd', lat_0 = latC, lon_0 = lonC )
NZ = params['NZ']


stationFile = open( "station.json"  )
stationDir = json.load( stationFile  )

stations = stationDir["station(point)"]


for key in Keys:
	X = stations[key][0]
	Y = stations[key][1]
	#print( X, Y )
	x = ( X - centerX ) * DH
	y = ( Y - centerY ) * DH
	
	lon, lat = proj1( x, y, inverse = "true" )
	print( "%f, %f" % ( lon, lat ) )
	
	


