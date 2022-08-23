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


stationDic = { }

file = open( 'SM_coor', 'r' )
LonLatTxt = file.readlines( )
file.close( )

stationNum = len( LonLatTxt )
stationName = []
Lon = np.zeros( stationNum  )
Lat = np.zeros( stationNum  )
for i in range( stationNum ):
	staStr = LonLatTxt[i].split( )
	stationName.append( staStr[0] )
	Lon[i] = float( staStr[1] )
	Lat[i] = float( staStr[2] )

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

DH = params["DH"]
centerX = params["centerX"]
centerY = params["centerY"]

nPML = params['nPML']
XRange1 = nPML + 1
XRange2 = params['NX'] - nPML

YRange1 = nPML + 1
YRange2 = params['NY'] - nPML


latC = params["centerLatitude"]
lonC = params["centerLongitude"]
proj = Proj( proj='aeqd', lat_0 = latC, lon_0 = lonC, ellps = "WGS84" )
NZ = params['NZ']

XCoord = np.zeros( stationNum )
YCoord = np.zeros( stationNum )
XCoord, YCoord = proj( Lon, Lat )

X = np.zeros( stationNum )
Y = np.zeros( stationNum )
n = 0
for i in range( stationNum ):
	X[i] = round( XCoord[i] / DH ) + centerX
	Y[i] = round( YCoord[i] / DH ) + centerY
	if X[i] >= XRange1 and X[i] < XRange2 and Y[i] >= YRange1 and Y[i] < YRange2:
		x = X[i]; y = Y[i]; z = NZ - 1;
		stationDic[stationName[i]] = [int(x),int(y),int(z)]



stationNum = 12
step = 100
offSet = 99

newStation = { }
for i in range( stationNum ):
	for j in range( stationNum ):
		index = j + i * stationNum
		xIndex = i * step + offSet
		yIndex = j * step + offSet
		print( xIndex, yIndex, NZ - 1 )
		newStation[index] = [ xIndex, yIndex,NZ - 1]


stationParam = {

	"point"   : 1,
	"lon_lat" : 0,

	"station(point)" : newStation
}

print( stationParam  )

stationJsonFile = open( "station.json", 'w' )

json_str = json.dumps( stationParam, indent = 4 )
stationJsonFile.write( json_str )
stationJsonFile.close( )


#json_str = json.dumps( stationParam )
#print( json_str )


