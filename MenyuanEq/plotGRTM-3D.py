#!/usr/bin/env python

import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import os
import json

#Keys = ['9', '10', '11', '12', '14', '15', '16']
#Keys = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
#		'10', '11', '12', '14', '15', '16', '17', '18', '19']
Keys = [ '10', '11', '12', '13', '14', '15']
Uvar = 'Ur'

GRTM3DDir = "GRTM3DDir"
GRTM3D = GRTM3DDir + "/grtcsgram"
file = open( GRTM3DDir + "/input.conf", 'r'  )
confStr = file.readlines(  )
file.close( )

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

DT = params["DT"]
DH = params["DH"]
TMAX = params["TMAX"]

sourceX = params['sourceX']
sourceY = params['sourceY']
sourceZ = params['sourceZ']

centerX = params['centerX']
centerY = params['centerY']


NT = int( params["TMAX"] / params["DT"] )
t = np.linspace( 0, params["TMAX"], NT )
rickerfc = params['rickerfc']

stationFile = open( "station.json" )
stationDir = json.load( stationFile )

station = stationDir["station(point)"] 

stationKeys     = station.keys( )
stationIndexes  = station.values( )

for key,value in  station.items( ):
	recvX = ( value[0] - sourceX ) * DH 
	recvY = ( value[1] - sourceY ) * DH 
	recvZ = ( grid.NZ - 1 - value[2] ) * DH 
	
	srcX = ( sourceX - centerX ) * DH 
	srcY = ( sourceY - centerY ) * DH 
	srcZ = ( grid.NZ - 1 - sourceZ ) * DH 

	recvList = confStr[20].split(  )
	
	recvList[2] = "%.2e," % recvX
	recvList[3] = "%.2e," % recvY
	recvList[4] = "%.2e"  % recvZ
	
	recvCat = ''
	recvCat += recvList[0]
	for iStr in recvList[1:]:
		recvCat += " " + iStr
	confStr[20] = recvCat + '\n'

	srcList = confStr[19].split(  )
	
	srcList[2] = "%.2e," % srcX
	srcList[3] = "%.2e," % srcY
	srcList[4] = "%.2e"  % srcZ

	srcCat = ''
	srcCat += srcList[0]
	for iStr in srcList[1:]:
		srcCat += " " + iStr
	confStr[19] = srcCat + '\n'
	
	confStr[7]  = 'model_file = \'%s/model.dat\'\n' % GRTM3DDir
	confStr[8]  = 'output_prefix = \'%s/outputDir/out%s\'\n' % ( GRTM3DDir, key )
	confStr[10] = 'record_time_length = %f\n' % TMAX
	confStr[11] = 'record_time_step = %e\n' % DT

	inputConf = "%s/inputDir/input%s.conf" % ( GRTM3DDir,  key )
	fileConf = open( inputConf, 'w' )
	for iStr in confStr:
		fileConf.write( iStr )
	fileConf.close( )
	os.system( "%s %s" % ( GRTM3D, inputConf )  )


stationNum = len( Keys )
t  = np.zeros( [ stationNum, NT ] )
U  = np.zeros( [ stationNum, NT ] )
Ur = np.zeros( [ stationNum, NT ] )
Ut = np.zeros( [ stationNum, NT ] )
Uz = np.zeros( [ stationNum, NT ] )


stationKeyNum = { }
i = 0
for key in Keys:
	stationKeyNum[key] = i
	Uout = '%s/outputDir/out%s.gu' % ( GRTM3DDir, key )
	txt = np.loadtxt( Uout, skiprows = 1)
	t [i] = txt[1:,0]
	Ur[i] = txt[1:,1]
	Ut[i] = txt[1:,2]
	Uz[i] = txt[1:,3]
	i += 1
	

#fit, axs = plt.subplots( stationNum, 1 )
plt.figure( 2 )

Ux = np.zeros( [stationNum, NT ] )
Uy = np.zeros( [stationNum, NT ] )

Ux = Ur
Uy = -Ut
Uz = -Uz

if Uvar == 'Ur':
	U = Ur
else:
	if Uvar == 'Ut':
		U = -Ut
	else:
		if Uvar == 'Uz':
			U = - Uz	
		else:
			U = U + 1


np.save( "CompareData/GRTM3D/GRTMUx.npy", Ux )
np.save( "CompareData/GRTM3D/GRTMUy.npy", Uy )
np.save( "CompareData/GRTM3D/GRTMUz.npy", Uz )


vmax = np.max( np.abs( U ) )
for key in station.keys( ):
	print( key )
	for iKey in Keys:
		if key == iKey:
			i = stationKeyNum[key]

			plt.plot( t[i], U[i] / vmax + i, color = 'r', ls = '--' )
			#axs[i].plot( t[i], Ur[i] )
			break

plt.show( )

