#!/usr/bin/env python

'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''

import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import  sys

var = 'Vy'
Uvar = 'Ux'
Keys = ['11', '12', '13']

Keys = [ '26', '29', '32', '62', '65', '68', '98', '101', '104' ]


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )


DT = params["DT"]
NT = int( params["TMAX"] / params["DT"] ) - 1


t = np.linspace( 0, params["TMAX"], NT )
stationFile = open( "station.json" )
stationDir = json.load( stationFile )

station = stationDir["station(point)"] 

fileName = "%s/station" % params["out"]

WAVESIZE = 9

mpiX = -1
mpiY = -1
mpiZ = -1

PX = grid.PX
PY = grid.PY
PZ = grid.PZ

stationNum = np.zeros( [PZ, PY, PX], dtype = 'int32' ) 

stationData = { }

seq = 0



varDir = { "Vx" : 0, "Vy" : 1, "Vz" : 2, "Txx" : 3, "Tyy" : 4, "Tzz" : 5, "Txy" : 6, "Txz" : 7, "Tyz" : 8 }

varId = varDir[var]



for index in station.values( ):
	X = index[0]
	Y = index[1]
	Z = index[2]

	for mpiZ in range( grid.PZ ):
		for mpiY in range( grid.PY ):
			for mpiX in range( grid.PX ):
				thisX = X - grid.frontNX[mpiX] + grid.halo
				thisY = Y - grid.frontNY[mpiY] + grid.halo
				thisZ = Z - grid.frontNZ[mpiZ] + grid.halo
				if thisX >= grid.halo and thisX < grid._nx[mpiX] and thisY >= grid.halo and thisY < grid._ny[mpiY] and thisZ >= grid.halo and thisZ < grid._nz[mpiZ]:
					stationNum[mpiZ, mpiY, mpiX] += 1

stationForDrawNum = len( Keys )
U  = np.zeros( [ stationForDrawNum, NT ] )
Ux = np.zeros( [ stationForDrawNum, NT ] )
Uy = np.zeros( [ stationForDrawNum, NT ] )
Uz = np.zeros( [ stationForDrawNum, NT ] )


stationKeyNum = { }



num = 0
for mpiZ in range( grid.PZ ):
	for mpiY in range( grid.PY ):
		for mpiX in range( grid.PX ):
			if stationNum[mpiZ, mpiY, mpiX] != 0:
				FileName = "%s_mpi_%d_%d_%d.bin" % ( fileName, mpiX, mpiY, mpiZ )
				File = open( FileName, "rb" )
				print( FileName )
				count = stationNum[mpiZ, mpiY, mpiX]
				XIndex = np.fromfile( File, dtype='int32', count = count )
				YIndex = np.fromfile( File, dtype='int32', count = count )
				ZIndex = np.fromfile( File, dtype='int32', count = count )

				count = NT * stationNum[mpiZ, mpiY, mpiX] * WAVESIZE
				data = np.fromfile( File, dtype='float32', count = count )
				dataRe = np.reshape( data, ( WAVESIZE, stationNum[mpiZ, mpiY, mpiX], NT ) )
				stationData[( mpiX, mpiY, mpiZ )] = dataRe
				for key in Keys:
					xidx = station[key][0]
					yidx = station[key][1]
					zidx = station[key][2]
					print( "key = %s, X = %d, Y = %d, Z = %d" % ( key, xidx, yidx, zidx ) )
					for i in range( stationNum[mpiZ, mpiY, mpiX] ):
						Ux_ = np.zeros( NT )
						Uy_ = np.zeros( NT )
						Uz_ = np.zeros( NT )

						'''
						UxSum = 0.0
						UySum = 0.0
						UzSum = 0.0
						for it in range( 1, NT ):
							UxSum += ( dataRe[0, i, it] + dataRe[0, i, it] ) * DT * 0.5
							UySum += ( dataRe[1, i, it] + dataRe[1, i, it] ) * DT * 0.5
							UzSum += ( dataRe[2, i, it] + dataRe[2, i, it] ) * DT * 0.5
							Ux_[it] = UxSum
							Uy_[it] = UySum
							Uz_[it] = UzSum
						'''

						if xidx == XIndex[i] and yidx == YIndex[i] and zidx == ZIndex[i]:
							Ux[num] = Ux_
							Uy[num] = Uy_
							Uz[num] = Uz_

							stationKeyNum[key] = num
							num += 1


if Uvar == 'Ux':
	U = Ux
else:
	if Uvar == 'Uy':
		U = Uy
	else:
		if Uvar == 'Uz':
			U = Uz	
		else:
			U = U + 1
vmax = np.max( np.abs( U ) )



plt.figure( 1 )
for key in station.keys( ):
	#print( key )
	for iKey in Keys:
		if key == iKey:
			print( key )
			i = stationKeyNum[key]
			plt.plot( t, U[i]/vmax + i , color = 'r', ls = '--' )
			#axs[i].plot( t[i], Ux[i] )
			break

plt.show( )
np.save( "GRTM3DDir/cgfd3d_Ux.npy", Ux )
np.save( "GRTM3DDir/cgfd3d_Uy.npy", Uy )
np.save( "GRTM3DDir/cgfd3d_Uz.npy", Uz )


#plt.plot( t, dataRe[varId, i, :] / np.max( dataRe[varId, i, :] ), '-' )
#plt.plot( t, U, '-' )

