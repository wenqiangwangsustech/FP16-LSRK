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

import os


StationLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 
				 'H', 'I', 'J', 'K', 'L', 'M', 'N',
				 'O', 'P', 'Q',		 'R', 'S', 'T',
				 'U', 'V', 'W', 	 'X', 'Y', 'Z']

var  = 'Vx' #Vy Vz
Uvar = 'Ux'
#Keys = [ '9', '11','15', '16' ]
#Keys = [ '9', '10', '11', '12']# '13', '14', '15', '16', '17' ]
#Keys = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17' ]
#Keys = [ '9', '10', '11', '12', '13', '14', '15', '16', '17' ]
#Keys = [ '10', '11', '12', '13', '14', '15']
Keys = [ '9', '11', '13', '15', '17' ]

#StationLabels = ['E'] 
#Keys = [ '17']

if len( Keys ) == 1:
	ampStr = "amplify"
else:
	ampStr = ''


StationTicks = np.arange( 0, len(Keys) )
lineWidth = 1
IsDiffer = 1


colors = [ 'k', 'r', 'g', 'b', 'y', 'm']
#lines  = [ '-', '--', '-.', ':', '*', ',']

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )


DT = params["DT"]

TMAX = params["TMAX"]

NT = int( TMAX / DT ) - 1


t = np.linspace( 0, TMAX, NT )
stationFile = open( "station.json" )
stationDir = json.load( stationFile )

station = stationDir["station(point)"] 


zoom = 10


NameMode = [ ]
NameMode.append( "CGFDM3D" )#params["out"]
NameMode.append( "CGFDM3D-LSRK" )#params["out"]
#NameMode.append( "CGFDM3D-CJM" )#params["out"]
NameMode.append( "CGFDM3D-CJMVS" )#params["out"]
#NameMode.append( "CGFDM3D-LSRK-CJM" )#params["out"]
NameMode.append( "CGFDM3D-LSRK-CJMVS" )#params["out"]







NVersion = len( NameMode )

dataPrefix = "CompareData/"
GRTM3DDir = dataPrefix + "GRTM3D/"


prefix = ""
fileName = [ ] 

for i in range( NVersion ):
	fileName.append( dataPrefix + NameMode[i] + "/station" ) 


NameMode = [ ]
NameMode.append( "CGFDM3D" )#params["out"]
NameMode.append( "LSRK" )#params["out"]
#NameMode.append( "CJM" )#params["out"]
NameMode.append( "CJMVS" )#params["out"]
#NameMode.append( "LSRK-CJM" )#params["out"]
NameMode.append( "LSRK-CJMVS" )#params["out"]

#fileName[0] = 'output/station'


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
U  = np.zeros( [ NVersion, stationForDrawNum, NT ] )
Ux = np.zeros( [ NVersion, stationForDrawNum, NT ] )
Uy = np.zeros( [ NVersion, stationForDrawNum, NT ] )
Uz = np.zeros( [ NVersion, stationForDrawNum, NT ] )

stationKeyNum = { }


for version in range( NVersion ):
	num = 0
	for mpiZ in range( grid.PZ ):
		for mpiY in range( grid.PY ):
			for mpiX in range( grid.PX ):
				if stationNum[mpiZ, mpiY, mpiX] != 0:
					FileName = "%s_mpi_%d_%d_%d.bin" % ( fileName[version], mpiX, mpiY, mpiZ )
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
						#print( "key = %s, X = %d, Y = %d, Z = %d" % ( key, xidx, yidx, zidx ) )
						for i in range( stationNum[mpiZ, mpiY, mpiX] ):
							Ux_ = np.zeros( NT )
							Uy_ = np.zeros( NT )
							Uz_ = np.zeros( NT )
	
							UxSum = 0.0
							UySum = 0.0
							UzSum = 0.0
							'''
							for it in range( 1, NT ):
								UxSum += ( dataRe[0, i, it - 1] + dataRe[0, i, it]) * DT * 0.5
								UySum += ( dataRe[1, i, it - 1] + dataRe[1, i, it]) * DT * 0.5
								UzSum += ( dataRe[2, i, it - 1] + dataRe[2, i, it]) * DT * 0.5
								Ux_[it] = UxSum
								Uy_[it] = UySum
								Uz_[it] = UzSum
							'''
							Ux_[:] = dataRe[0, i, :]
							Uy_[:] = dataRe[1, i, :]
							Uz_[:] = dataRe[2, i, :]

							if xidx == XIndex[i] and yidx == YIndex[i] and zidx == ZIndex[i]:
								print( np.shape( Ux ) )
								print( np.shape( Ux_ ) )
								print( "key = %s, X = %d, Y = %d, Z = %d" % ( key, xidx, yidx, zidx ) )
								for it in range( NT ):
									Ux[version, num, it] = Ux_[it]
									Uy[version, num, it] = Uy_[it]
									Uz[version, num, it] = Uz_[it]
	
								stationKeyNum[key] = num
								num += 1
	
if Uvar == 'Ux':
	U = Ux
else:
	if Uvar == 'Uy':
		U = Uy #/ np.sqrt( np.pi )# 2. * np.pi * np.sqrt( np.pi )
	else:
		if Uvar == 'Uz':
			U = Uz #/ np.sqrt( np.pi )#2. * np.pi *  np.sqrt( np.pi )	
		else:
			U = U + 1

vmax = np.max( np.abs( U ) )


if len( Keys ) > 1:
	vmaxUx = np.max( np.abs( Ux ) )
	vmaxUy = np.max( np.abs( Uy ) )
	vmaxUz = np.max( np.abs( Uz ) )
	vmax = np.max( [vmaxUx, vmaxUy, vmaxUz] )

print( "CGFDM3D Max = %f" % vmax  )
print( stationKeyNum )



if len( Keys ) > 1:
	plt.figure( figsize = ( 8,10) )
else:
	plt.figure( figsize = ( 8,2) )



for version in range( NVersion ):
	for key in station.keys( ):
		#print( key )
		for iKey in Keys:
			if key == iKey:
				print( key )
				i = stationKeyNum[key]
				if i == 0:
					plt.plot( t, U[version, i] / vmax + i, color = colors[version], ls = '-', label = NameMode[version], linewidth = lineWidth )
					if version == 0:
						plt.text( TMAX * 0.85, i + 0.05, "%.2fcm/s" % np.max( np.abs( U[version, i] * 100 ) ) )
					#if version != 0:
					#	plt.plot( t, (U[version, i] - U[0, i] ) * zoom / vmax + i, color = colors[version], ls = '--', label = NameMode[version] + ':Error $ \\times $ %d' % zoom , linewidth = lineWidth)

					if len( Keys ) > 1:
						plt.legend( loc = 1 )
				else:
					plt.plot( t, U[version, i] / vmax + i, color = colors[version], ls = '-', linewidth = lineWidth )
					if version == 0:
						plt.text( TMAX * 0.85, i + 0.05, "%.2fcm/s" % np.max( np.abs( U[version, i] * 100 ) ) )
					#if version != 0:
					#	plt.plot( t, ( U[version, i] - U[0, i] ) * zoom / vmax + i, color = colors[version], ls = '--',linewidth = lineWidth )
				#axis[i].plot( t[i], Ux[i] )
				break




GRTM3D_Ux = np.load( GRTM3DDir + "GRTMVx.npy" )
GRTM3D_Uy = np.load( GRTM3DDir + "GRTMVy.npy" )
GRTM3D_Uz = np.load( GRTM3DDir + "GRTMVz.npy" )
GRTM3D_U  = GRTM3D_Ux


if Uvar == 'Ux':
	GRTM3D_U = GRTM3D_Ux# * 2 * np.pi
else:
	if Uvar == 'Uy':
		GRTM3D_U = GRTM3D_Uy# * 2 * np.pi #/ ( 2 * np.sqrt( np.pi ) )
	else:
		if Uvar == 'Uz':
			GRTM3D_U = GRTM3D_Uz# * 2 * np.pi #/ ( 2 * np.sqrt( np.pi ) )

grtmvmax = np.max( np.abs( GRTM3D_U ) )

if len( Keys ) > 1:
	grtmvmaxUx = np.max( np.abs( GRTM3D_Ux ) )
	grtmvmaxUy = np.max( np.abs( GRTM3D_Uy ) )
	grtmvmaxUz = np.max( np.abs( GRTM3D_Uz ) )
	grtmvmax = np.max( [grtmvmaxUx, grtmvmaxUy, grtmvmaxUz] )



print( "GRTM Max = %f" % grtmvmax  )

T = np.linspace( 0, TMAX, NT + 1 )

print( np.shape( GRTM3D_U ) )
for iKey in range( stationForDrawNum ):
	print( iKey )
	if iKey == 0:
		plt.plot( T[:], GRTM3D_U[iKey, :] / grtmvmax + iKey, color = "c", ls = '--', label = "GRTM", linewidth = lineWidth )
		if len( Keys ) > 1:
			plt.legend( loc = 1 )
	else:
		plt.plot( T[:], GRTM3D_U[iKey, :] / grtmvmax + iKey, color = "c", ls = '--', linewidth = lineWidth )

if len( Keys ) == 1:
	plt.gca( ).set_ylim( [-1, StationTicks[-1] + 1] )
else:
	plt.gca( ).set_ylim( [-1, StationTicks[-1] + 1.2] )

plt.gca( ).set_xlim( [0, TMAX] )
plt.yticks( StationTicks[:len(Keys)],  StationLabels[:len(Keys)])
#plt.title( "%s on the surface with %s model" % ( Uvar, Model ) )
if var == 'Vx':
	plt.title( "Layered Media Model: $v_x$" )
if var == 'Vy':
	plt.title( "Layered Media Model: $v_y$" )
if var == 'Vz':
	plt.title( "Layered Media Model: $v_z$" )
plt.xlabel( 't(s)'  )

#plt.show( )
plt.savefig( "compGRTM_%s_%s.pdf" % ( var, ampStr )  )
#plt.plot( t, dataRe[varId, i, :] / np.max( dataRe[varId, i, :] ), '-' )
#plt.plot( t, U, '-' )
os.system("evince compGRTM_%s_%s.pdf" % ( var, ampStr ) )
