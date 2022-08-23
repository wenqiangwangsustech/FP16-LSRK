#!/usr/bin/env python

'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import sys


it = 2500

StationLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 
				 'H', 'I', 'J', 'K', 'L', 'M', 'N',
				 'O', 'P', 'Q',		 'R', 'S', 'T',
				 'U', 'V', 'W', 	 'X', 'Y', 'Z']

Keys = [ '7', '9', '11', '13', '15' ]# '13', '14', '15', '16', '17' ]
Keys = [ '19', '20', '21', '22' ]
Keys = [ '26', '29', '32', '62', '65', '68', '98', '101', '104' ]
Keys = [ '25', '27', '29', '31', '61', '63', '65', '67', '97', '99', '101', '103' ]

Keys = [ '14', '16', '18', '20',
		 '50', '52', '54', '56',
		 '74', '76', '78', '80',
		 '98', '100', '102', '104']
Keys = [ '26', '29', '32', '62', '65', '68', '98', '101', '104' ]

colors = [ 'k', 'r', 'g', 'b', 'y', 'm']

var = 'Vx' #Vx Vz
name = 0

DrawDiff = 0

FreeSurf =1

if len( sys.argv ) > 1:
	it = int( sys.argv[1] )
	var = str( sys.argv[2] )
	name = int( sys.argv[3] )
	DrawDiff = int( sys.argv[4] )



jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 5

outputPath = params["out"]
if name == 0:
	outputPath = 'CompareData/CGFDM3D' 
if name == 1:
	outputPath = 'CompareData/CGFDM3D-LSRK' 
if name == 2:
	outputPath = 'CompareData/CGFDM3D-CJMVS' 
if name == 3:
	outputPath = 'CompareData/CGFDM3D-LSRK-CJMVS' 

fileNameX = params["out"] + "/coordX"
fileNameY = params["out"] + "/coordY"

DT = params["DT"]

if FreeSurf == 1:
	fileNameTh = params["out"] + "/FreeSurf%s_%d"%( var, it )
	fileName = outputPath + "/FreeSurf%s_%d" % ( var, it )
else:
	fileNameTh = params["out"] + "/%s_%d"%( var, it )
	fileName = outputPath + "/%s_%d" % ( var, it )

print( "Draw " + fileName  )


FAST_AXIS = params["FAST_AXIS"]


'''
for Z in range( grid.PZ ):
	for Y in range( grid.PY ):
		for X in range( grid.PX ):
			print( "nx = %d, ny = %d, nz = %d\n" % ( grid.nx[X], grid.ny[Y], grid.nz[Z] ) )
'''

sliceX = params["sliceX"] - grid.frontNX
sliceY = params["sliceY"] - grid.frontNY
sliceZ = params["sliceZ"] - grid.frontNZ

for mpiSliceX in range( grid.PX ):
	if sliceX[mpiSliceX] >= 0 and sliceX[mpiSliceX] < grid.nx[mpiSliceX]:
		break

for mpiSliceY in range( grid.PY ):
	if sliceY[mpiSliceY] >= 0 and sliceY[mpiSliceY] < grid.ny[mpiSliceY]:
		break

for mpiSliceZ in range( grid.PZ ):
	if sliceZ[mpiSliceZ] >= 0 and sliceZ[mpiSliceZ] < grid.nz[mpiSliceZ]:
		break



if FAST_AXIS == 'Z':
	dataX = np.zeros( [grid.NX, grid.NY] )
	dataY = np.zeros( [grid.NX, grid.NY] )
	data  = np.zeros( [grid.NX, grid.NY] )
	dataTh= np.zeros( [grid.NX, grid.NY] )
else:
	dataX = np.zeros( [grid.NY, grid.NX] )
	dataY = np.zeros( [grid.NY, grid.NX] )
	data  = np.zeros( [grid.NY, grid.NX] )
	dataTh= np.zeros( [grid.NY, grid.NX] )



for mpiY in range( grid.PY ):
	for mpiX in range( grid.PX ):
		mpiZ = mpiSliceZ
		XFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		YFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		if FreeSurf:
			mpiZ = grid.PZ - 1
		else:
			mpiZ = mpiSliceZ
		File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		FileTh= open( "%s_Z_mpi_%d_%d_%d.bin" % (fileNameTh, mpiX, mpiY, mpiZ ), "rb" )
		ny = grid.ny[mpiY]
		nx = grid.nx[mpiX]

		print( "ny = %d, nx = %d" % ( nx, ny ) )
		datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
		datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
		data_ = np.fromfile(  File, dtype='float32', count = ny * nx )
		data_Th= np.fromfile(  FileTh, dtype='float32', count = ny * nx )

		J  = grid.frontNY[mpiY]
		J_ = grid.frontNY[mpiY] + ny
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx

		if FAST_AXIS == 'Z':
			dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
			dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
			data [I:I_, J:J_] = np.reshape( data_, (nx, ny) )
			dataTh[I:I_, J:J_] = np.reshape( data_Th, (nx, ny) )
		else:
			dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
			dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
			data [J:J_, I:I_] = np.reshape( data_, (ny, nx) )
			dataTh[J:J_, I:I_] = np.reshape( data_Th, (ny, nx) )


dpi = 5000
fig = plt.figure()
#plt.figure(figsize=(5, 5))
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
if DrawDiff:
	vm = np.max( np.abs( data - dataTh ) ) #* 1e-2
else:
	vm = np.max( np.abs( data ) ) #* 1e-2
	if it == 4000:
		vm = 5.8e-2
	

dataX = dataX/unit
dataY = dataY/unit





sourceX = params['sourceX']
sourceY = params['sourceY']
sourceZ = params['sourceZ']

centerX = params['centerX']
centerY = params['centerY']
DH = params['DH']

srcX = ( sourceX - centerX ) * DH 
srcY = ( sourceY - centerY ) * DH 
srcZ = - ( grid.NZ - 1 - sourceZ ) * DH 

#if DrawDiff:
#	plt.pcolormesh( dataX[::sample], dataY[::sample], np.abs( data[::sample] - dataTh [::sample]), vmax = vm, vmin = 0, cmap = "Reds" )
#else:                                           
#	plt.pcolormesh( dataX[::sample], dataY[::sample], data[::sample], vmax = vm, vmin = -vm, cmap = "seismic" )

if DrawDiff:
	plt.pcolor( dataX, dataY, np.abs( data - dataTh ), vmax = vm, vmin = 0, cmap = "Reds", rasterized = True )
else:                       
	plt.pcolor( dataX, dataY, data, vmax = vm, vmin = -vm, cmap = "seismic", rasterized = True )



#cb = plt.colorbar( orientation =  "horizontal", format = '%.1e', shrink = 0.5 )
cb = plt.colorbar( format = '%.1e', shrink = 0.8 )

plt.axis( "image" )

Xmin = np.min( dataX )
Xmax = np.max( dataX )
Ymin = np.min( dataY )
Ymax = np.max( dataY )
#xticksRange = np.linspace( -10.0, 10.0, 5)
# 
#zticksRange = [ -10, -8, -6, -4, -2, 0, 3 ] #np.linspace( -10, 2.5, 5 )
#
#plt.xticks(xticksRange)
#plt.yticks(zticksRange)




#plt.plot( dataX[:, -1], dataY[:, -1],  'k', linewidth = 2  )
#plt.plot( dataX[:,  0], dataY[:,  0],  'k', linewidth = 2  )
#plt.plot( dataX[-1, :], dataY[-1, :],  'k', linewidth = 2  )
#plt.plot( dataX[ 0, :], dataY[ 0, :],  'k', linewidth = 2  )






plt.tick_params(axis='both',which='major')

#plt.ylabel( "Latitude $(^{\circ})$" )
#plt.xlabel( "Longitude $(^{\circ})$" )
plt.ylabel( "Horizontal Distance-Y (km)" )
plt.xlabel( "Horizontal Distance-X (km)" )


if name == 0:
	fileName = 'CGFDM3D' 
if name == 1:
	fileName = 'LSRK' 
if name == 2:
	fileName = 'CJMVS' 
if name == 3:
	fileName = 'LSRK-CJMVS' 


if var == 'Vx':
	varLow = "$v_x$"
if var == 'Vy':
	varLow = "$v_y$"
if var == 'Vz':
	varLow = "$v_z$"




stationFile = open( "station.json" )
stationDir = json.load( stationFile )

station = stationDir["station(point)"] 

stationKeys     = station.keys( )
stationIndexes  = station.values( )


recvX = []
recvY = []
recvZ = []








for key in Keys:
	value = station[key]
	x = ( value[0] - centerX ) * DH
	y = ( value[1] - centerY ) * DH
	recvX.append(x)
	recvY.append(y)

lenRecv = len( recvX )

print( recvX )
print( recvY )



#plt.gca( ).set_ylim( [-125, 125] )
#plt.gca( ).set_xlim( [-125, 125] )

#cb.ax.tick_params( width = 0.25, length = 0.5  )


if DrawDiff:
	cb.ax.set_title( "Diff %s (m/s)" % varLow, loc = 'left' )
	cbRange = np.linspace( 0, vm, 5 )
	cb.set_ticks( cbRange )
	plt.title( fileName + " Diff: t = %.2fs" % ( it * DT ), fontsize = 14 )
	plt.savefig("PDF/Menyuan_%d_" % it +  fileName + var + "Diff" + ".pdf", dpi = dpi )
else:
	#plt.plot( srcX / unit, srcZ / unit, 'k*', markersize=10 )
	#for iRecv in range( lenRecv ):
	#	plt.text( recvX[iRecv] / unit, recvY[iRecv] / unit, StationLabels[iRecv] )
	#	plt.plot( recvX[iRecv] / unit, recvY[iRecv] / unit, colors[name] + 'v', markersize=5 )

	cb.ax.set_title( "%s (m/s)" % varLow , loc = 'center' )
	cbRange = np.linspace( -vm, vm, 5 )
	cb.set_ticks( cbRange )
	plt.title( fileName + ": t = %.2fs" % ( it * DT ), fontsize = 14 )
	plt.savefig( "PDF/Menyuan_%d_" % it +  fileName + var +  ".pdf", dpi = dpi )
