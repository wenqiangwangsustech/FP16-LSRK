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


it = 1000

Eq = "Original Equation"
#Eq = 'New Equation'
#Eq = 'Transformed Variables'


StationLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 
				 'H', 'I', 'J', 'K', 'L', 'M', 'N',
				 'O', 'P', 'Q',		 'R', 'S', 'T',
				 'U', 'V', 'W', 	 'X', 'Y', 'Z']

Keys = [ '7', '9', '11', '13', '15' ]# '13', '14', '15', '16', '17' ]

colors = [ 'k', 'r', 'g', 'b', 'y', 'm']

var = 'Vx' #Vx Vz
name = 1

DrawDiff = 0

if len( sys.argv ) > 1:
	it = int( sys.argv[1] )
	var = str( sys.argv[2] )
	name = int( sys.argv[3] )
	DrawDiff = int( sys.argv[4] )



jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 1

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
fileNameZ = params["out"] + "/coordZ"

DT = params["DT"]

varname = "./png_theorem/%s_%d" % ( var, it )

fileNameTh = params["out"] + "/%s_%d"%( var, it )
fileName = outputPath + "/%s_%d"%( var, it )
#fileName = params["out"] + "/Vs"
#fileName = "output_homo_original" +"/%s_%d"%( var, it )
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
	dataX = np.zeros( [grid.NX, grid.NZ] )
	dataZ = np.zeros( [grid.NX, grid.NZ] )
	data  = np.zeros( [grid.NX, grid.NZ] )
	dataTh= np.zeros( [grid.NX, grid.NZ] )
else:
	dataX = np.zeros( [grid.NZ, grid.NX] )
	dataZ = np.zeros( [grid.NZ, grid.NX] )
	data  = np.zeros( [grid.NZ, grid.NX] )
	dataTh= np.zeros( [grid.NZ, grid.NX] )


mpiY = mpiSliceY
for mpiZ in range( grid.PZ ):
	for mpiX in range( grid.PX ):
		fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		fileTh= open( "%s_Y_mpi_%d_%d_%d.bin" % (fileNameTh, mpiX, mpiY, mpiZ ), "rb" )
		print( fileName )
		print( fileNameTh )
		nx = grid.nx[mpiX]
		nz = grid.nz[mpiZ]
		print( "nx = %d, nz = %d" % ( nx, nz ) )
		datax = np.fromfile( fileX, dtype='float32', count = nx * nz )
		#print( np.shape( datax ) )
		dataz = np.fromfile( fileZ, dtype='float32', count = nx * nz )
		data_  = np.fromfile( file, dtype='float32', count = nx * nz )
		data_Th= np.fromfile(fileTh, dtype='float32', count = nx * nz )
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx
		K  = grid.frontNZ[mpiZ]
		K_ = grid.frontNZ[mpiZ] + nz

		if FAST_AXIS == 'Z':
			dataX[I:I_, K:K_] = np.reshape( datax, ( nx, nz ) )
			dataZ[I:I_, K:K_] = np.reshape( dataz, ( nx, nz ) )
			data [I:I_, K:K_] = np.reshape( data_, ( nx, nz ) )
			dataTh[I:I_, K:K_] = np.reshape( data_Th, ( nx, nz ) )
		else:
			dataX[K:K_, I:I_] = np.reshape( datax, ( nz, nx ) )
			dataZ[K:K_, I:I_] = np.reshape( dataz, ( nz, nx ) )
			data [K:K_, I:I_] = np.reshape( data_, ( nz, nx ) )
			dataTh[K:K_, I:I_] = np.reshape( data_Th, ( nz, nx ) )


#data = data / C
dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
if DrawDiff:
	vm = np.max( np.abs( data - dataTh ) ) #* 1e-2
else:
	vm = np.max( np.abs( dataTh ) ) #* 1e-2
	

dataX = dataX/unit
dataZ = dataZ/unit
#plt.pcolormesh( dataX, dataZ, data, cmap = "jet" )

plt.plot( dataX[:, -1], dataZ[:, -1],  'k', linewidth = 2  )
plt.plot( dataX[:,  0], dataZ[:,  0],  'k', linewidth = 2  )
plt.plot( dataX[-1, :], dataZ[-1, :],  'k', linewidth = 2  )
plt.plot( dataX[ 0, :], dataZ[ 0, :],  'k', linewidth = 2  )





sourceX = params['sourceX']
sourceY = params['sourceY']
sourceZ = params['sourceZ']

centerX = params['centerX']
centerY = params['centerY']
DH = params['DH']

srcX = ( sourceX - centerX ) * DH 
srcY = ( sourceY - centerY ) * DH 
srcZ = - ( grid.NZ - 1 - sourceZ ) * DH 

if DrawDiff:
	plt.pcolormesh( dataX, dataZ, np.abs( data - dataTh ), vmax = vm, vmin = 0, cmap = "Reds" )
else:
	plt.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" )
#plt.pcolormesh( data[::sample, ::sample] // unit, vmin = -vm / 2, vmax = vm, cmap = "jet" )
#plt.pcolormesh( data[5:grid.NZ -5:sample, 5:grid.NX - 5:sample] // unit, cmap = "seismic" )
#plt.pcolormesh( data, vmax = vm, vmin = -vm, cmap = "jet" )
#plt.imshow( data, cmap = "jet", origin= "lower" )

#cb = plt.colorbar( orientation =  "horizontal", format = '%.0e',  fraction = 0.05, pad = 0.15 )
cb = plt.colorbar( format = '%.1e', shrink = 0.8 )
#cb.ax.tick_params(labelsize=14)
#plt.colorbar( orientation = "horizontal" )
plt.axis( "image" )
#fileName = "%s: t = %.3fs" % (var, it * DT )
#fileName = "%s: %s" % ( Eq, var )


Xmin = np.min( dataX )
Xmax = np.max( dataX )
Zmin = np.min( dataZ )
Zmax = np.max( dataZ )
xticksRange = np.linspace( -10.0, 10.0, 5)
 
zticksRange = [ -10, -8, -6, -4, -2, 0, 3 ] #np.linspace( -10, 2.5, 5 )

plt.xticks(xticksRange)
plt.yticks(zticksRange)


plt.tick_params(axis='both',which='major')

plt.ylabel( "Vertical Distance-Z (km)" )
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

cal_depth = params["Depth(km)"] * unit
h = 0.2 * cal_depth
a = 0.1 * cal_depth
b = 0.1 * cal_depth


for key in Keys:
	value = station[key]
	recvX.append( ( value[0] - centerX ) * DH  )
	recvY.append( ( value[1] - centerY ) * DH  )
	x = ( value[0] - centerX ) * DH
	y = ( value[1] - centerY ) * DH
	print( x )
	z = h * np.exp( - 0.5 * ( x * x / ( a * a ) + y * y / ( b * b ) ) )
	recvZ.append( z ) 

lenRecv = len( recvX )

print( recvX )
print( recvZ )



plt.gca( ).set_ylim( [-10.5, 3.1] )
plt.gca( ).set_xlim( [-10.5, 10.5] )

if DrawDiff:
	cb.ax.set_title( "Diff %s (m/s)" % varLow )
	cbRange = np.linspace( 0, vm, 5 )
	cb.set_ticks( cbRange )
	plt.title( fileName + " Diff: t = %.2fs" % ( it * DT ), fontsize = 14 )
	plt.savefig("PDF/gauss" +  fileName + var + "Diff" + ".pdf", bbox_inches = 'tight' )
else:
	cb.ax.set_title( "%s (m/s)" % varLow )
	plt.plot( srcX / unit, srcZ / unit, 'k*', markersize=10 )
	for iRecv in range( lenRecv ):
		print( Keys[iRecv] + ":" + "x = %f, z = %f" % (recvX[iRecv] / unit, recvZ[iRecv] / unit )  )
		plt.text( recvX[iRecv] / unit - 0.3, recvZ[iRecv] / unit + 0.3, StationLabels[iRecv] )
		plt.plot( recvX[iRecv] / unit, recvZ[iRecv] / unit, colors[name] + 'v', markersize=5 )

	cbRange = np.linspace( -vm, vm, 5 )
	cb.set_ticks( cbRange )
	plt.title( fileName + ": t = %.2fs" % ( it * DT ), fontsize = 14 )
	plt.savefig( "PDF/gauss" + fileName + var +  ".pdf", bbox_inches = 'tight' )



#plt.savefig( varname + Eq + ".pdf" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


