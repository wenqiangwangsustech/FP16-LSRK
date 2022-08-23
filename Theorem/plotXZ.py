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


it = 600

Eq = "Original Equation"
#Eq = 'New Equation'
#Eq = 'Transformed Variables'

var = 'Txx' #Vx Vz
name = 0
Cv = 1e2
Cs = 1e-5

if var[0] == 'V':
	C = Cv
else:
	C = Cs

eq = 0

if len( sys.argv ) > 1:
	it = int( sys.argv[1] )
	var = str( sys.argv[2] )
	eq = int( sys.argv[3] )
    
    
if eq == 0:
	Eq = "Original Equation"
if eq == 1:
	Eq = 'New Equation'
if eq == 2:
	Eq = 'Transformed Variables'
	


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 1


fileNameX = params["out"] + "/coordX"
fileNameY = params["out"] + "/coordY"
fileNameZ = params["out"] + "/coordZ"

DT = params["DT"]

varname = "./png_theorem/%s_%d" % ( var, it )

fileName = params["out"] + "/%s_%d"%( var, it )
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
else:
	dataX = np.zeros( [grid.NZ, grid.NX] )
	dataZ = np.zeros( [grid.NZ, grid.NX] )
	data  = np.zeros( [grid.NZ, grid.NX] )


mpiY = mpiSliceY
for mpiZ in range( grid.PZ ):
	for mpiX in range( grid.PX ):
		fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		print( fileName )
		nx = grid.nx[mpiX]
		nz = grid.nz[mpiZ]
		print( "nx = %d, nz = %d" % ( nx, nz ) )
		datax = np.fromfile( fileX, dtype='float32', count = nx * nz )
		#print( np.shape( datax ) )
		dataz = np.fromfile( fileZ, dtype='float32', count = nx * nz )
		data_  = np.fromfile( file, dtype='float32', count = nx * nz )
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx
		K  = grid.frontNZ[mpiZ]
		K_ = grid.frontNZ[mpiZ] + nz

		if FAST_AXIS == 'Z':
			dataX[I:I_, K:K_] = np.reshape( datax, ( nx, nz ) )
			dataZ[I:I_, K:K_] = np.reshape( dataz, ( nx, nz ) )
			data [I:I_, K:K_] = np.reshape( data_, ( nx, nz ) )
		else:
			dataX[K:K_, I:I_] = np.reshape( datax, ( nz, nx ) )
			dataZ[K:K_, I:I_] = np.reshape( dataz, ( nz, nx ) )
			data [K:K_, I:I_] = np.reshape( data_, ( nz, nx ) )


#data = data / C
dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
vm = np.max( np.abs( data ) )
dataX = dataX/unit
dataZ = dataZ/unit
#plt.pcolormesh( dataX, dataZ, data, cmap = "jet" )

'''
plt.plot( dataX[:, -1], dataZ[:, -1],  'k', linewidth = 2  )
plt.plot( dataX[:,  0], dataZ[:,  0],  'k', linewidth = 2  )
plt.plot( dataX[-1, :], dataZ[-1, :],  'k', linewidth = 2  )
plt.plot( dataX[ 0, :], dataZ[ 0, :],  'k', linewidth = 2  )
'''



plt.pcolormesh( dataX, dataZ, data, vmax = vm, vmin = -vm, cmap = "seismic" )
#plt.pcolormesh( data[::sample, ::sample] // unit, vmin = -vm / 2, vmax = vm, cmap = "jet" )
#plt.pcolormesh( data[5:grid.NZ -5:sample, 5:grid.NX - 5:sample] // unit, cmap = "seismic" )
#plt.pcolormesh( data, vmax = vm, vmin = -vm, cmap = "jet" )
#plt.imshow( data, cmap = "jet", origin= "lower" )

cb = plt.colorbar( format='%.0e')
cbRange = np.linspace( -vm, vm, 5 )
cb.set_ticks( cbRange )
#cb.ax.tick_params(labelsize=14)
#plt.colorbar( orientation = "horizontal" )
plt.axis( "image" )
#fileName = "%s: t = %.3fs" % (var, it * DT )
#fileName = "%s: %s" % ( Eq, var )

if Eq == 'Original Equation':
	if var == 'Txz':
		fileName = "%s: $\sigma_{xz}$" % Eq
	if var == 'Vx':
		fileName = "%s: $v_x$" % Eq
	if var == 'Vz':
		fileName = "%s: $v_z$" % Eq


if Eq == 'New Equation':
	if var == 'Txz':
		fileName = "%s: $\Sigma_{xz}$" % Eq
	if var == 'Vx':
		fileName = "%s: $V_x$" % Eq
	if var == 'Vz':
		fileName = "%s: $V_z$" % Eq



if Eq == 'Transformed Variables':
	if var == 'Txz':
		fileName = "%s: $\Sigma_{xz}/C_v$" % Eq
	if var == 'Vx':
		fileName = "%s: $V_x/C_v$" % Eq
	if var == 'Vz':
		fileName = "%s: $V_z/C_v$" % Eq


Xmin = np.min( dataX )
Xmax = np.max( dataX )
Zmin = np.min( dataZ )
Zmax = np.max( dataZ )
xticksRange = np.linspace( Xmin, Xmax, 5)
 
zticksRange = np.linspace( Zmin, Zmax, 5 )

plt.xticks(xticksRange)
plt.yticks(zticksRange)

plt.tick_params(axis='both',which='major')

plt.ylabel( "Vertical Distance-Z (km)" )
plt.xlabel( "Horizontal Distance-X (km)" )


'''
if name == 0:
	fileName = 'CGFDM3D' 
if name == 1:
	fileName = 'CGFDM3D-LSRK' 
if name == 2:
	fileName = 'CGFDM3D-CJMVS' 
if name == 3:
	fileName = 'CGFDM3D-LSRK-CJMVS' 
'''

plt.title( Eq + ':' + var, fontsize = 14 )

plt.savefig( varname + Eq + ".pdf" )
#plt.savefig( fileName + var + ".pdf" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


