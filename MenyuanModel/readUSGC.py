#!/usr/bin/python
'''
Wenqiang Wang@SUSTech 

Date: 2021/12/21

'''


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import json as js
import struct
from pyproj import Proj
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator

def gaussFunc( t, t0, a ):
    source = np.exp(-( t - t0 )**2 / ( a ** 2 ) ) / ( np.sqrt( np.pi ) * a )
    return source



#right-hand rule
def rotateX( xTheta, x, y, z ):
	newX = x
	newY = x * 0 + y * np.cos( xTheta ) - z * np.sin( xTheta )
	newZ = x * 0 + y * np.sin( xTheta ) + z * np.cos( xTheta ) 
	return newX, newY, newZ

def rotateY( yTheta, x, y, z ):
	newX = x * np.cos( yTheta ) + z * np.sin( yTheta )
	newY = y
	newZ =-x * np.sin( yTheta ) + z * np.cos( yTheta )
	return newX, newY, newZ

def rotateZ( zTheta, x, y, z ):
	newX = x * np.cos( zTheta ) - y * np.sin( zTheta )
	newY = x * np.sin( zTheta ) + y * np.cos( zTheta )
	newZ = z
	return newX, newY, newZ
	

def localElementMesh( Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, n1, n2 ):
	xmesh = np.zeros( [n2, n1] ) 
	ymesh = np.zeros( [n2, n1] ) 
	zmesh = np.zeros( [n2, n1] ) 

	xmesh[:, 0] = np.linspace( Ax, Cx, n2 )
	ymesh[:, 0] = np.linspace( Ay, Cy, n2 )
	zmesh[:, 0] = np.linspace( Az, Cz, n2 )

	xmesh[:, -1] = np.linspace( Bx, Dx, n2 )
	ymesh[:, -1] = np.linspace( By, Dy, n2 )
	zmesh[:, -1] = np.linspace( Bz, Dz, n2 )

	for j in range( n2 ):
		xmesh[j, :] = np.linspace( xmesh[j,0], xmesh[j,-1], n1 ) 
		ymesh[j, :] = np.linspace( ymesh[j,0], ymesh[j,-1], n1 ) 
		zmesh[j, :] = np.linspace( zmesh[j,0], zmesh[j,-1], n1 )

	return xmesh, ymesh, zmesh 
	#return xmesh[:-1, :-1], ymesh[:-1, :-1], zmesh[:-1, :-1] 

def thickenMesh( nx, ny, thickenX, thickenY, X, Y, Z, XMesh, YMesh, ZMesh ):
	NX = ( nx - 1 ) * thickenX + 1
	NY = ( ny - 1 ) * thickenY + 1
	for j in range( 1, ny ):
		for i in range( 1, nx ):
			J_ = ( j - 1 ) * thickenY
			J  = j * thickenY
			I_ = ( i - 1 ) * thickenX
			I  = i * thickenX
			index = i - 1 + ( j - 1 ) * nx
			Ax, Ay, Az = X[index], Y[index], Z[index] 
			index = i + ( j - 1 ) * nx
			Bx, By, Bz = X[index], Y[index], Z[index] 
			index = i - 1 + j * nx
			Cx, Cy, Cz = X[index], Y[index], Z[index]
			index = i + j * nx
			Dx, Dy, Dz = X[index], Y[index], Z[index]
			n1 = thickenX
			n2 = thickenY

			XCoord, YCoord, ZCoord = localElementMesh( Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, n1 + 1, n2 + 1 )
			for jj in range( J_, J + 1 ):
				for ii in range( I_, I + 1 ):
					index = ii + jj * NX
					XMesh[index], YMesh[index], ZMesh[index] = XCoord[jj - J_, ii - I_], YCoord[jj - J_, ii - I_], ZCoord[jj - J_, ii - I_], 



NT = 1000
strike = 289
dip = 84


lengthPoint = 27
depthPoint = 12


thickenX = 10
thickenZ = 10


file = open( "FFM.geojson", "r" )
finiteFaultTest = file.read(  )

#print( finiteFaultTest  )

finiteFaultDir = js.loads( finiteFaultTest )

print( len( finiteFaultDir ) )

for key in finiteFaultDir.keys( ):
	print( key )

faultInfo = finiteFaultDir['features']

src_npts = len( faultInfo )

lonRect   = np.zeros( [ src_npts, 4] )
latRect   = np.zeros( [ src_npts, 4] )
depthRect = np.zeros( [ src_npts, 4] )

rakeCenter 	  = np.zeros( src_npts )

latRectCenter	= np.zeros( src_npts )
lonRectCenter	= np.zeros( src_npts )
depthRectCenter = np.zeros( src_npts )

timeFieldCenter	= np.zeros( src_npts )
halfDurationCenter = np.zeros( src_npts )
slipCenter = np.zeros( src_npts )

#print( faultInfo[0]["geometry"]["coordinates"][0][1] )


for j in range( src_npts ):
	for i in range( 4 ):
		lonRect  [j, i] = faultInfo[j]["geometry"]["coordinates"][0][i][0]
		latRect  [j, i] = faultInfo[j]["geometry"]["coordinates"][0][i][1]
		depthRect[j, i] = faultInfo[j]["geometry"]["coordinates"][0][i][2]
	timeFieldCenter[j] = faultInfo[j]["properties"]["trup"]
	rakeCenter[j] = faultInfo[j]["properties"]["rake"]
	halfDurationCenter[j] = faultInfo[j]["properties"]["rise"] * 0.5
	slipCenter[j] = faultInfo[j]["properties"]["slip"]

depthRect = -depthRect

lonRectCenter = ( lonRect[:, 0] + lonRect[:, 1] + lonRect[:, 2] + lonRect[:, 3] ) * 0.25
latRectCenter = ( latRect[:, 0] + latRect[:, 1] + latRect[:, 2] + latRect[:, 3] ) * 0.25
depthRectCenter = ( depthRect[:, 0] + depthRect[:, 1] + depthRect[:, 2] + depthRect[:, 3] ) * 0.25 

fig = plt.figure( 1 )
ax = Axes3D( fig )
ax.scatter( lonRectCenter, latRectCenter, depthRectCenter ) 

if ( src_npts != lengthPoint * depthPoint ):
	print( "src_npts = %d, lengthPoint = %d, depthPoint = %d" %( src_npts, lengthPoint, depthPoint ) )
	print( "The lengthPoint and depthPoint you set is wrong"  )
	sys.exit( )


nx = lengthPoint + 1
nz = depthPoint + 1
npts = nz * nx

lon  = np.zeros( npts )
lat  = np.zeros( npts )
depth = np.zeros( npts )



slip = np.zeros( npts )
timeField = np.zeros( npts )
halfDuration = np.zeros( npts )
Area = np.zeros( npts )
rake = np.zeros( npts  )

for j in range( nz - 1 ):
	for i in range( nx - 1 ):
		index = i + j * nx
		pos = i + j * lengthPoint
		lon[index]   =   lonRect[pos, 0]
		lat[index]   =   latRect[pos, 0]
		depth[index] = depthRect[pos, 0]
		

for j in range( 1, nz ):
	for i in range( 1, nx ):
		index = i + j * nx
		pos = i - 1 + ( j - 1 ) * lengthPoint
		lon[index]   =   lonRect[pos, 2]
		lat[index]   =   latRect[pos, 2]
		depth[index] = depthRect[pos, 2]

i = nx - 1
j = 0
index = i + j * nx
pos = i - 1 + j * lengthPoint
lon[index]   =   lonRect[pos, 1]
lat[index]   =   latRect[pos, 1]
depth[index] = depthRect[pos, 1]






i = 0
j = nz - 1
index = i + j * nx
pos = i + ( j - 1 ) * lengthPoint
lon[index]   =   lonRect[pos, 3]
lat[index]   =   latRect[pos, 3]
depth[index] = depthRect[pos, 3]

lon0 = ( np.max( lon ) + np.min( lon ) )* 0.5 
lat0 = ( np.max( lat ) + np.min( lat ) )* 0.5 
print( "Longitude and latitude center: lon0 = %f, lat = %f" % ( lon0, lat0 ) )



proj = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84" )
X, Y = proj( lon, lat )
Z = depth


XCenter, YCenter = proj( lonRectCenter, latRectCenter )
ZCenter = depthRectCenter

xTheta = 0#- 45 * np.pi / 180
yTheta = -dip * np.pi / 180#45 * np.pi / 180
#yTheta = 0#45 * np.pi / 180
zTheta = strike * np.pi / 180

#right-hand rule
X, Y, Z = rotateZ( zTheta, X, Y, Z )
X, Y, Z = rotateX( xTheta, X, Y, Z )
X, Y, Z = rotateY( yTheta, X, Y, Z )

XCenter, YCenter, ZCenter = rotateZ( zTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateX( xTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateY( yTheta, XCenter, YCenter, ZCenter )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), slipCenter )
slip = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), halfDurationCenter )
halfDuration = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), timeFieldCenter )
timeField = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), rakeCenter )
rake = interp( X, Y )


M0 = 0.0
miu = 3e10
for j in range( npts ):
	M0 += slip[i] * 2e6 * miu
Mw = 2/3 * np.log10( M0 ) - 6.06
print( "The M0 is %e"% M0  )
print( "The Mw is %e"% Mw  )


fig = plt.figure( 2 )
ax = Axes3D( fig )
ax.scatter( X, Y, Z ) 


NX = ( nx - 1 ) * thickenX + 1
NY = ( nz - 1 ) * thickenZ + 1
NPTS = NX * NY
XMesh = np.zeros( NPTS )
YMesh = np.zeros( NPTS )
ZMesh = np.zeros( NPTS )

SlipMesh = np.zeros( NPTS )
HalfDurationMesh = np.zeros( NPTS )
TimeFieldMesh = np.zeros( NPTS )
AreaMesh = np.zeros( NPTS )
rakeMesh = np.zeros( NPTS )
SlipRateMesh = np.zeros( [NPTS, NT] )
RakeTimeMesh = np.zeros( [NPTS, NT] )
StrikeMesh = np.zeros( NPTS ) + strike
DipMesh = np.zeros( NPTS  ) + dip




XCoord = np.zeros( [NY, NX]  )
YCoord = np.zeros( [NY, NX]  )
ZCoord = np.zeros( [NY, NX]  )

LonCoord = np.zeros( [NY, NX]  )
LatCoord = np.zeros( [NY, NX]  )

data = np.zeros( [NY, NX] )
dataTime = np.zeros( [NY, NX] )
dataSlip = np.zeros( [NY, NX] )
dataArea = np.zeros( [NY, NX] )
dataRake = np.zeros( [NY, NX] )
dataSlipRate = np.zeros( [NY, NX, NT] )



thickenMesh( nx, nz, thickenX, thickenZ, X, Y, Z, XMesh, YMesh, ZMesh )




for j in range( 1, NY ):
	for i in range( 1, NX ):
		index2 = i + j * NX
		index  = i - 1 + ( j - 1 ) * NX
		x = np.abs( XMesh[index2] - XMesh[index] )
		y = np.abs( YMesh[index2] - YMesh[index] )
		AreaMesh[index] = x * y

i = NX - 1
for j in range( NY ):
	index  = i + j * NX
	index2 = i - 1 + j * NX
	AreaMesh[index] = AreaMesh[index2]


j = NY - 1
for i in range( NX ):
	index  = i + j * NX
	index2 = i + ( j - 1 ) * NX
	AreaMesh[index] = AreaMesh[index2]

interp = LinearNDInterpolator( list(zip(X, Y)), slip)
SlipMesh = interp( XMesh, YMesh )
interp = LinearNDInterpolator( list(zip(X, Y)), halfDuration)
HalfDurationMesh = interp( XMesh, YMesh )
interp = LinearNDInterpolator( list(zip(X, Y)), timeField)
TimeFieldMesh = interp( XMesh, YMesh )
interp = LinearNDInterpolator( list(zip(X, Y)), rake )
rakeMesh = interp( XMesh, YMesh )

print( "AreaMesh = %f" % AreaMesh[0]  )


M0 = 0.0
miu = 3e10
for j in range( NPTS ):
	M0 += SlipMesh[i] * AreaMesh[i] * miu

Mw = 2/3 * np.log10( M0  - 6.06)
print( "The Mw is %f", Mw  )

frontT = np.max( HalfDurationMesh ) * 4
backT = np.max( HalfDurationMesh ) * 4
dt = ( np.max( TimeFieldMesh ) + frontT + backT ) / NT


t = np.linspace( 0, NT * dt, NT )
tmp = np.zeros( NT ) + 1.0

for i in range( NPTS ):
	t0 = TimeFieldMesh[i] + frontT
	source = gaussFunc( t, t0, HalfDurationMesh[i] )#
	SlipRateMesh[i, :] = source * SlipMesh[i]
	RakeTimeMesh[i, :] = rakeMesh[i] * tmp


for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = XMesh[index]
		YCoord[j, i] = YMesh[index]
		ZCoord[j, i] = ZMesh[index]
		data[j, i] = SlipMesh[index]
		dataTime[j, i] = TimeFieldMesh[index]
		dataSlip[j, i] = SlipMesh[index]
		dataArea[j, i] = AreaMesh[index]
		dataRake[j, i] = rakeMesh[index]
		dataSlipRate[j, i, :] = SlipRateMesh[index, :]



XMesh, YMesh, ZMesh = rotateY( -yTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateX( -xTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateZ( -zTheta, XMesh, YMesh, ZMesh )

lonMesh = np.zeros( NPTS )
latMesh = np.zeros( NPTS )

proj1 = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84")
lonMesh, latMesh = proj1( XMesh, YMesh, inverse = "true" )



jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/USGSsource.bin"
sourceFile = open( sourceFileName, "wb" )

value = struct.pack( "i", NPTS )
sourceFile.write( value )

value = struct.pack( "i", NT )
sourceFile.write( value )

value = struct.pack( "f", dt )
sourceFile.write( value )

for i in range( NPTS ):
	value = struct.pack( "f",  lonMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  latMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  ZMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  AreaMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  StrikeMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  DipMesh[i] )
	sourceFile.write( value )
	tvalue = struct.pack( "f" * NT,  *( RakeTimeMesh[i, :] ) )
	sourceFile.write( tvalue )
	tvalue = struct.pack( "f" * NT,  *( SlipRateMesh[i, :] ) )
	sourceFile.write( tvalue )

sourceFile.close( )



for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		LonCoord[j, i] = lonMesh[index]
		LatCoord[j, i] = latMesh[index]
		

fig = plt.figure( 3 )
plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

plt.savefig( "fault.png" )



'''
plt.ion( )
plt.figure( )
for it in range( 0, NT, 10 ):
	plt.clf( )
	vm = np.max( dataSlipRate[:, :, it] )
	#vm = 0.1
	plt.pcolormesh( dataSlipRate[:, :, it], vmin = - vm, vmax = vm, cmap = "seismic")
	plt.xlabel( "x" )
	plt.ylabel( "y" )
	plt.colorbar( )
	plt.axis( "image" )
	plt.show( )
	plt.pause( 0.002 )
'''
#fig = plt.figure( 1 )
#plt.pcolormesh( Slip, cmap = "jet" )
#plt.colorbar( )
#plt.axis( "image" )









#fig = plt.figure( 0 )
#plt.gca().set_box_aspect( ( 250, 130, 4)  )




#fig = plt.figure( 1 )
#plt.scatter( XCoord, YCoord, XCoord )
#plt.pcolormesh( XCoord, YCoord, data, cmap = "jet" )
#plt.pcolormesh( XCoord, YCoord, abs( dataRake ),  cmap = "jet" )

#plt.contourf( XCoord, YCoord, abs( dataArea ), label = 15,  cmap = "jet" )

#plt.colorbar( )
#plt.pcolormesh( XCoord[:-1, :-1], YCoord[:-1, :-1], XCoord[:-1, :-1], cmap = "jet" )
#plt.axis( "equal" )

#plt.scatter( X, Y, slip )


