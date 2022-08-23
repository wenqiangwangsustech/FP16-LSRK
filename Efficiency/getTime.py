#!/usr/bin/env python
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Mon 08 Aug 2022 03:55:34 PM CST
################################


import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import sys

LogDir = "LogDir"

Logs = [ 'logCGFDM3D', 'logCGFDM3D-CJMVS',  'logCGFDM3D-LSRK',  'logCGFDM3D-LSRK-CJMVS' ]

Model = [ 'CGFDM3D', 'CJMVS',  'LSRK',  'LSRK-CJMVS' ]

Grid = ["25", "50", "100", "200", "400"]


colors = [ 'k', 'r', 'g', 'b']

nl = len( Logs )
ng = len( Grid )

T = np.zeros( [nl, ng]  )

j = 0
for grid in Grid:
	i = 0
	for log in Logs:
		fileStr = LogDir + '_' + grid + "/" + log
		#print( fileStr )
		logFile = open( fileStr )
		lines = logFile.readlines( )
		#print( len( lines ) )
		str = lines[ -2 ].split( )
		print( str[-1][:-2] )
		T[i, j] = float( str[-1][:-2] )
		i += 1

	j += 1	



print( T[0, :] )

logT = T#np.log10( T  )

x = np.array( Grid  )

xtic = []

for ii in range( len(Grid ) ):
	xtic.append( "$%s^3$" % Grid[ii] )


print( xtic )


'''
for i in range( nl ):
	plt.plot( x, logT[i, :], label = Model[i] , c = colors[i]  )
	plt.legend( loc = 1 )
'''
plt.figure( figsize = ( 10,5) )

plt.plot( x, logT[0, :], label = Model[0], c = colors[0]  )
plt.legend( loc = 1 )

plt.plot( x, logT[-1, :], label = Model[-1], c = colors[-1]  )
plt.legend( loc = 1 )

plt.xlabel( 'Grid Scales: $NX * NY * NZ$'  )
plt.ylabel( 'Time Cost (s)'  )

plt.xticks( x, xtic )
plt.title( "Performance" )
plt.savefig( "PDF/Performance.pdf", bbox_inches = "tight"    )
#plt.axis( "image"  )
plt.show( )


