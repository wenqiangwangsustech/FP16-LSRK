#!/bin/env python

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Thu 14 Apr 2022 03:19:54 PM CST
#-------------------------------------------------------------------------------

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as pl

llr = np.loadtxt('LonLatRange.txt')
lon = np.load('lon.npy')
lat = np.load('lat.npy')
lomin = np.min( lon )
lomax = np.max( lon )
lamin = np.min( lat )
lamax = np.max( lat )

pgv = np.load('logPGV.npy')

nla, nlo = pgv.shape


print( np.max( pgv ) )


nPml = 12
ter1 = np.load('terrain.npy')
ter = ter1[nPml:nla+nPml, nPml:nlo+nPml]
print( np.shape( ter ) )


lo = np.linspace(lomin, lomax, nlo)
la = np.linspace(lamin, lamax, nla)

def saveasnc(fname, who, what):
  ncid = nc.Dataset(fname, 'w', format = 'NETCDF4')
  ncid.createDimension('lon', nlo)
  ncid.createDimension('lat', nla)
  loid = ncid.createVariable('lon', 'f4', ('lon',))
  laid = ncid.createVariable('lat', 'f4', ('lat',))
  vid = ncid.createVariable(who, 'f4', ('lat', 'lon',), fill_value = np.nan)
  loid.actual_range = [lomin, lomax]
  laid.actual_range = [lamin, lamax]
  vid.actual_range = [np.amin(what), np.amax(what)]
  ncid.Conventions = 'COARDS/CF-1.0'
  loid[:] = lo
  laid[:] = la
  vid[:] = what
  print( np.shape( what ) )
  ncid.close()

print('-R%.4f/%.4f/%.4f/%.4f' % (lomin, lomax, lamin, lamax))
saveasnc('PGV.nc', 'PGV', pgv)
saveasnc('ter.nc', 'terrain', ter)

#-------------------------------------------------------------------------------
fout = open('ann.txt', 'wt')
for i in range(-6, 1, 2):
  fout.write( '%3.1f a %5.3f\n' % (i, np.exp(i)) )
fout.close()

#-------------------------------------------------------------------------------
pgv = [ 0.02, 0.04, 0.09, 0.18, 0.35, 0.71, 1.41 ]

fout = open('cnt.lev', 'wt')
for pgvi in pgv:
  fout.write( '%f C\n' % ( np.log(pgvi) ) )
fout.close()

# vim:ft=python tw=80 ts=4 sw=2 et ai
