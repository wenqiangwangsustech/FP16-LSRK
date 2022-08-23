#!/usr/bin/env python
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Wed 27 Jul 2022 08:51:10 PM CST
################################


import numpy as np
import matplotlib.pyplot as plt
import os
import json

exeFile = [
'CGFDM3D',
'CGFDM3D-CJM',
'CGFDM3D-CJMVS', 
'CGFDM3D-LSRK', 
'CGFDM3D-LSRK-CJM', 
'CGFDM3D-LSRK-CJMVS']

jsonsFile = open( "params.json" )
params = json.load( jsonsFile )



for iExe in exeFile:
	params['out'] = 'CompareData/%s' % iExe
	
	print( params )
	paramsJsonFile = open( "paramsDir/params%s.json" % iExe, 'w' )
	
	json_str = json.dumps( params, indent = 4 )
	paramsJsonFile.write( json_str )
	paramsJsonFile.close( )
