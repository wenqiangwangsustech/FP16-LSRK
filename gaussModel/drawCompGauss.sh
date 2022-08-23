#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Fri 29 Jul 2022 02:56:14 PM CST
################################


n=4
var=Vx
it=1000

for((i=0;i<$n;i++))
do
	python plotMutiSeismos.py Vx
	python plotMutiSeismos.py Vz
	python compPlotXZ.py ${it} ${var} $i 1
	python compPlotXZ.py ${it} ${var} $i 0
done


