#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Fri 29 Jul 2022 02:56:14 PM CST
################################


n=4
var=Vx
it=2500

: '
for((i=0;i<$n;i++))
do
	python plotMutiSeismos.py Vx
	python plotMutiSeismos.py Vy
	python plotMutiSeismos.py Vz
	#python compPlotXY.py ${it} ${var} $i 1
	#python compPlotXY.py ${it} ${var} $i 0
done
'

python plotMutiSeismos.py Vx
python plotMutiSeismos.py Vy
python plotMutiSeismos.py Vz
#python compPlotXY.py ${it} ${var} $i 1

NT=4001
interval=2000

for((it=2000;it<$NT;it=it+$interval))
do
	echo $it
	python compPlotXY.py ${it} ${var} 0 0
	python compPlotXY.py ${it} ${var} 3 0
	python compPlotXY.py ${it} ${var} 0 1
	python compPlotXY.py ${it} ${var} 3 1
done

