#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Fri 29 Jul 2022 02:56:14 PM CST
################################


n=4
it=600

python compPlotXZ.py ${it} Vx   0 0
python compPlotXZ.py ${it} Vz   0 0
python compPlotXZ.py ${it} Txz  0 0


python compPlotXZ.py ${it} Vx   0 1
python compPlotXZ.py ${it} Vz   0 1
python compPlotXZ.py ${it} Txz  0 1


python compPlotXZ.py ${it} Vx   0 2
python compPlotXZ.py ${it} Vz   0 2
python compPlotXZ.py ${it} Txz  0 2


python compPlotXZ.py ${it} Vx   1 2
python compPlotXZ.py ${it} Vz   1 2
python compPlotXZ.py ${it} Txz  1 2


