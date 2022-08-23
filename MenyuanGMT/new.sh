#!/bin/bash

J1=m5c
J2=M2c
R1=99.94/102.71/36.66/38.84
R2=95/125/28/45
B="-Ba0.5f0.25 -BnWSe"
grd=/data0/home/wangwq/earth_relief_15s.grd

gmt begin map pdf

  gmt set FONT_ANNOT_PRIMARY=10p,4
  gmt set MAP_FRAME_WIDTH=0.1c

  # 底图
  gmt grdgradient $grd -A0/90 -R$R1 -Nt -Ginten.grd
  #gmt grdimage $grd -R$R1 -J$J1 $B -Cdem2 -Iinten.grd # elevation, dem2, grayC

  gmt makecpt -Cdem1 -T1300/5000/100 -Z 
  gmt grdimage $grd -R$R1 -J$J1 $B  -Iinten.grd # elevation, dem2, grayC

  # 比例尺
  gmt basemap -R$R1 -J$J1 -Lg102.30/38.75+c38.75+w40k+f+u

  # 色标
  gmt colorbar -DjCB+w10c/0.2c+o0/-1.5c+h+e+m -Bxaf

  # 断层
  gmt psxy -R$R1 -J$J1 -W0.5p,black -Gwhite  << EOF
	101.526040  37.717425
	100.966363  37.870343
	100.974772  37.889632
	101.534447	37.736710
	101.526040  37.717425
EOF


: '
  gmt text -R$R1 -J$J1 -F+f14p+j -Dj0.15c/0.15c << EOF
	101.526040  37.717425   BC   A
	100.966363  37.870343   TC   B
EOF
'

	
  # 震源
  echo 101.3 37.8 | gmt  psxy -Sa0.2 -W2.5p,green -Ggreen



  # 台站
  cat stationLonLat  | gmt  psxy -Si0.2 -W2.5p,red -Gred

  cat textStationLonLat | gmt text -R$R1 -J$J1 -F+f14p+j -Dj0.1c/0.1c 

  # 子图
  gmt coast -R$R2 -J$J2 -Ggray -Swhite -Na -Bnwse --MAP_FRAME_TYPE=inside
  gmt psxy -R$R2 -J$J2 -Wblack << EOF
     99.94   36.66
    102.71   36.66
    102.71   38.84
     99.94   38.84
     99.94   36.66
EOF



gmt end show

rm -f inten.grd
#lonStart = 101.526040, lonEnd = 100.966363, latStart = 37.717425, latEnd = 37.870343
