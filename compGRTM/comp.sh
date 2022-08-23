#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Wed 27 Jul 2022 04:32:24 PM CST
################################


MakeFileDir='./MakefileDir'
OBJDir=./obj

for mlist in `ls ${MakeFileDir}'/'Makefile_*`
do 
	if [ -f $mlist ]; then
		if [ $1 == "clean" ]; then
			echo "clean obj......."
			make -f $mlist  clean
			rm CompareData/* -rf
			mkdir CompareData/GRTM3D
		fi
		if [ $1 == "make"  ]; then
			echo "compiling......."
			mkdir ${OBJDir}'/'${mlist/${MakeFileDir}'/'Makefile_/}
			make -f $mlist -j
		fi
		#######mkdir ./obj/${mlist/Makefile_/}
	fi
done

