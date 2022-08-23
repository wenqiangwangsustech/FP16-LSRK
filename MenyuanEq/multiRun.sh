#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Thu 28 Jul 2022 02:11:52 PM CST
################################


python modifyParam.py


MakeFileDir='./MakefileDir'

n=0
#declare -A strArray

for cgfdm3dList in `ls ${MakeFileDir}'/'Makefile_* `
do
	if [ -f $cgfdm3dList ]; then
		echo ${cgfdm3dList/${MakeFileDir}'/'Makefile_/}
		strArray[n]=${cgfdm3dList/${MakeFileDir}'/'Makefile_/}
		#echo ${strArray[n]}
		n=$(($n+1))
		#echo $n
	fi
done

runDir="RunDir"
runShell="run"

: '
for((i=0;i<$n;i++))
do
	bash $runDir'/'$runShell${strArray[i]}
done
'
#bash RunDir/runCGFDM3D
bash RunDir/runCGFDM3D-LSRK-CJMVS

