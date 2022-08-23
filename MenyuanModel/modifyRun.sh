#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Wed 27 Jul 2022 10:11:10 PM CST
################################

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
#sed 's/params.json/params'

echo '============================='

runShell="run"
runDir="RunDir"
log="log"
logDir="LogDir"

main="main"
binDir="bin"
if [ ! -d $logDir ]; then
	mkdir $logDir
fi

for((i=0;i<$n;i++))
do
	echo $i
	#grep $log $runShell -rn > fileRun
	#sed s/$log/aaa/g $run > fileRun
	sed s/$log/'.\/'$logDir'\/'$log${strArray[i]}/g $runShell | sed s/$main/${strArray[i]}/g  > $runDir/$runShell${strArray[i]}
done
