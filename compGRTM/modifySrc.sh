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

file1=init_gpu.cpp
file2=readParams.cpp
file3=run.cpp
#echo ${strArray[4]}


SRC=./src
str="params"
suffix='.json'

prefix=".\/paramsDir"

#for cgfdm3dList in ${strArray[*]}
n=0
for srcList in `ls $SRC`
do
	#grep $str  $SRC'/'$srcList'/'$file1
	#grep $str  $SRC'/'$srcList'/'$file2
	#grep $str  $SRC'/'$srcList'/'$file3
	#grep $str${strArray[n]} $SRC'/'$srcList'/'$file1 
	#grep $str${strArray[n]} $SRC'/'$srcList'/'$file2
	#grep $str${strArray[n]} $SRC'/'$srcList'/'$file3

	sed -i s/$str${strArray[n]}/$prefix"\/"$str${strArray[n]}/g $SRC'/'$srcList'/'$file1 
	sed -i s/$str${strArray[n]}/$prefix"\/"$str${strArray[n]}/g $SRC'/'$srcList'/'$file2 
	sed -i s/$str${strArray[n]}/$prefix"\/"$str${strArray[n]}/g $SRC'/'$srcList'/'$file3 
	#sed s/$str${strArray[n]}/$prefix'/'$str${strArray[n]}/g $SRC'/'$srcList'/'$file1 > file1 
	#sed s/$str${strArray[n]}/$prefix'/'$str${strArray[n]}/g $SRC'/'$srcList'/'$file2 > file2
	#sed s/$str${strArray[n]}/$prefix'/'$str${strArray[n]}/g $SRC'/'$srcList'/'$file3 > file3
	#sed s/$str${strArray[n]}/$prefix$str${strArray[n]}/g >file1 #$SRC'/'$srcList'/'$file1 
	#sed s/$str${strArray[n]}/$prefix$str${strArray[n]}/g >file2 #$SRC'/'$srcList'/'$file2
	#sed s/$str${strArray[n]}/$prefix$str${strArray[n]}/g >file3 #$SRC'/'$srcList'/'$file3
	#sed -i s/$str$suffix/$str${strArray[n]}$suffix/g $SRC'/'$srcList'/'$file1 
	#sed -i s/$str$suffix/$str${strArray[n]}$suffix/g $SRC'/'$srcList'/'$file2
	#sed -i s/$str$suffix/$str${strArray[n]}$suffix/g $SRC'/'$srcList'/'$file3
	#echo s/${str}${suffix}/${str}$suffix${strArray[n]}/g #$SRC'/'$srcList'/'$file1 >> file1 
	#echo s/${str}${suffix}/${str}$suffix${strArray[n]}/g #$SRC'/'$srcList'/'$file2 >> file2
	#echo s/${str}${suffix}/${str}$suffix${strArray[n]}/g #$SRC'/'$srcList'/'$file3 >> file3
	n=$(($n+1))
	#echo $cgfdm3dList
done
