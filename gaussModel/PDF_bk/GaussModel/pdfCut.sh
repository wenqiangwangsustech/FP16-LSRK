#!/usr/bin/bash
################################
## Author: Wenqiang Wang
## Mail: 11849528@mail.sustech.edu.cn
## Created Time: Thu 11 Aug 2022 10:45:44 PM CST
################################


for pdfFile in `ls *.pdf`
do
	pdfcrop $pdfFile 
done

