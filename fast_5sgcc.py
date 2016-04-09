#!/usr/bin/env python
#file name: fast_5sgcc.py

from numpy import *

#fast_file="/home/niutj/python_pros/test_QC.fq"
fast_file="/1_disk/public_resources/jiankuihe-exome-1.fq"

rawdata=open(fast_file,"r")
rawdata.readline()
lnum=0
list1=zeros((100))

for line in rawdata:
    lnum+=1
    line=line.strip('\n')

    if(lnum % 4 ==1):
        line.upper()
        count=0
        for char in line:
            if char=='G' or char == 'C':
                count+=1
        list1[count]+=1
j=0
outfile=open('exome_1_5.txt','w')
outfile.write('percent\tnum\n')
for i in list1:
    j+=1
    i=int(i)
    outfile.write(str(j)+'\t'+str(i)+'\n')
outfile.close()
rawdata.close()


