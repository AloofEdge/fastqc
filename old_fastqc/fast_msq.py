#!/usr/bin/env python
#file name:fast_msq.py
from numpy import *

#fast_file="/home/niutj/python_pros/test_QC.fq"
fast_file="/1_disk/public_resources/jiankuihe-exome-1.fq"
rawdata=open(fast_file,"r")

lnum=0
rawdata.readline()
list1=zeros((80))

for line in rawdata:
    lnum+=1
    line=line.strip('\n')
    #print(line)
    allval=0
    if(lnum%4==3):
        for char in line:
            ascval=ord(char)-33
            allval+=ascval
        val=allval/100
        #print('val:',val)
        list1[val]+=1
j=0
outfile=open('exome_1_2.txt','w')
outfile.write('num\tqua\n')
for i in list1:
    i=int(i)
    j+=1
    #outfile.write(str(j))
    outfile.write(str(j)+'\t'+str(i)+'\n')
    #print(j,i)
    #print >> outfile, j i

outfile.close()
rawdata.close()
