#!/usr/bin/env python
#file name:fast_1.py
from numpy import *

#fast_file="/home/niutj/python_pros/test_QC.fq"
fast_file="/1_disk/public_resources/jiankuihe-exome-1.fq"
rawdata=open(fast_file,"r")
count=0

list1=zeros((100))
list2=zeros((100))
list3=([2])
lnum=0
rawdata.readline()

for line in rawdata:
    lnum +=1
    #print lnum;
    line=line.strip('\n')
    if(lnum%4==3):
        #print (line)
        #print(lnum)
        #print(line[0])
#        ls_line=line.split()
        list0=[]
        for char in line:
            #print (char)
            #print(ord(char))
            ascval=ord(char)-33
            #print(ascval)
            list0.append(ascval)
            '''
        if(len(list0)==100):
                #print(line)
            count+=1
            '''
        list2=add(list0,list1)
        list1=divide(list2,list3)
#print(count)
#print(list1)

outfile=open('exome_1_1.txt','w')
for i in list1:
    i=int(i)
    print >> outfile, i
#print >>outfile, list1
outfile.close()

rawdata.close()




