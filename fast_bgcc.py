from __future__ import division
#!/usr/bin/evn python
#file name: fast_bgcc.py
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

    if(lnum % 4 == 1):
        line.upper()
        count=0
        for char in line:
            if char=='G' or char=='C':
                list1[count]+=1
            count+=1

lnum1=lnum/4
j=0
outfile=open('exome_1_4.txt','w')
outfile.write('num\tgccon\n')
for i in list1:
    j+=1
    con=i/lnum1
    outfile.write(str(j)+'\t'+str(con)+'\n')
outfile.close()
rawdata.close()
