#!/usr/bin/env python
#file name:fast_bsc.py

from numpy import *

#fast_file="/home/niutj/python_pros/test_QC.fq"
fast_file="/1_disk/public_resources/jiankuihe-exome-1.fq"
rawdata=open(fast_file,"r")

lnum=0
rawdata.readline()
lista=zeros((100))
listt=zeros((100))
listg=zeros((100))
listc=zeros((100))

for line in rawdata:
    lnum+=1
    line=line.strip('\n')
    if(lnum%4==1):
        line.upper()
        count=0
        for char in line:
            if char=='A':
                lista[count]+=1

            if char=='T':
                listt[count]+=1

            if char=='G':
                listg[count]+=1

            if char=='C':
                listc[count]+=1

            count+=1

list_1=add(lista,listt)
list_2=add(listg,listc)
list_sum=add(list_1,list_2)

l=0
outfile=open('exome_1_3.txt','w')
outfile.write('num\ta\tt\tg\tc\tsum\n')
for h in lista:
    l+=1
    for i in listt:
        for j in listg:
            for k in listc:
                for m in list_sum:
                    outfile.write(str(l)+'\t'+str(h)+'\t'+str(i)+'\t'+str(j)+'\t'+str(k)+'\t'+str(m)+'\n')
                    break
                break
            break
        break
    
outfile.close()
rawdata.close()
