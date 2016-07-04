#!/usr/bin/env python
#file name:fast_1.py
# Compute per base sequence quality of fastq files
#

import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import linecache

fastq_file="sample1.fq"
# Compute line number
#length = len(linecache.getline(fastq_file,4))-1
#print(length)
FILE=open(fastq_file,"r")
lnum = 0

for line in FILE:
    lnum += 1
    line=line.strip('\n')

    if lnum == 2:
        #print(line+'\n')
        seq_length = len(line)
        bq_ls = np.zeros((seq_length,92)) # initialize array
        print('sequence length: '+ str(len(line)))

    # compute per base sequence quality
    if lnum % 4 == 0:
        #print(line)
        count = 0
        for char in line:
            bq = ord(char) - 33 -2
            bq_ls[count][bq]+=1
            count +=1
print('lnum: '+str(lnum))
TEM = open('tem.txt','w')
#box = np.array
ls = []
box_ls =[]
for i in bq_ls:
    #TEM.write(i)
    score_sum = 0
    sub = 2
    for j in i: # len(i) == 92
        #print(j)
        score_sum += j*sub
        sub +=1
    ls.append(score_sum*4/lnum)
    count = 0
    tem_ls = []
    #tem_ls.append(2)
    #tem_ls.append(93)
    for k in range(92):
        count += i[k]
        if count <= 0.1*lnum/4 and count+(i[k+1]) >= 0.1*lnum/4:
           #print('hi')
            tem_ls.append(k+2)
        if count <= 0.25*lnum/4 and count+(i[k+1]) >= 0.25*lnum/4:
            tem_ls.append(k+2)
            tem_ls.append(k+2)
        if count <= 0.5*lnum/4 and count+(i[k+1]) >= 0.5*lnum/4:
            tem_ls.append(k+2)
        if count <= 0.75*lnum/4 and count+(i[k+1]) >= 0.75*lnum/4:
            tem_ls.append(k+2)
            tem_ls.append(k+2)
        if count <= 0.9*lnum/4 and count+(i[k+1]) >= 0.9*lnum/4:
            #print('hi')
            tem_ls.append(k+2)
    box_ls.append(tem_ls)
#print(box_ls)
fs =10
plt.boxplot(box_ls)
plt.title('Test',fontsize = fs)
plt.axis([0,150, 0,50])
plt.savefig('box_test.jpg')

#print(ls[:20])

print('ls length: '+str(len(ls)))
x = np.linspace(1,seq_length,seq_length)
y = ls
plt.title('Quality scores across all bases', size = 14)
plt.xlabel('Position in read (bp)', size = 14)
#xtick.major.size(15)
#xtick.minor.size(1)
#ystick.major.size(10)
#ystick.minor.size(2)
plt.boxplot(x, y,box_ls, color = 'b', linestyle = '-', label = 'Mean Quality')
plt.axis([1, seq_length, 0, 50])

plt.legend(loc = 'upper right')
plt.savefig('perBaseSeqQua.jpg', format='jpg')
        
