#!/usr/bin/env python
# fn: box_plot.py
# plot boxplot and save it.

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# fake data
data = [[1, 2, 3, 4, 5],[1, 1.5, 2.5, 3, 3.5, 4.5, 5],[31, 34, 34, 37, 39, 39, 40]]
#labels = ['1','2']
fs = 10 # fontsize

plt.boxplot(data, whis = [0, 100] )
plt.title('Default', fontsize = fs)
plt.axis([0,5, 0, 6])
plt.savefig('box2.jpg')
