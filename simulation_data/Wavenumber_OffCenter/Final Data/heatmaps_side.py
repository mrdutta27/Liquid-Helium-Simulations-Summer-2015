# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\transmission\all_transmission.csv', 'rb'), delimiter=",", quotechar='|')

#Prompt Signal Rate vs All Signal Rate
z, time = [],[]

z0=float(0)

for row in data: 
    if float(row[1]) == float(z0) and float(row[5])>-.05: #phonon at depth 0
        z.append((float(row[5])))
        time.append(float(row[2]))
        
y,binEdges=np.histogram(z, range=[-.05,.05],bins=41)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters,
    y,
    yerr = y**0.5,
    marker = '.',
    drawstyle = 'steps-mid-'
)

plt.ylabel('# of Events')
plt.xlabel('Z')
plt.title("Side Z Heatmap: z0="+str(z0))
plt.savefig("C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\plots\heatmaps\Side Z Heatmap" + str(z0)+".jpg")
#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()
plt.clf()