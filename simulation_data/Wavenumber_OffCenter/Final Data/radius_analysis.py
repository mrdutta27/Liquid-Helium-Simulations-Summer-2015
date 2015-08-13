# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')

#Prompt Signal Rate vs All Signal Rate
radius, time = [],[]

wavenumber=float(2.2)
y0=float(0)

for row in data: 
    #if float(row[1]) == float(depth) and float(row[0])==float(wavenumber): #phonon at depth 0
    if float(row[1]) == float(y0): #phonon at depth 0
        radius.append((float(row[3])**2.0+float(row[4])**2.0))
        time.append(float(row[2]))
        
y,binEdges=np.histogram(radius, range=[0,.0025],bins=41)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters,
    y,
    yerr = y**0.5,
    marker = '.',
    drawstyle = 'steps-mid-'
)

plt.ylabel('# of Events')
plt.xlabel('Radius^2')
plt.title("Radius^2 Heatmap: y0=0")
plt.savefig("C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\plots\heatmaps\Radius^2 Heatmap" + str(y0)+".jpg")
#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()
plt.clf()

print y