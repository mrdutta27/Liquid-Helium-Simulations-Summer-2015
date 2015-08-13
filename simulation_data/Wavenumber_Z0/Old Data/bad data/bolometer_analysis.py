# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\small data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')
data2 = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\small data\transmission\all_transmission.csv', 'rb'), delimiter=",", quotechar='|')


#Prompt Signal Rate vs All Signal Rate
time, time_prompt, time_transmission = [], [], []
bin_count=50
range_count=1
wavenumber=float(2.2)
depth=float(0)
for row in data: 
    if float(row[1]) == float(depth) and float(row[0])==float(wavenumber): 
        time.append(float(row[2]))
        if float(row[10]) == 0: #prompt
            time_prompt.append(float(row[2]))
for row in data2: 
    if float(row[1]) == float(depth) and float(row[0])==float(wavenumber):
        time_transmission.append(float(row[2]))
                  

y,binEdges=np.histogram(time, range=[0,range_count],bins=bin_count)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters,
    y,
    yerr = y**0.5,
    marker = '.',
    drawstyle = 'steps-mid-'
)
z,binEdges2=np.histogram(time_prompt, range=[0,range_count],bins=bin_count)
bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
plt.errorbar(
    bincenters2,
    z,
    yerr = z**0.5,
    marker = '.',
    drawstyle = 'steps-mid-'
)
t,binEdges3=np.histogram(time_transmission, range=[0,range_count],bins=bin_count)
bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])
plt.errorbar(
    bincenters3,
    t,
    yerr = t**0.5,
    marker = '.',
    drawstyle = 'steps-mid-'
)

plt.ylabel('# of Events')
plt.xlabel('Time (Seconds)')
plt.legend(["Bolometer","Bolometer Prompt","Transmission"])
plt.title("Z="+str(depth)+", K="+str(wavenumber))
plt.savefig("C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\small data\plots\Prompt Large Scale\prompt" + str(wavenumber) + str(depth)+".jpg")
#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()
plt.clf()