import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\bad data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')

#Prompt Signal Rate vs All Signal Rate
time, time_prompt = [], [] 
bin_count=200
range_count=.05
for row in data: 
    if float(row[1]) == -.45 and 1.1< float(row[0])< 2.0: #phonon at depth 0
        time.append(float(row[2]))
        if float(row[10]) == 0: #prompt
            time_prompt.append(float(row[2]))
            

y,binEdges=np.histogram(time, range=[0,range_count],bins=bin_count)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
z,binEdges2=np.histogram(time_prompt, range=[0,range_count],bins=bin_count)
bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
p.plot(bincenters,y,'-',label="All")
p.plot(bincenters2,z,'-',label="Prompt")
p.legend()
p.show()

#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()