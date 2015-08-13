# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')

#Prompt Signal Rate vs All Signal Rate
x_pos, z_pos, time = [],[],[]

wavenumber=float(.7)
y0=float(-0.025)
for row in data: 
    if float(row[1]) == float(y0) and .000<float(row[2])<.001 and float(row[10])==1: #and float(row[0])==float(wavenumber): #phonon at depth 0
        x_pos.append(float(row[3]))
        z_pos.append(float(row[4]))
        time.append(float(row[2]))
        
hist, xedges, yedges = np.histogram2d(x_pos, y_pos, range=[[-.05,.05],[-.05,.05]], bins=21) 


import plotly.plotly as py
from plotly.graph_objs import *

data = Data([
    Heatmap(
        z=hist
    )
])
#plot_url = py.plot(data, filename='basic-heatmap')
py.image.save_as({'data': data}, 'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\plots\heatmaps\Time Lapse 0.025\heatmap_y0_0.025_0.001_noprompt.png')
