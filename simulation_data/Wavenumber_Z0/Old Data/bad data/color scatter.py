# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\Particle Sim\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')

#Prompt Signal Rate vs All Signal Rate
x_pos, y_pos, time = [],[],[]

wavenumber=float(2.2)
depth=float(.45)

for row in data: 
    #if float(row[1]) == float(depth) and float(row[0])==float(wavenumber): #phonon at depth 0
    if float(row[1]) == float(depth): #phonon at depth 0
        x_pos.append(float(row[3]))
        y_pos.append(float(row[4]))
        time.append(float(row[2]))
        
hist, xedges, yedges = np.histogram2d(x_pos, y_pos, range=[[-.5,.5],[-.5,.5]], bins=21) 


import plotly.plotly as py
from plotly.graph_objs import *

data = Data([
    Heatmap(
        z=hist
    )
])

py.image.save_as({'data': data}, r"C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\Particle Sim\plots\Heatmaps\Heatmap_z.45.png")

print x_pos