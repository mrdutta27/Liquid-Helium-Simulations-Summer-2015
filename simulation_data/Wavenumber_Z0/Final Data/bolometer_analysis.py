# -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')
data2 = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\transmission\all_transmission.csv', 'rb'), delimiter=",", quotechar='|')
#Prompt Signal Rate vs All Signal Rate
time, time_prompt, time_transmission, energy_weights, energy_weights_prompt, energy_weights_transmission = [], [], [], [], [], []
energy_weights2, energy_weights_prompt2 = [],[]
bin_count=50
range_count=2
#wavenumber=float(2.2)
z0=float(0)
mass_4He=float(6.646478)*10**(-27)
quas_scale_factor_array=([[.945618065925,.966629528029,.975561818179],[.37241779175,.329633447528,.325155566599],[1.10550599649,1.11341028803,1.10453839248]]) #phonon, roton R+, roton R- | 0.045,0,-0.045
quas_scale_factor=float(1.60217646*(10**-19))/float(quas_scale_factor_array[0][1]*(10**-20))

lines=[0.00,0.00,18.20,45.11,-163.09,211.67,-136.12,43.27,-5.39]
disp_curve=np.poly1d([float(lines[8]),float(lines[7]),float(lines[6]),float(lines[5]),float(lines[4]),float(lines[3]),0.0,float(lines[2]),0])#dispersion curve

for row in data: 
    if float(row[1]) == float(z0) and 1.2 > float(row[0]): 
        time.append(float(row[2])*1000.0)
        energy_weights.append((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))*quas_scale_factor/(range_count*1000/bin_count))
        energy_weights2.append(((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))+(1.441959*(10**-21)))*quas_scale_factor/(range_count*1000/bin_count))
        if float(row[10]) == 0: #prompt
            time_prompt.append(float(row[2])*1000.0)
            energy_weights_prompt.append((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))*quas_scale_factor/(range_count*1000/bin_count))
            energy_weights_prompt2.append(((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))+(1.441959*(10**-21)))*quas_scale_factor/(range_count*1000/bin_count))
for row in data2: 
    if float(row[1]) == float(z0) and 1.2 > float(row[0]):
        time_transmission.append(float(row[2])*1000.0)
        energy_weights_transmission.append((1.38065*(10**-23))*disp_curve(float(row[0]))*quas_scale_factor/(range_count*1000/bin_count))
                  
x,binEdges=np.histogram(time, range=[0,range_count],bins=bin_count, weights=energy_weights2)#quantumevaporation+BE
bincenters4 = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters4,
    x,
    #yerr = 0#(y**0.5)/quas_number,
    marker = '.',
    drawstyle = 'steps-mid-',
    color="b",
    linestyle = "--"
)
y,binEdges=np.histogram(time, range=[0,range_count],bins=bin_count, weights=energy_weights)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters,
    y,
    #yerr = 0#(y**0.5)/quas_number,
    marker = '.',
    drawstyle = 'steps-mid-',
    color="b"
)
p,binEdges=np.histogram(time_prompt, range=[0,range_count],bins=bin_count, weights=energy_weights_prompt2)#quantumevaporation+BE
bincenters5 = 0.5*(binEdges[1:]+binEdges[:-1])
plt.errorbar(
    bincenters5,
    p,
    #yerr = 0#(y**0.5)/quas_number,
    marker = '.',
    drawstyle = 'steps-mid-',
    color="g",
    linestyle = "--"
)

z,binEdges2=np.histogram(time_prompt, range=[0,range_count],bins=bin_count, weights=energy_weights_prompt)
bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
plt.errorbar(
    bincenters2,
    z,
    #yerr = 0#(y**0.5)/quas_number,
    marker = '.',
    drawstyle = 'steps-mid-',
    color="g",
)
t,binEdges3=np.histogram(time_transmission, range=[0,range_count],bins=bin_count, weights=energy_weights_transmission)
bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])
plt.errorbar(
    bincenters3,
    t,
    #yerr = 0#(t**0.5)/quas_number,
    marker = '.',
    drawstyle = 'steps-mid-',
    color="r"
)

plt.ylabel('Power (Watts)')
plt.xlabel('Time (Milliseconds)')
plt.legend(["QE + BE","QE","QE Prompt + BE","QE Prompt","Wall Absorption"])
plt.title("Z0="+str(z0)+", Phonon")
plt.ylim([0,.75*10**-18])
plt.savefig("C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\plots\Categories\Energy Small Scale\EnergyPhonon" + str(z0)+".jpg")
#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()
plt.clf()
