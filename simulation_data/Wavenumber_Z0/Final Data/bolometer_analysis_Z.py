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
wavenumber=float(2.2)
z0=float(0.045)
mass_4He=float(6.646478)*10**(-27)
quas_scale_factor_array=([[0.955405029142,1.16736603371,1.27286665057,1.69892723614,1.81554714974,1.85405693655,0.692011623406,1.85160570536,1.8363215644,1.6682221161,1.29237372069,1.29464481074,1.20687773426,0.986676220476,0.918338092856,0.94826458046,1.28013891291,1.49577442413],[1.04339485096,1.26798837521,1.27641717357,1.69436763932,1.81774848002,1.86102990937,0.70534885183,1.85482720746,1.83465638616,1.67077534623,1.26402613194,1.33305098802,1.26663441032,0.995499545481,0.914632864674,0.93817675079,1.08816446883,1.26999325566],[1.08927735425,1.30212573031,1.29525277549,1.7020711679,1.8162809265,1.85688381742,0.693726409917,1.85748020918,1.83373128713,1.66986782141,1.24749044125,1.30722604454,1.24689578974,0.9702385393,0.91245379221,0.922040139124,1.08696745452,1.24254807235]]
) #phonon, roton R+, roton R- | 0.045,0,-0.045
quas_scale_factor=float(1.60217646*(10**-19))/float(quas_scale_factor_array[0][int(wavenumber*10-5)]*(10**-20))

lines=[0.00,0.00,18.20,45.11,-163.09,211.67,-136.12,43.27,-5.39]
disp_curve=np.poly1d([float(lines[8]),float(lines[7]),float(lines[6]),float(lines[5]),float(lines[4]),float(lines[3]),0.0,float(lines[2]),0])#dispersion curve

for row in data: 
    if float(row[1]) == float(z0) and float(row[0])==wavenumber: 
        time.append(float(row[2])*1000.0)
        energy_weights.append((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))*quas_scale_factor/(range_count*1000/bin_count))
        energy_weights2.append(((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))+(1.441959*(10**-21)))*quas_scale_factor/(range_count*1000/bin_count))
        if float(row[10]) == 0: #prompt
            time_prompt.append(float(row[2])*1000.0)
            energy_weights_prompt.append((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))*quas_scale_factor/(range_count*1000/bin_count))
            energy_weights_prompt2.append(((.5*mass_4He*(float(row[6])**2+float(row[7])**2+float(row[8])**2))+(1.441959*(10**-21)))*quas_scale_factor/(range_count*1000/bin_count))
for row in data2: 
    if float(row[1]) == float(z0) and float(row[0])==wavenumber:
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
plt.title("Z0="+str(z0)+", K="+str(wavenumber))
plt.ylim([0,1*10**-19])
plt.savefig("C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_Z0\Final Data\plots\Z Values\Energy Small Scale\Energy" + str(wavenumber)+str(z0)+".jpg")
#plt.hist(time, range=[0,range_count],bins=bin_count, color='#000000', alpha=0.5)
#plt.hist(time_prompt, range=[0,range_count],bins=bin_count, color='#0000aa', alpha=0.5)
#plt.show()
plt.clf()
