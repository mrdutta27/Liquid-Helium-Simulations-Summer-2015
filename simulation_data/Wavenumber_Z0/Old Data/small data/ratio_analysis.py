import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')
data2 = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber_OffCenter\Final Data\transmission\all_transmission.csv', 'rb'), delimiter=",", quotechar='|')

transmission_count=0
absorption_count=0
wavenumber=2.2

for row in data: 
    if float(row[0])==float(wavenumber):
        absorption_count+=1
for row in data2: 
    if float(row[0])==float(wavenumber):
        transmission_count+=1 
print str(wavenumber) + ", " + str((float(absorption_count)/float(transmission_count)))