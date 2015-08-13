import csv
import matplotlib.pyplot as plt
import pylab as p
import numpy as np

data = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\bad data\bolometer\all_bolometer.csv', 'rb'), delimiter=",", quotechar='|')
data2 = csv.reader(open(r'C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\SimulationData\Wavenumber\bad data\transmission\all_transmission.csv', 'rb'), delimiter=",", quotechar='|')
ratio = []
transmission_count=0
absorption_count=0
for i in range(5,23):
    for row in data: 
        if float(row[0])==float(i)/10.0:
            absorption_count+=1
    for row in data2: 
        if float(row[0])==float(i)/10.0:
            transmission_count+=1 
    ratio.append(float(absorption_count))
    transmission_count=0
    absorption_count=0

print ratio