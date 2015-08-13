from matplotlib.pylab import *
from numpy import *

with open(r"C:\Users\Suryabrata\Dropbox\Yale\Summer 2015\McKinsey Research\Python\roton.txt") as textFile:
    lines = [line.split() for line in textFile]

idx=raw_input("Index(Temp and Pressure)= ")
idx=int(idx)

    
p0=.00001264 #Critical Momentum: 1.264 * 10^-5 ev*s/m
mu=621.3669 #m4/6
k_b = 8.62*(1/100000.0) #boltzmann's constant in ev/k
delta=8.62 #in k
hbar=4.135667516*(1/(10.0**15))
status=False

def energy(p):
    return float(lines[idx][2])*p + float(lines[idx][3])*p**3 + float(lines[idx][4])*p**4 + float(lines[idx][5])*p**5 + float(lines[idx][6])*p**6 + float(lines[idx][7])*p**7 +float(lines[idx][8])*p**8
 
def velocity(p,status):
    if status==False:
        return float(lines[idx][2]) + 3*float(lines[idx][3])*p**2 + 4*float(lines[idx][4])*p**3 + 5*float(lines[idx][5])*p**4 + 6*float(lines[idx][6])*p**5 + 7*float(lines[idx][7])*p**6 + 8*float(lines[idx][8])*p**7 
    else:
        return float(lines[idx][2])
    
p = linspace(0, 2.5, 10000) # 31 points between 0 and 3
e = zeros(len(p)) # allocate y with float elements
v = zeros(len(p)) # allocate y with float elements
p_squared = zeros(len(p)) # allocate y with float elements
e_squared = zeros(len(p)) # allocate y with float elements
m = zeros(len(p))
del_p=p[2]-p[1]

for i in xrange(len(p)):
    
    if p[i]>1:
        if float(velocity(p[i],status))>float(lines[idx][2]):
            status=True
            print p[i]
    v[i] = velocity(p[i],status)
    e[i] = e[i-1]+del_p*v[i]
    m[i] = p[i]/v[i]
    if e[i-1]<e[i] and e[i-1]<e[i-2]:
        print p[i]
        

plot(p,e)
title("energy and velocity relation")
xlabel("energy")
ylabel("group velocity")
show()

for i in xrange(0,10000,4000):
    print e[i]