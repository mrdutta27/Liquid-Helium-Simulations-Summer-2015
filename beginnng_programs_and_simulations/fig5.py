from numpy import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

with open("roton.txt") as textFile:
    lines = [line.split() for line in textFile]

idx=raw_input("Index(Temp and Pressure)= ")
idx=int(idx)

def array_maker(omega):
    array_omega=zeros(32)
    array_omega[31]=-1.0*omega**2
    return array_omega
    
theta=3*pi/8

def rootsolver(array_omega):
    disp_curve=poly1d([float(lines[idx][8]),0.0,float(lines[idx][7]),0.0,float(lines[idx][6]),0.0,float(lines[idx][5]),0.0,float(lines[idx][4]),0.0,float(lines[idx][3]),0.0,0.0,0.0,float(lines[idx][2]),0.0,0.0])
    disp_curve_squared=disp_curve*disp_curve + poly1d(array_omega)
    disp_curve_squared_roots=disp_curve_squared.r
    real_roots=isreal(disp_curve_squared_roots)
    final_roots=zeros(1)
    for i in xrange(len(disp_curve_squared_roots)):
        if real_roots[i]==True and disp_curve_squared_roots[i]>0:
            a=array([disp_curve_squared_roots[i]])
            final_roots=sort(append(final_roots,a))
    return array([float(real(final_roots[1])),-1.0*float(real(final_roots[2])), float(real(final_roots[3]))])
    
def z(k):
    return k*cos(theta)
    
def refl(i,j,k):
    return (abs(4*z(i)*z(j))/((z(i)-z(j))**2))*abs(((z(i)+z(k))*(z(j)+z(k)))/((z(i)-z(k))*(z(j)-z(k))))*abs((i+j+k-(z(i)+z(j)))/(i+j+k))**2

def refl2(i,j,k):
    return abs(((z(i)+z(j))*(z(i)+z(k)))/((z(i)-z(j))*(z(i)-z(k))))**2*abs((i+j+k-2*z(i))/(i+j+k))**2
    
energy=linspace(8.5,14.5,350)
angles=linspace(0,pi/2,250)
refl_coef_12=zeros((len(energy),len(angles)))
refl_coef_23=zeros((len(energy),len(angles)))
refl_coef_13=zeros((len(energy),len(angles)))
refl_coef_11=zeros((len(energy),len(angles)))
refl_coef_22=zeros((len(energy),len(angles)))
refl_coef_33=zeros((len(energy),len(angles)))

for i in xrange(len(energy)):
    for j in xrange(len(angles)):
        theta=angles[j]
        k=rootsolver(array_maker(energy[i]))
        refl_coef_12[i,j]=refl(k[0],k[1],k[2])
        refl_coef_23[i,j]=refl(k[1],k[2],k[0])
        refl_coef_13[i,j]=refl(k[0],k[2],k[1])
        refl_coef_11[i,j]=refl2(k[0],k[1],k[2])
        refl_coef_22[i,j]=refl2(k[1],k[2],k[0])
        refl_coef_33[i,j]=refl2(k[2],k[0],k[1])
        
def reflector12(i,j):
    i=int(i)
    j=int(j)
    return refl_coef_23[i,j]

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def fun(x, y):
  return x + y

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 10
xs = [i for i in range(len(energy)) for _ in range(len(angles))]
ys = range(len(angles))*len(energy)
x2s=[(xs[i]/50.0)+8.5 for i in range(len(xs))]
y2s=[ys[i]*(pi/500) for i in range(len(xs))]
zs = [reflector12(x, y) for x,y in zip(xs,ys)]


ax.scatter(x2s, y2s, zs,s=.5)

ax.set_xlabel('Energy')
ax.set_ylabel('Angles')
ax.set_zlabel('Reflection Coefficient 23')

plt.title("23 Reflection Coefficients")

plt.show()