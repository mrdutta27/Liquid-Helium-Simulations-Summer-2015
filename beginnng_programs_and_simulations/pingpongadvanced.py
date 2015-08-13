# -*- coding: utf-8 -*-
from numpy import *
from mayavi.mlab import *

#initial conditions
n_pong=raw_input("Number of Ping Pong Balls = ")
n_pong=int(n_pong)
v0=raw_input("Initial Velocity = ")
v0=float(v0)
theta=raw_input("Theta(Degrees) = ")
theta=float(theta)*(pi/180)
phi=raw_input("Phi(Degrees) = ")
phi=float(phi)*(pi/180)
bbox=raw_input("Bounding Box Dimension = ")
bbox=float(bbox)
sim_time=raw_input("Simulation Duration = ")
sim_time=float(sim_time)

#array initialization
t=linspace(0,sim_time,sim_time*1000) #1000 points per second
a=ones((n_pong, 7, len(t))) #3-D array with N rows, 
del_t=t[2]-t[1]

#initial column values
"""
columns are 
0: # of ball
1: x-pos 
2: x-vel 
3: y-pos 
4: y-vel
5: z-pos
6 z-vel
"""
def x_pos(ball_n,time):
    return a[ball_n,1,time]
def x_vel(ball_n,time):
    return a[ball_n,2,time]
def y_pos(ball_n,time):
    return a[ball_n,3,time]
def y_vel(ball_n,time):
    return a[ball_n,4,time]
def z_pos(ball_n,time):
    return a[ball_n,5,time]
def z_vel(ball_n,time):
    return a[ball_n,6,time]


for i in xrange(n_pong):
    v0_rand=random.normal(loc=v0) #gaussian random speed
    theta_rand=random.normal(loc=theta, scale=pi) #gaussian random theta
    phi_rand=random.normal(loc=phi, scale=pi) #gaussian random phi
    a[i,0]=i+1 #fill in column 0
    a[i,1,0]=0 #initial x position
    a[i,2,0]=v0_rand*sin(phi_rand)*cos(theta_rand)
    a[i,3,0]=0 #initial y position
    a[i,4,0]=v0_rand*sin(phi_rand)*sin(theta_rand)
    a[i,5,0]=0 #initial z position
    a[i,6,0]=v0_rand*cos(phi_rand)
        
#filling array values
for j in xrange(1,len(t)-1): #iterating time values
    for i in xrange(n_pong): #iterating ping pong values
        # X Direction
        a[i,2,j]=a[i,2,j-1] 
        a[i,1,j+1]=x_pos(i,j)+del_t*x_vel(i,j)
        if a[i,1,j+1]>bbox or a[i,1,j+1]<-bbox: #x boundary condition
            a[i,2,j]=-1*a[i,2,j]  
        a[i,1,j+1]=x_pos(i,j)+del_t*x_vel(i,j)
        # Y Direction
        a[i,4,j]=a[i,4,j-1] 
        if a[i,3,j]>bbox or a[i,3,j]<-bbox: #y boundary condition
            a[i,4,j]=-1*a[i,4,j]      
        a[i,3,j+1]=y_pos(i,j)+del_t*y_vel(i,j) 
        # Z Direction
        a[i,6,j]=a[i,6,j-1] 
        if a[i,5,j]>bbox or a[i,5,j]<-bbox: #z boundary condition
            a[i,6,j]=-1*a[i,6,j] 
        a[i,5,j+1]=z_pos(i,j)+del_t*z_vel(i,j)
        
for i in xrange(n_pong):
    color=float(i)/n_pong
    plot3d(a[i,1],a[i,3],a[i,5],color=(0,0,color),tube_radius=0.125)