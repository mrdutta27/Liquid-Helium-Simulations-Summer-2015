from matplotlib.pylab import *
from numpy import *

#initial conditions
v0=raw_input("Initial Velocity= ")
v0=float(v0)
y_pos=raw_input("Initial Height= ")
y_pos=float(y_pos)
alpha=raw_input("Angle Thrown from Horizontal(Degrees)= ")
alpha=float(alpha)*(pi/180)
v0x=v0*cos(alpha)
v0y=v0*sin(alpha)
accel=-9.8065
bbox=raw_input("Bounding Box Size= ")
bbox=float(bbox)
sim_time=raw_input("Simulation Duration= ")
sim_time=float(sim_time)

#array initialization
t=linspace(0,sim_time,sim_time*1000) #10000 points per second
pos_y = zeros(len(t)) # allocate Y positions with float elements
vel_y = zeros(len(t)) # allocate Y velocities with float elements
pos_x = zeros(len(t)) # allocate X positions with float elements
vel_x = zeros(len(t)) # allocate X velocities with float elements

#time definitions
t_bounce_x=bbox/(v0x)
t1_bounce_y=(-v0y-sqrt(v0y**2-2*accel*y_pos))/accel
t2_bounce_y=(v0y-sqrt(v0y**2-2*accel*y_pos))/accel
t_bounce_y=t1_bounce_y+t2_bounce_y
del_t=t[2]-t[1]
pos_y[0]=y_pos

#filling array values
for i in xrange(len(t)-1): 
    # X Direction
    if (t[i]%(2*t_bounce_x))<=t_bounce_x:
        vel_x[i]=v0x
    else:
        vel_x[i]=-v0x
    pos_x[i+1]=pos_x[i]+del_t*vel_x[i]
    # Y Direction
    vel_y[i]=-0.5*accel*t_bounce_y+accel*((t[i]+t2_bounce_y)%t_bounce_y)
    pos_y[i+1]=pos_y[i]+del_t*vel_y[i]
    
plot(pos_x,pos_y)
show()
