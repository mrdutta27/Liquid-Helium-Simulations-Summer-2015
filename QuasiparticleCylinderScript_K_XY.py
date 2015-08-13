from numpy import *
#Definitions

def z(k):
    return k*cos(theta_i)
        
def refl(i,j,k): #reflection coefficient of R(ij)
    return (abs(4*z(i)*z(j))/((z(i)-z(j))**2))*abs(((z(i)+z(k))*(z(j)+z(k)))/((z(i)-z(k))*(z(j)-z(k))))*abs((i+j+k-(z(i)+z(j)))/(i+j+k))**2

def refl2(i,j,k): #reflection coefficient of R(ii)
    return abs(((z(i)+z(j))*(z(i)+z(k)))/((z(i)-z(j))*(z(i)-z(k))))**2*abs((i+j+k-2*z(i))/(i+j+k))**2    
    
def neg_array(array1): #turning an array negative
    array2=zeros(len(array1))
    for i in xrange(len(array1)):
        array2[i]=-1*array1[i]
    return array2
    
def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return sqrt(dotproduct(v, v))

def angle(v1, v2):
    if (length(v1) * length(v2) == 0):
        return 0
    return arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))
        
def newstater(state1):
    state1=int(state1)
    prob_state1=refl2(k[state1],k[(state1+1)%3],k[(state1+2)%3])
    prob_state2=refl(k[state1],k[(state1+1)%3],k[(state1+2)%3])
    prob_state3=refl(k[state1],k[(state1+2)%3],k[(state1+1)%3])
    rand=random.random_sample()
    if rand<prob_state1:
        return state1
    elif rand>=prob_state1 and rand<(prob_state1+prob_state2):
        return (state1+1)%3
    elif rand>=(prob_state1+prob_state2) and rand<(prob_state1+prob_state2+prob_state3):
        return (state1+2)%3
    else:
        return 3
        
def newstater2(state1):
    state1=int(state1)
    rand=random.random_sample()
    if rand<QEprob(final_zeroes[state1+1]):
        return 3
    else:
        return state1

def bounce(vel1,normal,state1,state2):
    state1=int(state1)
    state2=int(state2)
    ratio=(velocities[state2]*k[state1])/(velocities[state1]*k[state2])
    mag_normal=linalg.norm(normal)
    unit_normal=divide(normal,mag_normal)
    vel1_proj_normal=multiply(unit_normal, inner(neg_array(vel1),unit_normal))
    vel1_proj_t=add(vel1, vel1_proj_normal)
    vel2_proj_t=multiply(vel1_proj_t,ratio)
    vel2_proj_normal=multiply(unit_normal, sqrt(abs((vel(state2)**2)-(linalg.norm(vel2_proj_t)**2))))
    vel2=add(vel2_proj_t,vel2_proj_normal)
    return vel2

def quantumevaporate(vel1,normal):
    ratio=(final_zeroes[a[i,7,j-1]+1]/k_a)
    mag_normal=linalg.norm(normal)
    unit_normal=divide(normal,mag_normal)
    vel1_proj_normal=multiply(unit_normal, inner(neg_array(vel1),unit_normal))
    vel1_proj_t=add(vel1, vel1_proj_normal)
    vel2_proj_t=multiply(vel1_proj_t,ratio)
    vel2_proj_normal=multiply(neg_array(unit_normal), sqrt(abs(((h_bar_velocity*k_a/mass_helium)**2)-(linalg.norm(vel2_proj_t)**2))))
    vel2=add(vel2_proj_t,vel2_proj_normal)
    return vel2

def QEprob(wavenumber_current):
    with open("/global/homes/e/epease/quantumevaporationprobabilities.txt") as textFile:
        QEprobs = [line.split() for line in textFile]
    for i in range(1610):
        if wavenumber_current<float(QEprobs[i][0]):
            return float(QEprobs[i][1])

def vel(state):
        return velocities[state]*13.1 #get better value
    
def array_maker(omega): #array of zeroes plus omega for final polynomial
    array_omega=zeros(32)
    array_omega[31]=-1.0*omega**2
    return array_omega
    
def rootsolver(array_omega): #solving the roots for omega^2(k^2)
    disp_curve_2=poly1d([float(lines[8]),0.0,float(lines[7]),0.0,float(lines[6]),0.0,float(lines[5]),0.0,float(lines[4]),0.0,float(lines[3]),0.0,0.0,0.0,float(lines[2]),0.0,0.0])
    disp_curve_squared=disp_curve_2*disp_curve_2 + poly1d(array_omega)
    disp_curve_squared_roots=disp_curve_squared.r
    real_roots=isreal(disp_curve_squared_roots)
    final_roots=zeros(1)
    for i in xrange(len(disp_curve_squared_roots)):
        if real_roots[i]==True and disp_curve_squared_roots[i]>0:
            a=array([real(disp_curve_squared_roots[i])])
            final_roots=sort(append(final_roots,a))
    return array([float(final_roots[1]),-1.0*float(final_roots[2]), float(final_roots[3])])
    
    
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
    
#____________________________________________________________________

#Constants
mass_helium=float(2.658591*(10**-26))
reduced_plancks_constant_squared=float(8.0551*(10**-26))
bindingenergy=float(7.16)
h_bar_velocity=float(1.054572*10**-24)
n_particle=int(10)
#n_particle=int(10000)
x0 = float(.025)
z0 = float(0)
cyl_rad=float(0.05)
cyl_height=float(.1)
sim_time=float(.5)
vacuum_gap=float(0.005)


#___________________________________________________________________

#array initialization
t=linspace(0,sim_time,sim_time*1000000) #1000000 points per second
transmission=[]
bolometer=[]
top=[]
atom_side=[]
del_t=t[2]-t[1]

#__________________________________________________________________
randloopnumber=int(random.rand()*10**16)
for zz in xrange(5,23): #iterating over wavenumbers
        #print "Currently Processing Wavenumber: " + str(float(zz)/10)
        #initial conditions
        wavenumber = float(zz)/10
        y0 = float(.025)
    
        #dispersion curve analysis
        lines=[0.00,0.00,18.20,45.11,-163.09,211.67,-136.12,43.27,-5.39]
        disp_curve=poly1d([float(lines[8]),float(lines[7]),float(lines[6]),float(lines[5]),float(lines[4]),float(lines[3]),0.0,float(lines[2]),0])#dispersion curve
        omega=disp_curve(wavenumber)
        disp_curve_minusomega=poly1d([float(lines[8]),float(lines[7]),float(lines[6]),float(lines[5]),float(lines[4]),float(lines[3]),0.0,float(lines[2]),-1.0*omega]) #dispersion curve minus omega
        disp_curve_diff=poly1d([8.0*float(lines[8]),7.0*float(lines[7]),6.0*float(lines[6]),5.0*float(lines[5]),4.0*float(lines[4]),3*float(lines[3]),0.0,float(lines[2])]) #first derivative dispersion curve
        disp_curve_diff2=poly1d([56.0*float(lines[8]),42.0*float(lines[7]),30.0*float(lines[6]),20.0*float(lines[5]),12.0*float(lines[4]),6*float(lines[3]),0.0]) #second derivative dispersion curve
        #roton minimum and maximum
        zeroes_diff=disp_curve_diff.r #first derivative roots, 2 and 3 are min and max values respectively
        roton_min=disp_curve(zeroes_diff[1]).real
        roton_max=disp_curve(zeroes_diff[2]).real
        if wavenumber<zeroes_diff[2].real:
            state=int(0) #phonon
        elif wavenumber<zeroes_diff[1].real:
            state=int(1) #R- roton
        else:
            state=int(2) #R+ roton
    
        #finding velocities at given energy
        zeroes=disp_curve_minusomega.r
        real_zeroes=isreal(zeroes)
        final_zeroes=zeros(1)
        for i in xrange(len(zeroes)):
            if real_zeroes[i]==True and real_zeroes[i]>0:
                a_array=array([real(zeroes[i])])
                final_zeroes=sort(append(final_zeroes,a_array))
        velocities=array([disp_curve_diff(final_zeroes[1]),disp_curve_diff(final_zeroes[2]),float(lines[2])])
        final_zeroes=append(final_zeroes,(omega+29.382)/18.2)
        if len(final_zeroes)>3:
            velocities[2]=disp_curve_diff(final_zeroes[3])
        
        #quantum evaporation
        k_a=sqrt(((omega-bindingenergy)*2*mass_helium)/(reduced_plancks_constant_squared))
            
        k=rootsolver(array_maker(omega)) #array of roots for omega^2(k^2)
        
        """
        columns are 
        0: # of ball
        1: x-pos 
        2: x-vel 
        3: y-pos 
        4: y-vel
        5: z-pos
        6: z-vel
        7: state
        8: bounce?
        """
        
        a=zeros((n_particle, 9, len(t))) #3-D array with N rows, 
        for i in xrange(n_particle): #filling initial conditions
            theta_rand=random.rand()*pi*2 #random theta
            phi_rand=random.rand()*pi #random phi
            a[i,0]=i+1 #fill in column 0
            a[i,1,0]=x0 #initial x position
            a[i,1,1]=x0
            a[i,2,1]=vel(state)*sin(phi_rand)*cos(theta_rand)
            a[i,3,0]=y0 #initial y position
            a[i,3,1]=y0
            a[i,4,1]=vel(state)*sin(phi_rand)*sin(theta_rand)
            a[i,5,0]=z0 #initial z position
            a[i,5,1]=z0
            a[i,6,1]=vel(state)*cos(phi_rand)
            a[i,7,0]=state #fill in column 7
            a[i,7,1]=state #fill in column 7
            
        #filling array values
        for j in xrange(2,len(t)): #iterating time values
            for i in xrange(n_particle): #iterating particle values
                vel1=array([a[i,2,j-1],a[i,4,j-1],a[i,6,j-1]])
                pos1=array([a[i,1,j-1],a[i,3,j-1],a[i,5,j-1]])
                vel2=vel1           # new velocity is same as before
                newstate=a[i,7,j-1] # newstate is same as before
                a[i,8,j]=a[i,8,j-1] # bounce status same as before
                if float(a[i,5,j-1])>cyl_height/2.0+vacuum_gap and float(a[i,5,j-2])<cyl_height/2.0+vacuum_gap: #bolometer
                    if newstate==3:
                        bolometer.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    vel2=array([0,0,0])
                if float(a[i,5,j-1])>cyl_height/2.0 and float(a[i,5,j-2])<cyl_height/2.0: #top boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    if newstate!=3:
                        newstate=newstater2(a[i,7,j-1])
                        if newstate==3:
                            vel2=quantumevaporate(vel1,normal)
                            top.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                        else:
                            vel2=bounce(vel1,normal, a[i,7,j-1], newstate)
                            a[i,8,j]=1
                    else:
                        vel2=array([0,0,0])
                if float(a[i,1,j-1])**2.0 + float(a[i,3,j-1])**2.0 > cyl_rad**2.0 and float(a[i,1,j-2])**2.0 + float(a[i,3,j-2])**2.0 < cyl_rad**2.0: #side boundary
                    normal=neg_array(array([a[i,1,j-1],a[i,3,j-1],-0]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    else:
                        vel2=array([0,0,0])
                        atom_side.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                if float(a[i,5,j-1])<-1.0*cyl_height/2.0: #bottom boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                
                a[i,2,j]=vel2[0]
                a[i,4,j]=vel2[1]
                a[i,6,j]=vel2[2]
                a[i,7,j]=newstate
                a[i,1,j]=x_pos(i,j-1)+del_t*x_vel(i,j)
                a[i,3,j]=y_pos(i,j-1)+del_t*y_vel(i,j)
                a[i,5,j]=z_pos(i,j-1)+del_t*z_vel(i,j)
                

        #initial conditions
        y0 = float(0)
        
        a=zeros((n_particle, 9, len(t))) #3-D array with N rows, 
        for i in xrange(n_particle): #filling initial conditions
            theta_rand=random.rand()*pi*2 #random theta
            phi_rand=random.rand()*pi #random phi
            a[i,0]=i+1 #fill in column 0
            a[i,1,0]=x0 #initial x position
            a[i,1,1]=x0
            a[i,2,1]=vel(state)*sin(phi_rand)*cos(theta_rand)
            a[i,3,0]=y0 #initial y position
            a[i,3,1]=y0
            a[i,4,1]=vel(state)*sin(phi_rand)*sin(theta_rand)
            a[i,5,0]=z0 #initial z position
            a[i,5,1]=z0
            a[i,6,1]=vel(state)*cos(phi_rand)
            a[i,7,0]=state #fill in column 7
            a[i,7,1]=state #fill in column 7
            
        #filling array values
        for j in xrange(2,len(t)): #iterating time values
            for i in xrange(n_particle): #iterating particle values
                vel1=array([a[i,2,j-1],a[i,4,j-1],a[i,6,j-1]])
                pos1=array([a[i,1,j-1],a[i,3,j-1],a[i,5,j-1]])
                vel2=vel1           # new velocity is same as before
                newstate=a[i,7,j-1] # newstate is same as before
                a[i,8,j]=a[i,8,j-1] # bounce status same as before
                if float(a[i,5,j-1])>cyl_height/2.0+vacuum_gap and float(a[i,5,j-2])<cyl_height/2.0+vacuum_gap: #bolometer
                    if newstate==3:
                        bolometer.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    vel2=array([0,0,0])
                if float(a[i,5,j-1])>cyl_height/2.0 and float(a[i,5,j-2])<cyl_height/2.0: #top boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    if newstate!=3:
                        newstate=newstater2(a[i,7,j-1])
                        if newstate==3:
                            vel2=quantumevaporate(vel1,normal)
                            top.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                        else:
                            vel2=bounce(vel1,normal, a[i,7,j-1], newstate)
                            a[i,8,j]=1
                    else:
                        vel2=array([0,0,0])
                if float(a[i,1,j-1])**2.0 + float(a[i,3,j-1])**2.0 > cyl_rad**2.0 and float(a[i,1,j-2])**2.0 + float(a[i,3,j-2])**2.0 < cyl_rad**2.0: #side boundary
                    normal=neg_array(array([a[i,1,j-1],a[i,3,j-1],-0]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    else:
                        vel2=array([0,0,0])
                        atom_side.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                if float(a[i,5,j-1])<-1.0*cyl_height/2.0: #bottom boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                
                a[i,2,j]=vel2[0]
                a[i,4,j]=vel2[1]
                a[i,6,j]=vel2[2]
                a[i,7,j]=newstate
                a[i,1,j]=x_pos(i,j-1)+del_t*x_vel(i,j)
                a[i,3,j]=y_pos(i,j-1)+del_t*y_vel(i,j)
                a[i,5,j]=z_pos(i,j-1)+del_t*z_vel(i,j)
                
    
        #initial conditions
        y0 = float(-.025)
        
        a=zeros((n_particle, 9, len(t))) #3-D array with N rows, 
        for i in xrange(n_particle): #filling initial conditions
            theta_rand=random.rand()*pi*2 #random theta
            phi_rand=random.rand()*pi #random phi
            a[i,0]=i+1 #fill in column 0
            a[i,1,0]=x0 #initial x position
            a[i,1,1]=x0
            a[i,2,1]=vel(state)*sin(phi_rand)*cos(theta_rand)
            a[i,3,0]=y0 #initial y position
            a[i,3,1]=y0
            a[i,4,1]=vel(state)*sin(phi_rand)*sin(theta_rand)
            a[i,5,0]=z0 #initial z position
            a[i,5,1]=z0
            a[i,6,1]=vel(state)*cos(phi_rand)
            a[i,7,0]=state #fill in column 7
            a[i,7,1]=state #fill in column 7
            
        #filling array values
        for j in xrange(2,len(t)): #iterating time values
            for i in xrange(n_particle): #iterating particle values
                vel1=array([a[i,2,j-1],a[i,4,j-1],a[i,6,j-1]])
                pos1=array([a[i,1,j-1],a[i,3,j-1],a[i,5,j-1]])
                vel2=vel1           # new velocity is same as before
                newstate=a[i,7,j-1] # newstate is same as before
                a[i,8,j]=a[i,8,j-1] # bounce status same as before
                if float(a[i,5,j-1])>cyl_height/2.0+vacuum_gap and float(a[i,5,j-2])<cyl_height/2.0+vacuum_gap: #bolometer
                    if newstate==3:
                        bolometer.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    vel2=array([0,0,0])
                if float(a[i,5,j-1])>cyl_height/2.0 and float(a[i,5,j-2])<cyl_height/2.0: #top boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    if newstate!=3:
                        newstate=newstater2(a[i,7,j-1])
                        if newstate==3:
                            vel2=quantumevaporate(vel1,normal)
                            top.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                        else:
                            vel2=bounce(vel1,normal, a[i,7,j-1], newstate)
                            a[i,8,j]=1
                    else:
                        vel2=array([0,0,0])
                if float(a[i,1,j-1])**2.0 + float(a[i,3,j-1])**2.0 > cyl_rad**2.0 and float(a[i,1,j-2])**2.0 + float(a[i,3,j-2])**2.0 < cyl_rad**2.0: #side boundary
                    normal=neg_array(array([a[i,1,j-1],a[i,3,j-1],-0]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    else:
                        vel2=array([0,0,0])
                        atom_side.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                if float(a[i,5,j-1])<-1.0*cyl_height/2.0: #bottom boundary
                    normal=neg_array(array([0,0,a[i,5,j-1]]))
                    theta_i=angle(vel1,neg_array(normal))
                    if newstate != 3:
                        newstate=newstater(a[i,7,j-1])
                        if newstate==3:
                            vel2=array([0,0,0])
                            transmission.append([wavenumber,y0,t[j],a[i,1,j-1],a[i,3,j-1],a[i,5,j-1],a[i,2,j-1],a[i,4,j-1],a[i,6,j-1],a[i,7,j-1],a[i,8,j-1]])
                    if newstate != 3:
                        vel2=bounce(vel1,normal,a[i,7,j-1],newstate)
                        a[i,8,j]=1
                
                a[i,2,j]=vel2[0]
                a[i,4,j]=vel2[1]
                a[i,6,j]=vel2[2]
                a[i,7,j]=newstate
                a[i,1,j]=x_pos(i,j-1)+del_t*x_vel(i,j)
                a[i,3,j]=y_pos(i,j-1)+del_t*y_vel(i,j)
                a[i,5,j]=z_pos(i,j-1)+del_t*z_vel(i,j)
        
#bolometer_save="bolometer_"+str(randloopnumber)+".csv"
#transmission_save="transmission_"+str(randloopnumber)+".csv"
#top_save="top_"+str(randloopnumber)+".csv"
#atom_side_save="atomside_"+str(randloopnumber)+".csv"
bolometer_save="/global/homes/e/epease/bolometer2/bolometer_"+str(randloopnumber)+".csv"
transmission_save="/global/homes/e/epease/transmission2/transmission_"+str(randloopnumber)+".csv"
top_save="/global/homes/e/epease/top2/top_"+str(randloopnumber)+".csv"    
atom_side_save="/global/homes/e/epease/atomside2/atomside_"+str(randloopnumber)+".csv"       
savetxt(transmission_save, transmission, delimiter=",")
savetxt(top_save, top, delimiter=",")
savetxt(bolometer_save, bolometer, delimiter=",")
savetxt(atom_side_save, atom_side, delimiter=",")