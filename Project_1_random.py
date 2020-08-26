#PHYS7280 Project 1. Simulation
#By Fangjun Zhao
#Code based on files from Prof. Tamara Davis
#Code from file drifers.py and panels.py
#This code is for the planet-orbit simulation. It makes a cluster of
#50 particles, and a single particle orbiting around it.

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Set simulation parameters
Nd = 3 #number of spatial dimensions, minimum 2
Np = 100 #number of particles
M = 0.00001 #mass of a single particle
l = 5.0 #half of the length of the box


# Set animation parameters
Nt = 1200 #total number of timesteps for animation
frame_duration = 30 #number of milliseconds displayed for each timestep

# Choose projection for the main panel
project_3d = True

# Generate initial conditions
np.random.seed(4080) #generate the random seed. Use set value for reproducibility.
position = l-2*l*np.random.random((Nd,Np))

v_max= 0.1 #maximum drift velocity of particle
velocity = np.zeros((Nd, Np)) + v_max*(1-2*np.random.random((Nd,Np))) 
    #generate initial velocity of each particle at random
accel = np.zeros((Nd, Np))

# Define Parameter for Power Spectrum
klim = 100 #define size of k space we're interested in
x2 = np.array(range(1,klim+1)) #generate a list of k


# Function to apply boundary conditions
def apply_boundary(p):
    for i in range(Nd):    
        p[i] = [x - 2*l*((x+l)//(2*l)) for x in p[i]] #periodic boundary
    return p

# Function to find separations from position vectors
def separation(p): 
    s = p[:,None,:] - p[:,:,None] # find N x N x Nd matrix of particle separations
    return np.sum(s**2,axis=0)**0.5 # return N x N matrix of scalar separations

sft = 0.1 #softening parameter of gravity

# Function to find acceleration of particle at each timestep
def acceleration (p):
    F = np.zeros((Nd, Np))
    sep = p[:,None,:] - p[:,:,None] # find N x N x Nd matrix of particle separations
    s = (np.sum(sep**2,axis=0)+sft)**0.5 # return N x N matrix of scalar separations
    for k in range (Nd):    #for each coordinate
        for i in range (Np): #for all particles
            totF = 0
            for j in range (Np): #for each particle pair
                    totF += sign(p[k][j] - p[k][i])/s[i][j] #calculate Force and add it to the sum
            F[k][i] = totF #sum of the force
    result = M*F #apply mass to get acceleration from force
    return result



# Set the axes on which the points will be shown
plt.ion() # Set interactive mode on
fig = figure(figsize=(16,8)) # Create frame and set size
subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.15,hspace=0.2)
# Create one set of axes as the left hand panel in a 1x2 grid
if project_3d:
    ax1 = subplot(121,projection='3d') # For very basic 3D projection
    ax1.set_zlim(-l,l)  # Set z-axis limits
else:
    ax1 = subplot(121) # For normal 2D projection
xlim(-l,l)  # Set x-axis limits
ylim(-l,l)  # Set y-axis limits


# Create command which will plot the positions of the particles
if project_3d:
    points, = ax1.plot([],[],[],'o',markersize=4)  ## For 3D projection
else:
    points, = ax1.plot([],[],'o',markersize=4) ## For 2D projection

ax2 = subplot(222) # Create second set of axes as the top right panel in a 2x2 grid
ax2.grid(True)
xmax = floor((12*l*l)**(1/2))+2 # Set xaxis limit of correlation function
xlim(0,xmax) # Apply limit
xlabel('Separation')
ylabel('Correlation Function')
dx=0.2 # Set width of x-axis bins
ylim(-1.5, 5) #Normalized correlation function yaxis scale
xb = np.arange(0,xmax+dx,dx)  # Set x-axis bin edges
xb[0] = 1e-6 # Shift first bin edge by a fraction to avoid showing all the zeros (a cheat, but saves so much time!)
fine, = ax2.plot([],[],drawstyle='steps-post') # Define a command that plots a line in this panel

ax4 = plt.subplot(224) # Create last set of axes as the bottom right panel in a 2x2 grid
xmax = 4 # Set xaxis limit
xlim(1,klim) # Apply limit
xlabel('k')
ylabel('Power Spectrum')
ylim(-1,10) # Reasonable guess for suitable yaxis scale    
ax4.xb = np.arange(0,xmax+dx,dx)  # Set x-axis bin edges
ax4.xb[0] = 1e-6 # Shift first bin edge by a fraction to avoid showing all the zeros (a cheat, but saves so much time!)
line, = ax4.plot([],[],drawstyle='steps-post') # Define a command that plots a line in this panel

# Find correlation function
def RR(n): 
    #find the distribution of random separation. Generate n position distributions,
    #then find the separation histogram for each of them,
    #then find the average for each bin
    tot = 0 
    for i in range(n):
         h =  histogram(np.ravel(tril((separation(l-2*l*np.random.random((Nd,Np)))))),bins=xb)[0]
         tot += h
    RR = tot/n
    return RR 

# Calculate Power Spectrum by integration
def powerSpec(klim, h, x):
    p = []
    for k in range (1,klim+1):
        tot = 0
        for i in range(h.size):
            tot += dx*x[i]*h[i]*sin(k*x[i])/k #rectangular integration
        p = np.append(p, tot)
    return 2*pi*p

RRD = RR(100) #random distribution of separation

# Define procedure to update positions at each timestep
def update(i):
    global position,velocity # Get positions and velocities
    accel = acceleration(position) #update accelerations
    velocity += accel #increase velocity according to their accelerations
    position += velocity # Increment positions according to their velocites
    position = apply_boundary(position) # Apply boundary conditions
    points.set_data(position[0,:], position[1,:]) # Show 2D projection of first 2 position coordinates
    if project_3d:
        points.set_3d_properties(position[2,:])  ## For 3D projection
    r,x1 = histogram(np.ravel(tril((separation(position)))),bins=xb) 
        #generate histogram array of particle-particle separation
    h = r.astype(double) #change type to double for the continuous correlation function
    for i in range(h.size):
        if RRD[i] > 1:
            h[i] = (double(h[i])/RRD[i])-1 #Peebles-Hauser estimator
        else:
            h[i] = -1 #we ignore bins without enough data
        # Make histogram of the lower triangle of the seperation matrix. Density = True normalize it.
    fine.set_data(x1[:-1],h) # Set the new data for the fine in the 2nd panel
    p = powerSpec(klim, h, x1)
    line.set_data(x2[1:], p[1:]) # Set the new data for the line in the 3rd panel       
    return points,line,fine, # Plot the points and the line
    
# Create animation
# https://matplotlib.org/api/_as_gen/matplotlib.animation.FuncAnimation.html
ani = animation.FuncAnimation(fig, update, frames=Nt,interval = frame_duration)
#plt.show()

ani.save("random_default.mp4")