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


# For reproducibility, set a seed for randomly generated inputs. Change to your favourite integer.

# Set simulation parameters
Nd = 3 #number of spatial dimensions, minimum 2
Np = 51 #number of particles
M = 0.000001 #mass of a single particle
l = 1.0 #half of the length of the box


# Set animation parameters
Nt = 900 #total number of timesteps for animation
frame_duration = 30 #number of milliseconds displayed for each timestep

# Choose projection for the main panel
project_3d = True

# Generate initial conditions
np.random.seed(4080) #generate the random seed. Use set value for reproducibility.
#position = 1-2*np.random.random((Nd,Np)) ##generate position at random within box
position = np.zeros((Nd, Np)) #initial position of all particle at (0,0,0)
position[0][0] += 0.4 #initial position of the orbiting particle


v_max= 0.01 #maximum drift velocity of particle
velocity = np.zeros((Nd, Np))  #initial velocity of all particles is zero
#+ v_max*(1-2*np.random.random((Nd,Np))) 
    #generate initial velocity of each particle at random
velocity[1][0] += 0.01 #initial velocity of the orbiting particle
accel = np.zeros((Nd, Np))



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
fig = figure(figsize=(12,6)) # Create frame and set size
subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.15,hspace=0.2)
# Create one set of axes as the left hand panel in a 1x2 grid
if project_3d:
    ax1 = subplot(121,projection='3d') # For very basic 3D projection
    ax1.set_zlim(-1,1)  # Set z-axis limits
else:
    ax1 = subplot(121) # For normal 2D projection
xlim(-1,1)  # Set x-axis limits
ylim(-1,1)  # Set y-axis limits


# Create command which will plot the positions of the particles
if project_3d:
    points, = ax1.plot([],[],[],'o',markersize=4)  ## For 3D projection
else:
    points, = ax1.plot([],[],'o',markersize=4) ## For 2D projection

ax2 = subplot(222) # Create second set of axes as the top right panel in a 2x2 grid
xmax = 4 # Set xaxis limit
xlim(0,xmax) # Apply limit
xlabel('Separation')
ylabel('No. of particle pairs')
dx=0.2 # Set width of x-axis bins
ylim(0,dx*Np*10) # Reasonable guess for suitable yaxis scale    
xb = np.arange(0,xmax+dx,dx)  # Set x-axis bin edges
xb[0] = 1e-6 # Shift first bin edge by a fraction to avoid showing all the zeros (a cheat, but saves so much time!)
line, = ax2.plot([],[],drawstyle='steps-post') # Define a command that plots a line in this panel

ax4 = plt.subplot(224) # Create last set of axes as the bottom right panel in a 2x2 grid

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
    h,x = histogram(np.ravel(tril(separation(position))),bins=xb)
        # Make histogram of the lower triangle of the seperation matrix. Density = True normalize it.
    line.set_data(x[:-1],h) # Set the new data for the line in the 2nd panel      
    return points,line, # Plot the points and the line
    
# Create animation
# https://matplotlib.org/api/_as_gen/matplotlib.animation.FuncAnimation.html
ani = animation.FuncAnimation(fig, update, frames=Nt,interval = frame_duration)
#plt.show()

ani.save("orbits.mp4")