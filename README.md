# N-body-Simulation
N-body gravitational simulation
This is the code for PHYS7280 Project 1, by Fangjun ZHAO 44362805

The code runs in Python, and it is a N-body simulation that calculates the gravitational interactions between N identical particles in a periodic space.

Periodic boundary condition means when the particle crosses a boundary, it will automatically appear on the other side.

The simulation outputs an animation of the particle's evolution, as well as the separation correlation function and the corresponding power spectrum.

You can freely change the parameters of the simulation located at the start of the file.


There are two codes incuded:
1. Project_1_Orbital.py is the code that shows the gravitational interaction leads to orbital motion. The default parameter shows a single particle 
with pre-defined velocity been gravitationally attracted to cluster of 50 particles. The single particle enters an orbit. 

2. Project_1_random.py is the code that simulates N particles uniformly distributed in a periodic space with some random drift velocity. The default
paramter shows 100 particles interacting in a 10x10x10 3D space.

Below is the list of parameters that can be changed, and their effect:
Nd - changes the dimension particle lives in
Np - number of particles
M - mass of the particle. Primarily affect strength of gravity.
l - Half of the length of the box. Determines the size of the periodic space.
Nt - total number of time step for animation. Affects the length of the animation.
frame_duration - in ms, affects the speed of animation.
v_max - magnitude of the drift velocity.
klim - size of the k-space of power spectrum we're interested in.
