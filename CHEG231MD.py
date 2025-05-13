"""
CHEG231MD.py
Simple MD code in Python (it's slow!) for a Lenard-Jones fluid
based on QuickBASIC code by Richard L. Rowley [1]
and FORTRAN code by Smit and Frenkel [2]

Eric M. Furst
November 2023

[1] Richard L. Rowley (1994). Statistical Mechanics for Thermophysical Property Calculations. Prentice Hall, New York.
[2] Daan Frenkel and Berend Smit (2002). Understanding Molecular Simulation, 2nd ed. Academic Press, New York.

Revisions
October 2024 - sped up with numba jit of accel routine and LJ calc.

Requires:
    python (tested with ver. 3.12.4 build h99e199e_1 channel defaults)
    numba  (tested with ver. 0.60.0 build py312hd77ebd4_0 channel defaults)
    numpy  (tested with ver. 1.26.4 build py312h7f4fdc5_0 channel defaults)

Needs:
    nn - total number of particles per dimension (thus, N = nn^3)
    rho - dimensionless number density of particles (units of rho/sigma^3)
    vmax - maximum velocity (units of sqrt[epsilon/m])

Next to do:
    write a 2D animated version
    a challenge for students: dimensionless variables
    can a fast r_cutoff be implemented based on 
    from scipy.spatial.distance import pdist, squareform?
"""

import numpy as np
from numba import jit

X, Y, Z = 0, 1, 2 # helpful for indexing

class MDSimulation:
    """
    x - new position coordinate array (N x 3)
    v - velocity array (N x 3)
    """

    def __init__(self, nn, rho, vmax):
        """
        Initialize the simulation
        - particles start on a cubic lattice
        - velocities are assigned a random fraction of the maximum
        - velocities are corrected to eliminate net flow
        """
    
        # initial variables
        self.N = nn**3  # total number of particles
        self.rho = rho
        self.vol = self.N/rho # volume of the simulation 
        self.L = self.vol**(1/3) # size of the simulation
        self.dL = self.L/nn # initial cubic array lattice size

        self.dt = 0.005 # dimensionless time step (0.005 is common)
        self.dt2 = self.dt*self.dt
        
        # initialize position and acceleration arrays 
        self.x = np.zeros((self.N,3)) # init N by 3 array of positions
        self.xnew = np.zeros((self.N,3)) # init N by 3 array of positions
        self.a = np.zeros((self.N,3)) # init N by 3 array of accelerations

        # state of the simulation
        self.pe = 0 # total potential energy
        self.ke = 0 # total kinetic energy
        self.vir = 0 # virial coef
        self.T = 1 # temperature
        self.P = 0 # pressure

        # create a cubic array of particles
        particle = 0
        for i in range(0, nn):
            for j in range(0,nn):
                for k in range(0,nn):
                    self.x[particle,X]=(i+0.5)*self.dL
                    self.x[particle,Y]=(j+0.5)*self.dL
                    self.x[particle,Z]=(k+0.5)*self.dL
                    particle += 1
        
        # assign random velocities from uniform distribution (units?)
        self.v = vmax*np.random.rand(self.N,3) # init N x 3 array of velocities

        # normlize the velocities s/t net velocity is zero
        vcm = np.sum(self.v,axis=0)/self.N # 1x3 array of <vx>, <vy>, <vz>    
        for particle in range(0,self.N):
            self.v[particle,:] -= vcm
 
        # set the initial accelerations
        self.accel()

    def move(self):
        """
        Verlet velocity algorithm
        """
        # find new positions
        self.xnew = self.x + self.v*self.dt + self.a*self.dt2/2.

        # apply periodic boundary conditions
        beyond_L = self.xnew[:,:] > self.L # N x 3 Boolean array
        self.xnew[beyond_L] -= self.L
        beyond_zero = self.xnew[:,:] < 0
        self.xnew[beyond_zero] += self.L

        self.x = self.xnew   # update positions
        self.v = self.v + self.a*self.dt/2 #half update velocity

        # call accelerate 
        self.accel()

        # finish velocity update
        self.v = self.v + self.a*self.dt/2 

        # update properties
        self.ke = 0.5*np.sum(self.v*self.v)
        self.T = 2*self.ke/self.N/3
        self.P = self.T*self.rho + self.rho/self.N*self.vir/3


    # Non-JIT-compiled class method
    def accel(self):
        """                                                                     
        Acceleration calculated using F = ma                                    
        """
        self.a = np.zeros((self.N, 3))  # Zero out all accelerations
        self.pe = 0
        self.vir = 0

        # Call the JIT-compiled function to calculate acceleration
        self.pe, self.vir = accel_jit(self.N, self.x, self.a, self.L)

# Functions are below
# We're using numba jit to make the loop of the calculation faster

# JIT-compiled function to calculate acceleration
# Note: numba works by "pass in reference" with numpy arrays,
# so self.a values are modified and used in self.move()
# in addition to returning the values of potential energy and
# the virial factor
@jit(nopython=True, parallel=False)
def accel_jit(N, x, a, L):
    pe = 0
    vir = 0

    for i in range(0, N - 1):
        for j in range(i + 1, N):
            r2 = 0
            dx = x[i, :] - x[j, :]  # 1x3 array
            # Apply minimum image convention over periodic boundary conditions
            for k in range(3):  # Iterate over X, Y, Z
                if abs(dx[k]) > 0.5 * L:
                    dx[k] -= np.sign(dx[k]) * L
                r2 += dx[k] * dx[k]

            # Calculate the force, potential, and virial contribution
            forceri, potential, virij = LJ(r2)
            a[i, :] += forceri * dx
            a[j, :] -= forceri * dx
            pe += potential
            vir += virij

    return pe, vir

# JIT-compiled function to calculate interaction
@jit(nopython=True)
def LJ(r2):
    """
    Compute force and potential from LJ interaction 
    r2 is r^2, forceri returns force/r
    """
    r2i = 1/r2
    r6i = r2i**3
    potential = 4*(r6i*r6i-r6i)
    virij = 48*(r6i*r6i-0.5*r6i)
    forceri = virij*r2i
    return forceri, potential, virij