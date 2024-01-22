import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob


def GenerateInitialConfig(no_particles, density, ndim, vSTDEV, Temp):
    #Initialise position and velocity matrices
    coords = np.zeros((no_particles, ndim))
    vels = np.zeros((no_particles, ndim))

    #XYZ folder
    os.mkdir('XYZFiles')
    os.chdir('XYZFiles')

    

    #Get box side length
    len = (no_particles/density)**(1/ndim)

    #Get smallest perfect cube/square/^1 number greater or equal to number of particles
    ncube = 2
    while ncube**ndim < no_particles:
        ncube = ncube + 1

    #Position particles in the box
    index = np.zeros((ndim))
    adder = np.resize(np.tile(np.array([0.5]), ndim), (ndim))
    multiplier = len/ncube

    for i in range(no_particles):
        coords[i] = (index + adder) * multiplier
        vels[i] = np.random.normal(0, vSTDEV, (1, ndim))

        #Advance index
        index[0] = index[0] + 1
        if index[0] == ncube:
            index[0] = 0
            index[1] = index[1] + 1
            if index[1] == ncube:
                index[1] = 0
                index[2] = index[2] + 1
    #Save initial config to XYZ folder
    WriteXYZFile(coords, no_particles, 0)


    #Set initial momentum to 0
    totV = np.resize(np.sum(vels) / no_particles, (1, ndim))
    vels = vels - totV

    #Set initial KE
    msv = np.sum(np.multiply(vels, vels)) / no_particles  #Mean squared velocity
    velscale = np.sqrt(ndim * Temp / msv)
    vels = vels * velscale

    #Get initial forces
    forces = CalcForces(coords, len, no_particles, ndim)

    return vels, coords, len, forces
    
    
def CalcForces(coords, len, no_particles, ndim):
    #Initialise forces
    forces = np.zeros(coords.shape)

    #Loop over all particle pairs
    for particleA in range(no_particles - 1):
        for particleB in range(particleA + 1, no_particles):
            dr = coords[particleA] - coords[particleB]
            dr = distPBC3D(dr, len, ndim)
            dr2 = np.dot(dr, dr)

            #LJ potential
            forceFact = 1/(dr2)**4 * (1/(dr2)**3-0.5)
            forces[particleA] = forces[particleA] + dr*forceFact
            forces[particleB] = forces[particleB] - dr*forceFact
    
    forces = forces * 48
    
    return forces


def distPBC3D(vec, len, ndim):
    hL = len / 2
    for i in range(ndim):
        if vec[i] > hL:
            vec[i] = vec[i] - len
        elif vec[i] < -hL:
            vec[i] = vec[i] + len
    return vec

def WriteXYZFile(x_coords, no_particles, iteration):
    xyzfile = open('config' + str(iteration) + '.xyz', 'w')
    #xyzfile.write(str(no_particles) + '\n')
    #xyzfile.write("\n")
    for i in range(no_particles):
        xyzfile.write(str(x_coords[i,0]) + " " + str(x_coords[i,1]) + " " + str(x_coords[i,2]) + '\n')
    xyzfile.close()

def RunSimulation(no_particles, density, ndim, vSTDEV, temp, no_steps, dt):
    m = 40

    #Initialise velocities and positions
    vels, coords, len, forces = GenerateInitialConfig(no_particles, density, ndim, vSTDEV, temp)
    t = 0

    #Run loop
    for i in range(no_steps):
        next_coords = PositionAlgorithm(coords, vels, dt, forces, m)
        '''for j in range(no_particles):
            next_coords[j] = distPBC3D(next_coords[j], len, ndim)'''

        next_forces = CalcForces(next_coords, len, no_particles, ndim)
        next_vels = VelocityAlgorithm(vels,forces, next_forces, dt, m)
        t = t + dt

        

        vels = next_vels
        coords = next_coords
        forces = next_forces


        #Take measurements
        WriteXYZFile(coords, no_particles, i+1)
        #Put functions for taking measurments/saving coordinates here


    
    print('Simulation Complete')
        
    


def PositionAlgorithm(coords, vels, dt, forces, m):
    return coords + vels*dt + 1/(2*m) * forces * dt**2

def VelocityAlgorithm(vels, forces, next_forces, dt, m):
    return vels + dt/(2*m) * (forces + next_forces)





RunSimulation(27, 3, 3, 5, 298, 1000, 0.0001)

