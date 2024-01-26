import numpy as np
import matplotlib.pyplot as plt
import datetime as datetime

#Programme to calculate how pressure and enthalpy evolves over 50 configurations

#Constants
no_configs = 50
m = 1E-27
e = 1.67E-21 #In J
lama = 49
lamr = 50
LJsigma = 3.405
PHSsigma = 3.305

#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
        self.velocities = None
        self.V = None
        self.a = None

#Create confiiguration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    
    
    configurations = []
    
    #Generate configuration class
    for i in range(no_configurations):
        config = configuration()
        config.no_particles = int(lines[3])
        config.timestep = int(lines[(config.no_particles + 9) *i + 1])
        config.coords = np.zeros((config.no_particles, 3))
        config.velocities = np.zeros((config.no_particles, 3))

        #Get volume of box
        a = np.zeros((3))
        for axis in range(3):
            x = lines[5 + (config.no_particles + 9) * i + axis].split()
            a[axis] = float(x[1]) - float(x[0])
        vol = np.prod(a)

        config.a = a
        config.V = vol

        for j in range(config.no_particles):
            xyzline = lines[i * (config.no_particles + 9) + 9 + j].split()
            config.coords[j] = xyzline[4:7]
            config.velocities[j] = xyzline[7:10]

        #Scale positions and vels into Angstrom
        config.coords *= a
        config.velocities *= a
        
        configurations.append(config)

    return configurations


def GetEnthalpy(config, type):
    #(N,N,3) tensor containing all interparticle distances in the x, y, z directions
    rij = config.coords[:, np.newaxis, :] - config.coords
    #(N,N) matrix containing the norm of all interparticle distances
    r = np.linalg.norm(rij, axis = 2)
    forces, energy = GetForceEnergy(r, rij, type)
    #Volume in m^3
    vol_SI = config.V * 10E-30

    #To calculate the total work the axis matrices of the rij and forces tensors must be multplied element wise and their sum taken, gives work in J
    force_term = np.sum(rij * forces)
    
    #Want mv^2 => Multiply velocities element wise and take sum, 1E+10 conversion factor to convert to J
    kin_term = np.sum(config.velocities * config.velocities) * m * 1E+10

    pressure = (force_term + kin_term) / (3 * vol_SI)
    enthalpy = pressure * vol_SI * + energy


    return enthalpy, pressure, energy #Now in units J, Pa, J

#Get all forces and energies of the system, this is done by creating a creating an (N,N) matrix containing all the forces according to the potential and the 
#norm of particle distances, to get the forces for every particle along every direction this force matrix is broadcasted into a (N,N,3) tensor where the
#(i, j, k) element is given by F(i,j) * rij(i, j, k) / r(i, j) where i and j indicate the particle pair and k represents the axis, so essentially have 3 matrices each containing all the forces along a specified axis, these matrices
#are then stacked, as the energies are only dependent on the norm the sum of all the elements in the energy matrix can be taken
def GetForceEnergy(r, rij, type):
    
    if type == 'LJ':
        np.fill_diagonal(r, 3 * LJsigma) #To prevent dividing by 0 errors the r matrix is filled with the cutoff value along i=j
        f = GetLJForce(r, LJsigma)[:, :, np.newaxis] * rij / r[:, :, np.newaxis]
        energy = np.sum(GetLJEnergy(r, LJsigma))
    
    elif type == 'PHS':
        np.fill_diagonal(r, 1.5 * PHSsigma)
        f = GetPHSForce(r, PHSsigma)[:, :, np.newaxis] * rij / r[:, :, np.newaxis]
        energy = np.sum(GetPHSEnergy(r, PHSsigma))
    return f, energy

#Functions for getting all forces and energies according to LJ potential    
def GetLJEnergy(r, sigma):
    return np.where(r < 3 * sigma, 4 * e * ((sigma/r)**12 -  (sigma/r)*6), 0)
def GetLJForce(r, sigma):
    return np.where(r < 3 * sigma , 24 * e * 1/r * sigma * (2*(sigma/r)*12 - (sigma/r)**6), 0)

#Functions for getting all forces and energies according to PHS potential
def GetPHSEnergy(r, sigma):
    return np.where(r < 1.5 * sigma, lamr * (lamr / lama) ** lama * e * ((sigma / r) ** (lamr) - (sigma / r)**(lama) ) + e, 0)
def GetPHSForce(r, sigma):
    b = - lamr * ((lamr / lama) ** lama) * e
    c1 = lama * (sigma ** lama) * ((1 / r) ** (lama + 1))
    c2 = lamr * (sigma ** lamr) * ((1 / r) ** (lamr + 1))

    return np.where(r < 1.5 * sigma, b * (c1 - c2), 0)


def InstPressure(filename, no_configurations, type):
    configs = InitConfigs(filename, no_configurations)
    plist = []
    elist= []
    hlist = []
    time = []

    for i in range(len(configs)):
        h, p, e = GetEnthalpy(configs[i], type)
        plist.append(p)
        hlist.append(h)
        elist.append(e)
        time.append(configs[i].timestep)
        print('Config ' + str(i+1) + ' of ' + str(len(configs)) + ' Complete' )
    
    plt.plot(time, elist)
    plt.show()
    plt.plot(time, hlist)
    plt.show()
    plt.plot(time, plist)
    plt.show()


    return







InstPressure('pres.xyz', no_configs, 'LJ')