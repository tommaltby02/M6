import numpy as np
import matplotlib.pyplot as plt

#Programme to calculate energy of conf.xyz using different models

#Constants
kb = 1.38E-23 
temp = 179.81
mass = 6.63E-26
sigma1 = 3.405
sigma2 = 3.05
er = 119.87
lama = 49
lamr = 50
k = 5
no_configs = 1


#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
        self.particle_types = None
        self.velocities = None
        self.V = None
        self.a = None
        self.type_mat = None

#Create confiiguration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    configurations = []

    for i in range(no_configurations):
        config = configuration()
        config.no_particles = int(lines[3])
        config.timestep = int(lines[i * (config.no_particles + 9) + 1])
        config.coords = np.zeros((config.no_particles, 3))
        config.particle_types = np.zeros((config.no_particles))

        #Get volume of box
        a = np.zeros((3))
        for axis in range(3):
            x = lines[5 + (config.no_particles + 9) * i + axis].split()
            a[axis] = float(x[1]) - float(x[0])
        V = np.prod(a)

        config.V = V
        config.a = a


        for j in range(config.no_particles):
            config.coords[j] = lines[i * (config.no_particles + 9) + 9 + j].split()[4:7]
            config.particle_types[j] = lines[i * (config.no_particles + 9) + 9 + j].split()[2]
        
        #Scale coordinates into Angstrom
        config.coords *= a
        
        #Generate velocities
        config = GenerateVelocities(config)
        #Get type matrix
        config = GetTypeMatrix(config)

        #Get histogram
        #GenerateVelHistogram(config)


        configurations.append(config)

    return configurations

#Use normal distribution with STD derived from Maxwell-Boltzmann distribution
def GenerateVelocities(config):
    std = np.sqrt(kb * temp / mass)
    config.velocities = np.random.normal(0, std, size = (config.no_particles, 3))
    return config

#Generate matrix of overlaps of particle types, entry will be 1 for particles of the same type and -1 for particles of differing types, element multiplying this 
#by the matrix generated according to the YDH potential will give the energy for each particle pair
def GetTypeMatrix(config):
    vec = np.reshape(np.where(config.particle_types == 2, -1, 1),(config.no_particles, 1))
    config.type_mat = np.matmul(vec, vec.transpose())
    return config


#Generate histogram of velocities
def GenerateVelHistogram(config):
    for i in range(3):
        plt.hist(config.velocities[:, i], bins = 50, density = True, alpha = 0.7, label = f'Axis {i+1}')
        plt.xlabel('Velocity')
        plt.ylabel('Prob Density')
        plt.legend()
        plt.show()

#Fill diagonals to prevent dividing by 0 errors also divide by 2 to prevent double counting of particle pairs
def Potential(config, type):
    r = np.linalg.norm(config.coords[:, np.newaxis, :] - config.coords, axis = 2)
    if type == 'LJ':
        np.fill_diagonal(r, 3 * sigma1)
        pot = GetLJEnergy(r) / 2
    elif type == 'PHS':
        np.fill_diagonal(r, 1.5 * sigma1)
        pot = GetPHSEnergy(r) / 2
    elif type == 'HSYDH':
        np.fill_diagonal(r, 3.5 * sigma2)
        pot =  (GetHSEnergy(r) + GetYDHEnergy(r, config)) / 2
    return pot

def Kinetic(config):
    v = np.linalg.norm(config.velocities, axis = 1)
    return 0.5 * mass * np.sum(v**2) 

#Calculate pair potential using Lennard-Jones model
def GetLJEnergy(r):
    return np.sum(np.where(r < 3 * sigma1, 4 * er * kb * ((sigma1/r)**12 - (sigma1/r)**6), 0))

#Calculate pair potential using Pseudo Hard Sphere model
def GetPHSEnergy(r):
    return np.sum(np.where(r < 1.5 * sigma1, lamr * (lamr / lama) ** lama * er * kb * ((sigma1/r)**lamr - (sigma1/r)**lama) + er * kb, 0))
   
#Calculate pair potential using Hard Sphere model
def GetHSEnergy(r):
    return np.sum(np.where(r < sigma2, float('inf'), 0))

#Calculate pair potential using Yukawa Debye Huckel model
def GetYDHEnergy(r, config):
    return np.sum(config.type_mat * np.where(r < 3.5 * sigma2, er * kb * sigma2 / r * np.exp(-k * (r - sigma2)), 0))
    
def CheckTemp(ke, config):
    avgke = ke / config.no_particles
    return 2/3 * avgke/kb


def GetEnergy(xyzfile, no_configurations, type):
    configs = InitConfigs(xyzfile, no_configurations)
    for config in configs:
        pot = Potential(config, type)
        ke = Kinetic(config)
        energy = pot + ke
        temp = CheckTemp(ke, config)
    print('Energies for ' + type + ' model')
    print('PE: ' + str(pot) + 'J')
    print('KE: ' + str(ke) + 'J')
    print('E: ' + str(energy) + 'J')
    print('T: ' + str(temp) + 'K')
    return

GetEnergy('conf.xyz', no_configs, 'PHS')