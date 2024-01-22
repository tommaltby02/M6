import numpy as np
import matplotlib.pyplot as plt

#Programme to calculate energy of conf.xyz using different models

#Constants
kb = 1.38E-23
temp = 179.81
mass = 6.63E-26

#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
        self.particle_types = []
        self.velocities = None

#Create confiiguration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    #Get volume of box
    x = lines[5].split()
    a = float(x[1]) - float(x[0])
    V = a**3
    
    configurations = []

    for i in range(no_configurations):
        config = configuration()
        config.timestep = int(lines[i * 2009 + 1])
        config.no_particles = int(lines[i * 2009 + 3])
        config.coords = np.zeros((config.no_particles, 3))
        for j in range(config.no_particles):
            xyzline = lines[i * 2009 + 9 + j].split()
            xyz = np.array([xyzline[4], xyzline[5], xyzline[6]])
            config.coords[j] = xyz
            config.particle_types.append(xyzline[2])
        
        #Scale coordinates into Angstrom
        config.coords *= a
        
        #Generate velocities
        config = GenerateVelocities(config)

        #Get histogram
        GenerateVelHistogram(config)


        configurations.append(config)

    

    return configurations, V

#Use normal distribution with STD derived from Maxwell-Boltzmann distribution
def GenerateVelocities(config):

    std = np.sqrt(kb * temp / mass)

    config.velocities = np.random.normal(0, std, size = (config.no_particles, 3))

    return config

#Generate histogram of velocities
def GenerateVelHistogram(config):
    for i in range(3):
        plt.hist(config.velocities[:, i], bins = 50, density = True, alpha = 0.7, label = f'Axis {i+1}')
        plt.xlabel('Velocity')
        plt.ylabel('Prob Density')
        plt.legend()
        plt.show()


def Potential(config, type):
    Utot = 0
    for i in range(config.no_particles - 1):
        for j in range(i + 1, config.no_particles):
            r = np.linalg.norm(config.coords[i] - config.coords[j])
            if type == 'LJ':
                Utot += GetLJPairEnergy(r)
            elif type == 'PHS':
                Utot += GetPHSPairEnergy(r)
            elif type == 'HSYDH':
                Utot += GetHSPairEnergy(r) + GetYDHPairEnergy(r, config.particle_types[i], config.particle_types[j])

    return Utot

#Calculate pair potential using Lennard-Jones model
def GetLJPairEnergy(r):
    #Constants
    sigma = 3.405 
    e = 119.87
    if r < 3 * sigma:
        return 4 * e * kb * ((sigma/r)**12 - (sigma/r)**6)
    else:
        return 0
    
#Calculate pair potential using Pseudo Hard Sphere model
def GetPHSPairEnergy(r):
    #Constants
    sigma = 3.405
    lamr = 50
    lama = 49
    er = 119.87
    if r < (lamr/lama) * sigma:
        return lamr * (lamr / lama) ** lama * er * kb * ((sigma/r)**lamr - (sigma/r)**lama) + er * kb
    else:
        return 0
    
#Calculate pair potential using Hard Sphere model
def GetHSPairEnergy(r):
    #Constants
    sigma = 3.05
    
    if r < sigma:
        return float('inf')
    else:
        return 0

#Calculate pair potential using Yukawa Debye Huckel model
def GetYDHPairEnergy(r, type_a, type_b):
    #Constants
    sigma = 3.05
    e = 119.87
    k = 5

    if r < 3.5 * sigma:
        return  GetSign(type_a, type_b) * e * kb * sigma / r * np.exp(-k * (r - sigma))
    else:
        return 0

    

def GetSign(type_a, type_b):
    if type_a == type_b:
        return 1
    else:
        return -1


def Kinetic(config):
    ke = 0
    for i in range(config.no_particles):
        speed = np.linalg.norm(config.velocities[i])
        ke += GetKE(speed)
    return ke 


def GetKE(speed):
    return 0.5 * mass * speed **2

def CheckTemp(ke, config):
    avgke = ke / config.no_particles
    return 2/3 * avgke/kb


def GetEnergy(xyzfile, no_configurations, type):
    configs, V = InitConfigs(xyzfile, no_configurations)

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

GetEnergy('conf.xyz', 1, 'HSYDH')