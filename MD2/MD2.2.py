import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###2.2###
no_configs = 50


#Configuration class to represent each configuration, require identity of particle as data option
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = 4500
        self.no_O = 1500
        self.no_H = 3000
        self.coords_O = None
        self.coords_H = None
        self.type = None
        self.V = None
        self.a = None

#Create confiiguration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines =[]
    for line in file:
        lines.append(line[:-1])
    file.close()
    
    configurations = []

    for i in range(no_configurations):
        #Initialise config
        config = configuration()
        config.timestep = int(lines[i * (config.no_particles + 9) + 1].split()[0])

        #Get volume of box
        a = np.zeros((3))
        for axis in range(3):
            bound = lines[5 + axis + (config.no_particles + 9) * i].split()
            a[axis] = float(bound[1]) - float(bound[0])
            vol  = np.prod(a)

        config.V = vol
        config.a = a

        #Empty arrays for type of each atom and coordinates of H and O atoms 
        config.type = np.zeros((config.no_particles))
        config.coords_O = np.zeros((config.no_O, 3))
        config.coords_H = np.zeros((config.no_H, 3))

        #Initialise coordinate arrays for O and H atoms
        O_index = 0
        H_index = 0
        for j in range(config.no_particles):
            config.type[j] = lines[i * (config.no_particles + 9) + 9 + j].split()[2]
            if config.type[j] == 1:
                config.coords_O[O_index] = lines[i * (config.no_particles + 9) + 9 + j].split()[4:7]
                O_index += 1
            elif config.type[j] == 2:
                config.coords_H[H_index] = lines[i * (config.no_particles + 9) + 9 + j].split()[4:7]
                H_index += 1

        #Scale coordinates into Angstrom
        config.coords_H *= a
        config.coords_O *= a

        configurations.append(config)
    
    return configurations


#Calculate RDF for a configuration
def RDF(config, r_edges, r_centres, dr, type):
    if type == 'OO':
        n_atoms = config.no_O ** 2
        rij = np.linalg.norm(PBConditions(config.coords_O[:, np.newaxis, :] - config.coords_O, config.a), axis = 2)
        np.fill_diagonal(rij, np.nan)
    elif type == 'OH':
        n_atoms = config.no_H * config.no_O
        rij = np.linalg.norm(PBConditions(config.coords_O[:, np.newaxis, :] - config.coords_H, config.a), axis = 2)
    n = np.histogram(rij, bins = r_edges)[0]
    g = n * config.V / (4 * np.pi * (r_centres + 0.5 * dr) ** 2 * dr * n_atoms)
    return g

def PBConditions(mat, box_length):
    return mat % (box_length / 2)



#Function to calculate RDF for every configuration from the XYZ file and average them, OH for OH RDF or OO for OO RDF
def GetAverageRDF(filename, maxr, dr, no_configurations, type):
    configs = InitConfigs(filename, no_configurations)

    r_edges = np.arange(0, maxr, dr)
    r_centres =  0.5 * (r_edges[:-1] + r_edges[1:])
    g_mat = np.zeros((no_configurations, len(r_edges) - 1))

    for config in range(len(configs)):
        g_mat[config] = RDF(configs[config], r_edges, r_centres, dr, type)
        print('Config ' + str(config+1) + ' of '+ str(no_configurations) + ' Complete' )

    g_average = np.sum(g_mat, axis = 0) / no_configurations
    plt.plot(r_centres, g_average)
    plt.xlabel('r / Angstrom')
    plt.ylabel('Average g(r)')
    plt.title(type + ' RDF')
    plt.show()

