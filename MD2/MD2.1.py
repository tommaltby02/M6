import numpy as np
import matplotlib.pyplot as plt

#Programme to process XYZ file containing different configurations and calculate the average RDF

###2.1 and 2.3###

#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
        self.V = None
        self.a = None

#Create configuration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines =[]
    for line in file:
        lines.append(line[:-1])
    file.close()
    
    configurations = []
    for i in range(no_configurations):
        config = configuration()
        config.no_particles = int(lines[3])
        config.timestep = int(lines[i * (config.no_particles + 9) + 1])
        
        #Get volume and length of box
        a = np.zeros((3))
        for axis in range(3):
            bound = lines[5 + axis + (config.no_particles + 9) * i].split()
            a[axis] = float(bound[1]) - float(bound[0])
        vol  = np.prod(a)
        config.V = vol
        config.a = a
        
        #Get coordinates of every particle and store in (N, 3) matrix
        config.coords = np.zeros((config.no_particles, 3))
        for j in range(config.no_particles):
            config.coords[j] = lines[i * (config.no_particles + 9) + 9 + j].split()[4:7]
        config.coords *= a #Scale coords from box length units to Angstrom
        
        #Add config to list
        configurations.append(config)
    
    return configurations


#Calculate RDF for a configuration
def RDF(config, r_edges, r_centres, dr):
    rij = config.coords[:, np.newaxis, :] - config.coords
    #Generate (N,N) matrix containing the norm of all inter particle distances, also apply periodic boundary conditions
    rij = np.linalg.norm(PBConditions(rij, config.a), axis = -1)
    #To prevent self counting fill diagonal with nan
    np.fill_diagonal(rij, np.nan)
    #For each bin in r_edges count number of rij values that lie within the boundary
    n = np.histogram(rij, bins = r_edges)[0]
    return n  * config.V / (4 * np.pi * (r_centres + 0.5 * dr) ** 2 * dr * config.no_particles**2) #Return the RDF function calculated for every r_center

#Use modulo to get smallest distance between any two particles
def PBConditions(mat, box_length):
    mat -= np.rint(mat / box_length) * box_length
    return mat


#Function to calculate RDF for every configuration from the XYZ file and average them 
def GetAverageRDF(filename, maxr, dr, no_configurations):
    #Initialise all configurations with relevant parameters
    configs = InitConfigs(filename, no_configurations)

    #Generate bins based on dr
    r_edges = np.arange(0, maxr, dr)
    #Get centre of each r bin
    r_centres =  0.5 * (r_edges[:-1] + r_edges[1:])
    #Empty g matrix for averaging over configurations
    g_mat = np.zeros((no_configurations, len(r_edges) - 1))

    for config in range(len(configs)):
        g_mat[config] = RDF(configs[config], r_edges, r_centres, dr)
        print('Config ' + str(config+1) + ' of '+ str(no_configurations) + ' Complete' )

    #Get average RDF over all configurations
    g_average = np.sum(g_mat, axis = 0) / no_configurations
    
    plt.plot(r_centres, g_average)
    plt.xlabel('r / Angstrom')
    plt.ylabel('Average g(r)')
    plt.title(filename[:-4] + ' RDF')
    plt.show()

    
    return

GetAverageRDF("XYZ Files/set2.xyz", 20, 0.5, 50)





