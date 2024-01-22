import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Programme to process XYZ file containing 50 different configurations and calculate the RDF for each one

#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None

#Create confiiguration class for every confirguration storing timestep and x, y, z coords of every particle
def InitConfigs(xyz_file, no_configurations):
    file = open(xyz_file, 'r')
    lines =[]
    for line in file:
        lines.append(line[:-1])
    file.close()
    
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
        configurations.append(config)
    

    #Get volume of box, here coordinates are scaled so the box volume is 1
    V = 1
     

    return configurations, V

#Calculate RDF for a configuration
def CalcRDF(config, r, dr, V):
    n_sum = 0
    for i in range(config.no_particles):
        n_in_shell = GetAtoms(r, dr, config, i)
        n_sum  += n_in_shell

    return V/(4* np.pi * r**2 * dr * config.no_particles**2) * n_sum

    
#Get number of atoms in shell defined by r and dr
def GetAtoms(r, dr, config, centre_particle):
    rijvec  = np.linalg.norm(config.coords - config.coords[centre_particle], axis = 1) 
    df = pd.DataFrame({'r': rijvec})
    n = sum((df.r < r) & (df.r >= r - dr))
    return n


#Function to calculate RDF for every configuration from the XYZ file and average them 
def GetAverageRDF(filename, maxr, dr, rstep, no_configurations, desired_config):
    configs, V = InitConfigs(filename, no_configurations)

    
    GetRDF(configs[desired_config], 0.005, maxr, dr, rstep, V)
    
    return

def GetRDF(config, minr, maxr, dr, rstep, V):
    glist = []
    rlist = np.arange(minr, maxr, rstep)
    num = 0
    for r in rlist:
        glist.append(CalcRDF(config, r, dr, V))
        num += 1
        print(str((num / len(rlist)) * 100) + "% Complete")
    plt.plot(rlist, glist)
    plt.show()




GetAverageRDF('ideal.xyz', 0.5, 0.001, 0.01, 50, 15)