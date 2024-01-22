import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Constants
sigma = 3.4
no_confs = 40

###6.2###

#Configuration class to contin configuration info, here we only care about polymer beads so polymer identity is irrelevant
class config:
    def __init__(self):
        self.timestep = None
        self.no_atoms = None
        self.coords = None

def InitConfigs(xyzfile, no_configurations):

    #Get file data
    file = open(xyzfile, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    #Get box size
    a, V = GetBoxVolume(lines)

    #Initiate polymers for all configs
    confs = []

    #Initiate coords of every atom in every config
    for i in range(no_configurations):
        conf = config()
        conf.no_atoms = int(lines[3])
        conf.timestep = int(lines[1 + i * (conf.no_atoms + 9)])
        conf.coords = np.zeros((conf.no_atoms, 3))
        for j in range(conf.no_atoms):
            conf.coords[j] = lines[9 + i * (conf.no_atoms + 9) + j].split()[3:6]
        
        #Scale coordinates into Angstrom
        for axis in range(3):
            conf.coords[:,axis] *= a[axis]


        confs.append(conf)
    
    return confs, a , V 


#Function to get box volume as system is not a cube
def GetBoxVolume(lines):
    xyz = np.zeros((3,2))
    for i in range(3):
        xyz[i] = np.array(lines[5 + i].split())  
    a = xyz[:,1] - xyz[:,0]
    V = np.prod(a)
    return a, V



def GetConfigDensity(config, bins, v):
    num_density = []

    for bin_count in range(len(bins) - 1):
        n = GetNumberParticles(config, bins[bin_count], bins[bin_count + 1])
        density = n * sigma**3 / v
        num_density.append(density)
    return num_density

def PlotDensity(num_density, bin_list, name):
    plt.plot(bin_list, num_density)
    plt.title(name)
    plt.xlabel('z')
    plt.ylabel('Number Density')
    plt.show()


def GetNumberParticles(config, bin_start, bin_end):
    df = pd.DataFrame({'z': config.coords[:,2]})
    n = sum((df.z >= bin_start) & (df.z < bin_end))
    return n


def GetDensityProfile(xyzfile, no_confs, no_bins):
    confs, a, V = InitConfigs(xyzfile, no_confs)

    dz = a[2] / no_bins
    v = dz * a[0] * a[1]
    bins = np.arange(0, a[2], dz)


    density_mat = np.zeros((no_confs, no_bins - 1))

    for config in range(len(confs)):
        num_density = GetConfigDensity(confs[config], bins, v)
        density_mat[config] = num_density
    
    average_density = np.sum(density_mat, 0) / no_confs
    average_bin = []
    for i in range(no_bins - 1):
        average_bin.append((bins[i] + bins[i+1])/2)

    PlotDensity(average_density, average_bin, xyzfile[:-4])

    return
    







GetDensityProfile('300K.xyz', no_confs, 40)