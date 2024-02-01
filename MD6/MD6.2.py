import numpy as np
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
        self.a = None
        self.V = None

def InitConfigs(xyzfile, no_configurations):

    #Get file data
    file = open(xyzfile, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    #Initiate polymers for all configs
    confs = []

    #Initiate coords of every atom in every config
    for i in range(no_configurations):
        conf = config()
        conf.no_atoms = int(lines[3])
        conf.timestep = int(lines[1 + i * (conf.no_atoms + 9)])
        conf.coords = np.zeros((conf.no_atoms, 3))

        #Get box volume
        a = np.zeros((3))
        for axis in range(3):
            x = lines[5 + (conf.no_atoms + 9) * i + axis].split()
            a[axis] = float(x[1]) - float(x[0])
        vol = np.prod(a)
        conf.a = a
        conf.V = vol

        for j in range(conf.no_atoms):
            conf.coords[j] = lines[9 + i * (conf.no_atoms + 9) + j].split()[3:6]
        
        #Scale coordinates into Angstrom
        conf.coords *= a


        confs.append(conf)
    
    return confs


def GetConfigDensity(config, z_edges, vol):
    z = config.coords[:,2] 
    return np.histogram(z, bins = z_edges)[0] / vol

def GetDensityProfile(xyzfile, no_confs, no_bins):
    confs = InitConfigs(xyzfile, no_confs)

    #Get relevant z parameters
    z = confs[0].a
    dz = z[2] / no_bins
    vol = dz * z[0] * z[1]
    z_edges = np.arange(0, z[2], dz)
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])

    average_density = np.zeros((no_bins - 1))

    for conf in confs:
        average_density += GetConfigDensity(conf, z_edges, vol)
    
    average_density /= no_confs

    PlotDensity(average_density, z_centres, xyzfile)


    return
    
def PlotDensity(num_density, bin_list, name):
    plt.plot(bin_list, num_density)
    plt.title(name[:-4])
    plt.xlabel('z')
    plt.ylabel('Number Density')
    plt.show()



GetDensityProfile('300K.xyz', no_confs, 60)