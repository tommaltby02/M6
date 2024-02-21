import numpy as np
import matplotlib.pyplot as plt
#Program to calculate MSD 

#Constants
m = 6.63E-26 #kg
e = 6.047E-20 #J
sigma = 3.405 #Angstrom
tau = np.sqrt(m * (sigma * 1E-10)**2  / e) #s
no_configs = 50



#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
        self.V = None
        self.a = None
        self.id = None

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
        config.id = np.zeros((config.no_particles, 1))

        #Get volume of box
        a = np.zeros((3))
        for axis in range(3):
            x = lines[5 + (config.no_particles + 9) * i + axis].split()
            a[axis] = float(x[1]) - float(x[0])
        vol = np.prod(a)

        config.a = a
        config.V = vol

        for j in range(config.no_particles):
            config.coords[j] = lines[i * (config.no_particles + 9) + 9 + j].split()[4:7]
            config.id[j] = int(lines[i * (config.no_particles + 9) + 9 + j].split()[0])

        configurations.append(config)
    return FixMatrices(configurations)

#Get COM drift between the frames
def GetCOM(config1, config2):
    com = np.mean(config2.coords, axis = 0) - np.mean(config1.coords, axis = 0)
    return com

def GetMSD(config1, config2):
    msd = np.sum((config2.coords - config1.coords - GetCOM(config1, config2))**2) / config2.no_particles
    return msd / sigma**2 #In units of sigma

#XYZ file is written such that the order of particle identity is not consistent between frames, fix matrices to be in same order as first frame
def FixMatrices(configs):
    master_id = configs[0].id
    for i in range(1,len(configs)):
        new_mat = np.zeros((configs[i].no_particles,3))
        for id in range(len(master_id)):
            current_id = master_id[id][0]
            index = np.where(configs[i].id == current_id)[0][0]
            new_mat[id] = configs[i].coords[index]
        configs[i].coords = new_mat
    return configs


def MSD(filename, no_configurations, frame_skip):
    configs = InitConfigs(filename, no_configurations)

    msd = []
    t = []

    n = 0
    for i in range(1, len(configs)):
        if n % frame_skip == 0:
            msd.append(GetMSD(configs[0], configs[i]))
            t.append(configs[i].timestep * 10E-15 / tau) #In units of tau
        n += 1
    plt.plot(t, msd)
    plt.xlabel("t / tau")
    plt.ylabel("MSD / $σ^2$")
    plt.title(filename[10:-4])
    plt.show()



MSD('XYZ Files/diffusionA.xyz', no_configs, 1)