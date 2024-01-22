import numpy as np
import matplotlib.pyplot as plt
import datetime as datetime

#Programme to calculate how pressure and enthalpy evolves over 50 frames of configurations
no_configs = 50
m = 1

#Configuration class to represent each configuration
class configuration:
    def __init__(self):
        self.timestep = None
        self.no_particles = None
        self.coords = None
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
    
    #Generate configuration class
    for i in range(no_configurations):
        config = configuration()
        config.no_particles = int(lines[3])
        config.timestep = int(lines[(config.no_particles + 9) *i + 1])
        config.coords = np.zeros((config.no_particles, 3))
        config.velocities = np.zeros((config.no_particles, 3))
        for j in range(config.no_particles):
            xyzline = lines[i * (config.no_particles + 9) + 9 + j].split()
            coords = np.array([xyzline[4], xyzline[5], xyzline[6]])
            vels = np.array([xyzline[7], xyzline[8], xyzline[9]])
            config.coords[j] = coords
            config.velocities[j] = vels

        
        #Scale positions and vels
        config.coords *= a
        config.velocities *= a
        
        configurations.append(config)

    

    return configurations, V


def GetPressure(config, V, type):
    sum = 0
    for i in range(config.no_particles):
        sum += GetMomentum(config.velocities[i]) ** 2 / m + np.linalg.norm(config.coords[i]) * GetForce(config, i, type)
    return 1/(3 * V) * sum

def GetMomentum(vel):
    return np.linalg.norm(vel) * m

def GetForce(config, particleA, type):
    f = 0
    for particleB in range(config.no_particles):
        if particleB == particleA:
            pass
        else:
            r = np.linalg.norm(config.coords[particleA] - config.coords[particleB])
            if type == 'LJ':
                f += GetLJForce(r)
            elif type == 'PHS':
                f += GetPHSForce(r)
    return f
            
def GetLJForce(r):
    sigma = 3.405
    e = 0.1
    if r < 3 * sigma:
        return - 24 * e * ((sigma/r)**7 - 2 * (sigma/r)*13)  
    else:
        return 0
    
def GetPHSForce(r):
    sigma = 3.305
    e = 0.1
    lama = 49
    lamr = 50

    if r < (lamr/lama) * sigma:
        return - lamr * (lamr / lama) ** lama * e * (lama * (sigma / r) ** (lama + 1) - lamr * (sigma / r)**(lamr + 1) )
    else:
        return 0

def InstPressure(filename, no_configurations, type):
    configs, V = InitConfigs(filename, no_configurations)
    plist = []
    time = []

    for i in range(len(configs)):
        t1 = datetime.datetime.now()
        p = GetPressure(configs[i], V, type)
        plist.append(p)
        time.append(configs[i].timestep)
        t2 = datetime.datetime.now()
        duration = t2.second - t1.second
        print('Config ' + str(i+1) + ' of ' + str(len(configs)) + ' Complete' )
        print('Time taken ' + str(duration) + 's')
    
    plt.plot(time, plist)
    plt.show()


    return







InstPressure('pres.xyz', no_configs, 'LJ')