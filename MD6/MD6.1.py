import numpy as np

#Constants
no_configurations = 40
no_pols = 90
rc = 4.1

###6.1###

#Configuration class to contin configuration info
class config:
    def __init__(self):
        self.no_polymers = None
        self.polymers = []
        self.timestep = None
        self.no_atoms = None

#Configuration class to contain polymer info
class polymer:
    def __init__(self):
        self.no_segments = None
        self.coords = None
        self.no = None

#Initiate configurations, each configuration contains a list of all polymers in the system, the polymer objects contain the coordinates of the respective polymer segments
def InitConfigs(xyzfile, no_configs, no_polymers):

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

    #For every frame in the xyz file
    for i in range(no_configs):
        conf = config()
        conf.no_polymers = no_polymers
        conf.no_atoms = int(lines[3])

        #For every polymer in the frame
        for j in range(no_polymers):
            pol = polymer()
            pol.no_segments = int(conf.no_atoms / no_polymers)
            pol.no = int(lines[9 + (conf.no_atoms + 9) * i + j * pol.no_segments].split()[1])
            pol.coords = np.zeros((pol.no_segments, 3))

            #For every segment of the polymer
            for k in range(pol.no_segments):
                coordline = lines[9 + (conf.no_atoms + 9) * i + j * pol.no_segments + k ].split()
                pol.coords[k] = np.array([coordline[3], coordline[4], coordline[5]])
            
            
            #Scale polymer coordinates into Angstrom
            for axis in range(3):
                pol.coords[:,axis] *= a[axis]
            
            #Add polymer to current config
            conf.polymers.append(pol)
        
        #Add config to list of configs
        confs.append(conf)
    
    return confs, a, V

#Function to get box volume as system is not a cube
def GetBoxVolume(lines):
    xyz = np.zeros((3,2))
    for i in range(3):
        xyz[i] = np.array(lines[5 + i].split())  
    print(xyz)
    a = xyz[:,1] - xyz[:,0]
    V = np.product(a)
    return a, V

#Functions to get intermolecular contacts, for every pair of polymers every unique pair of segments are compared, if pair are within rc contact is counted
def GetInterContacts(polymers):
    n = 0
    for polymerA in range(len(polymers) - 1):
        for polymerB in range(polymerA + 1, len(polymers)):
            n += Get2PolymerContacts(polymers[polymerA], polymers[polymerB])
    return n

def Get2PolymerContacts(polymerA, polymerB):
    n = 0
    for segmentA in range(polymerA.no_segments):
        for segmentB in range(polymerB.no_segments):
            r = np.linalg.norm(polymerA.coords[segmentA] - polymerB.coords[segmentB])
            if r < rc:
                n += 1
    return n


#Function to get intramolecular contacts, for every polymer every unique pair of segments within the same polymer are compared, if pair are within rc and not adjacent contact is counted
def GetIntraContacts(polymers):
    n = 0
    for polymer in polymers:
        for segmentA in range(polymer.no_segments - 1):
           for segmentB in range(segmentA + 1, polymer.no_segments):
            r = np.linalg.norm(polymer.coords[segmentA] - polymer.coords[segmentB])
            if segmentB != segmentA + 1 and segmentB != segmentA - 1 and r < rc:
                n += 1
    return n

def GetConfigContacts(conf):
    ninter = GetInterContacts(conf.polymers)
    nintra = GetIntraContacts(conf.polymers)
    n = nintra + ninter
    return n 

#Function to average contacts over all polymers in all configs
def AverageOverConfigs(filename, no_confs, no_polymers):
    confs, a, V = InitConfigs(filename, no_confs, no_polymers)
    msum = 0

    for conf in confs:
        msum += GetConfigContacts(conf) / conf.no_polymers
    
    return msum / no_confs
        



AverageOverConfigs('264K.xyz', no_configurations, no_pols)
        

