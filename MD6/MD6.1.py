import numpy as np
import matplotlib.pyplot as plt

###6.1###
#Program calculates all contacts for all temps in ~2mins

#Constants
no_configurations = 40
no_pols = 90
rc = 4.1
no_atoms = 3600
temp_list = [264, 276, 288, 300]


#Configuration class to contin configuration info
class config:
    def __init__(self):
        self.no_polymers = None
        self.polymers = []
        self.timestep = None
        self.no_atoms = None
        self.V = None
        self.a = None
        self.coords = None
        self.overlap_mat = None

#Initiate configurations, each configuration contains a list of all polymers in the system, the polymer objects contain the coordinates of the respective polymer segments
def InitConfigs(xyzfile, no_configs, no_polymers):

    #Get file data
    file = open(xyzfile, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    file.close()

    #Initiate all configs
    confs = []

    #Get overlap matrix
    overlap_mat = GetOverlapMatrix(no_atoms, no_pols)

    #For every frame in the xyz file
    for i in range(no_configs):
        conf = config()
        conf.no_polymers = no_polymers
        conf.no_atoms = no_atoms

        #Get box volume
        a = np.zeros((3))
        for axis in range(3):
            x = lines[5 + (conf.no_atoms + 9) * i + axis].split()
            a[axis] = float(x[1]) - float(x[0])
        vol = np.prod(a)

        conf.a = a
        conf.V = vol

        conf.coords = np.zeros((conf.no_atoms, 3))
        for j in range(conf.no_atoms):
            conf.coords[j] = lines[i * (conf.no_atoms + 9) + 9 + j].split()[3:6]
        
        conf.coords *= a
        conf.overlap_mat = overlap_mat

        #Add config to list of configs
        confs.append(conf)
    
    return confs

#Generate overlap matrix to account for directly joined beads in each polymer, this matrix multiplied element wise by the rij matrix containg the 
#distances between every bead will account for removing the close contacts
def GetOverlapMatrix(no_atoms, no_pols):
    overlap_mat = np.ones((no_atoms, no_atoms))
    no_beads = int(no_atoms / no_pols)
    mat = np.ones((no_beads, no_beads))
    
    #Generate mini matrix which has 1s everywhere but i = j + 1, i = j - 1, these elements have 0
    for bead1 in range(no_beads - 1):
        for bead2 in range(bead1, no_beads):
            if bead2 == bead1 + 1 or bead2 == bead1 - 1:
                mat[bead1][bead2] = 0
                mat[bead2][bead1] = 0
    #Modify overlap_matrix to have the mini matrix where polymer[i] = polymer[i]
    for i in range(no_pols):
        overlap_mat[i * no_beads : no_beads + i * no_beads , i * no_beads : no_beads + i * no_beads] = mat
    return overlap_mat

#Get all contacts at once
def GetContacts(conf):
    rij = np.linalg.norm(conf.coords[:, np.newaxis, :] - conf.coords, axis = 2) 
    rij_overlap = rij * conf.overlap_mat #Disregards all nearest neighbour distances
    return int(np.sum(np.where((rij_overlap > 0) & (rij_overlap < rc), 1, 0)) / 2) #Divide by 2 to account for double counting of interactions

#Function to average contacts over all polymers in all configs
def AverageOverConfigs(filename, no_confs, no_polymers):
    confs = InitConfigs(filename, no_confs, no_polymers)
    msum = 0
    for i in range(len(confs)):
        n = GetContacts(confs[i])
        print("Contacts for conformer " + str(i+1) + " = " + str(n))
        msum += n / confs[i].no_polymers
    
    average_contacts = msum / no_confs
    print("Average contacts at " + filename[:-4] + " = " + str(average_contacts))
    return average_contacts
        
def GetContactsvsTemp(temps, no_confs, no_polymers):
    average_n = []
    for temp in temps:
        average_n.append(AverageOverConfigs(str(temp) + "K.xyz", no_confs, no_polymers))
    
    plt.plot(temps, average_n)
    plt.xlabel("T / K")
    plt.ylabel("Average contacts per polymer")
    plt.show()



GetContactsvsTemp(temp_list, no_configurations, no_pols)



