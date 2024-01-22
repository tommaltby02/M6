import numpy as np
import matplotlib.pyplot as plt

possible_states = [1, -1]
only_up = [1]
only_down = [-1]
kB = 1.38E-23
def InitialiseSpins(lattice_size):
    spins = np.random.choice(only_up, (lattice_size, lattice_size))
    return spins

def MC(spins, no_steps, H, J, T):
    current_spins = np.array(spins)
    magnetism = []
    steps = []
    mag = CalcMagnetism(current_spins)
    magnetism.append(mag)
    steps.append(0)
    for i in range(no_steps):
        new_spins = MCAlgorithm(current_spins, H, J, T)
        mag = CalcMagnetism(new_spins)
        magnetism.append(mag)
        steps.append(i+1)
        current_spins = new_spins

    return magnetism, steps

def MCAlgorithm(spins, H, J, T):
    new_spins = np.array(spins)
    (i , j) = (np.random.randint(0, spins.shape[0]), np.random.randint(0, spins.shape[1]))
    
    new_spins[i][j] = -1 * new_spins[i][j]


  


    energy_change = CalcEnergyChange(new_spins, i, j, H, J, T)
    beta = 1/(kB * T)
    

    if np.random.uniform(0,1) < np.exp(-beta * energy_change):
        return new_spins
    else:
        return spins


def CalcMagnetism(spins):
    return 1/spins.shape[0]**2 * np.sum(spins)

def CalcEnergyChange(new_spins, i, j, H, J, T):
    nearest_sum = 0
    for k in range(-1,2):
        if k == 0:
            pass
        else:
            x_coord = 0
            y_coord = 0
            if i + k == new_spins.shape[0]:
                x_coord = 0
            elif i + k < 0:
                x_coord = new_spins.shape[0] - 1
            else:
                x_coord = i + k
            
            if j + k  == new_spins.shape[1]:
                y_coord = 0
            elif j + k < 0:
                y_coord = new_spins.shape[1] - 1
            else:
                y_coord = j +k 
            

            nearest_sum += new_spins[x_coord, j]
            nearest_sum += new_spins[i, y_coord]
    
    return -2 * new_spins[i,j] * (-J/(kB * T) * nearest_sum + H)
    
            
def RunSimulation(lattice_size, no_steps, J, H, T):
    initial_spins = InitialiseSpins(lattice_size)
    mag, steps = MC(initial_spins, no_steps, H, J, T)
    plt.plot(steps, mag)
    plt.show()
    return




RunSimulation(40, 10000, 0.2E-20, -2, 100000)