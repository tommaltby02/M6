import numpy as np
import matplotlib.pyplot as plt

def MetropolisAlgorithm(MCsteps, initial_nj, be, change_type):    
    currentnj = Initialisenj(initial_nj)
    sumnj = 0 
    for i in range(MCsteps):
        rand_num = np.random.uniform(0,1) - 0.5
        trialnj, a = GetTrialnj(rand_num, currentnj, change_type)
        if np.random.uniform(0,1) < np.exp(-a * be):
            currentnj = trialnj
        sumnj += currentnj 
    return sumnj / MCsteps

def Initialisenj(initial_nj):
    if initial_nj == "1":
        return 1
    elif initial_nj =="rand":
        return np.random.randint(0,50)

def GetTrialnj(rand_num, currentnj, change_type):
    if change_type == "1":
        if rand_num > 0:
            return currentnj + 1, 1
        else:
            if currentnj == 0:
                return 0, 0
            else:
                return currentnj - 1, -1
    elif change_type == "3":
        if rand_num > 0:
            return currentnj + 3, 3
        else:
            if currentnj - 3 < 0:
                return 0, 0
            else:
                return currentnj - 3, -3
    elif change_type == "rand5":
        rand_int = np.random.randint(-5, 6)
        if currentnj + rand_int < 0:
            return 0, 0
        else:
            return currentnj + rand_int, rand_int

def AnalyticSolution(be):
    return 1 / (np.exp(be) - 1)

def RunAlgorithm(MCsteps, initial_nj, change_type, no_runs):
    be_list = np.arange(0.1,2,0.05)
    nj_analytical = []
    nj_numerical = np.zeros((len(be_list), no_runs))

    for i, be in enumerate(be_list):
        for j in range(no_runs):
            nj_numerical[i][j] = MetropolisAlgorithm(MCsteps, initial_nj, be, change_type)
        nj_analytical.append(AnalyticSolution(be))
    
    average_numerical = np.sum(nj_numerical, axis = 1) / no_runs
    errors = np.abs(average_numerical - nj_analytical)

    plt.plot(be_list, average_numerical, label = "Numerical")
    plt.plot(be_list, nj_analytical, label = "Analytical")
    plt.xlabel("βε")
    plt.ylabel("<n>")
    plt.legend()
    plt.text(1.5, 7, f"MAE = {np.around(np.mean(errors), 3)}")
    plt.title(f"Initial nj = {initial_nj}, nj changed by {change_type}, averaged over {no_runs} runs of {MCsteps} steps")
    plt.show()
    return 
    


RunAlgorithm(10000, "rand", "rand5", 50)


    