import numpy as np
import matplotlib.pyplot as plt

#Initial variables
k = 1
m = 1
x0 = 0
v0 = 1
t0 = 0
dt = 0.3

nsteps = 20000
omega = np.sqrt(k/m)



def Algorithm(type):
    x = x0
    v = v0
    t = t0
    f = CalcForce(x, k)
    e = CalcEnergy(x, k, type) 



    x_values = [x]
    p_values = [v*m]
    f_values = [f]
    e_values = [e]
    t_values = [t]
    for i in range(nsteps):
        last_x = 0
        if len(x_values) > 1:
            last_x = x_values[-2]

        next_x = CalcPosition(x, last_x, v, dt, f, type)
        next_f = CalcForce(next_x, k) #Get next force
        next_v = CalcVelocity(last_x, next_x, v, dt, f, next_f, type)
        next_e = CalcEnergy(next_x, k, type) #Get next energy
        next_t = t + dt #Next time

        x_values.append(next_x)
        p_values.append(next_v*m)
        f_values.append(next_f)
        e_values.append(next_e)
        t_values.append(next_t)

        f = next_f
        x = next_x
        v = next_v
    
    plt.plot(x_values, p_values)
    plt.show()

    plt.plot(e_values, t_values)
    plt.show()

    

def CalcPosition(x, last_x, v, dt, f, type):
    if type == 'euler' or type == 'vverlet':
        return x + v*dt + 0.5 * (f/m) * dt**2  #Taylor expansion for position
    elif type == 'pverlet':
        return 2 * x - last_x + 1/m * f * dt**2 #Position verlet expansion


def CalcVelocity(last_x, next_x, v, dt, f, next_f, type):
    if type == 'euler':
        return v + 1/m * f * dt #Taylor expansion for velocity
    elif type == 'pverlet':
        return 1/(2*dt) * (next_x - last_x) #Position verlet velocity
    elif type == 'vverlet':
        return v + dt/(2*m) * (f + next_f) #Velocity verlet velocity

def CalcForce(x, k):
    return -k*x

def CalcEnergy(x, k, type):
    if type == 'euler':
        return 0
    else:
        return 0.5 * k * x**2

Algorithm('vverlet')