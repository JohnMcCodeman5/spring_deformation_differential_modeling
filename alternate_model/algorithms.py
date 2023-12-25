import numpy as np
import matplotlib.pyplot as plt
import math

Wall = False

#ovo cemo da koristimo da switchujemo verziju problema
def change_wall_status():
    global Wall
    Wall = not Wall


#def differential_equation_wall(v, x, t, m, k, mu_s, mu_d, g):
#    dxdt = v

#    dvdt = (-k * x - mu_s * np.sign(v) - mu_d * v) / m 

#    if x < 0:
#        e = 0.8   #(koeficijent restitucije)
#        v = -e * v
    
#    return dvdt, dxdt


def differential_equation(v, x, t, m, k, mu_s, mu_d, g):
    dxdt = v

    dvdt = -k * x / m - mu_s * g * np.sign(v) - mu_d * g * v 
    
    return dvdt, dxdt

def runge_kutta_method(f, v0, x0, t, m, k, mu_s, mu_d, g):
    h = t[1] - t[0]
    v = np.zeros_like(t)
    x = np.zeros_like(t)
    v[0] = v0
    x[0] = x0
    counter = 0

    for i in range(1, len(t)): 

        k1v, k1x = f(v[i-1], x[i-1], t[i-1], m, k, mu_s, mu_d, g)
        k1v, k1x = h * k1v, h * k1x

        k2v, k2x = f(v[i-1] + 0.5*k1v, x[i-1] + 0.5*k1x, t[i-1] + 0.5*h, m, k, mu_s, mu_d, g)
        k2v, k2x = h * k2v, h * k2x

        k3v, k3x = f(v[i-1] + 0.5*k2v, x[i-1] + 0.5*k2x, t[i-1] + 0.5*h, m, k, mu_s, mu_d, g)
        k3v, k3x = h * k3v, h * k3x
        
        k4v, k4x = f(v[i-1] + k3v, x[i-1] + k3x, t[i-1] + h, m, k, mu_s, mu_d, g)
        k4v, k4x = h * k4v, h * k4x

        v[i] = v[i-1] + (k1v + 2*k2v + 2*k3v + k4v) / 6
        x[i] = x[i-1] + (k1x + 2*k2x + 2*k3x + k4x) / 6

        if Wall:
            if x[i] < -100:
                x[i] = -100
                e = 0.65   #(koeficijent restitucije)
                v[i] = -e * v[i]

        if v[i] == v[i-1]:
            #print(counter)
            counter += 1
            if counter == 10:
                break

    return v, x


def euler_method(f, v0, x0, t, m, k, mu_s, mu_d, g):
    h = t[1] - t[0]
    v = np.zeros_like(t)
    x = np.zeros_like(t)
    v[0] = v0
    x[0] = x0

    for i in range(1, len(t)):
        dvdt, dxdt = f(v[i-1], x[i-1], t[i-1], m, k, mu_s, mu_d, g)
        v[i] = v[i-1] + dvdt * h
        x[i] = x[i-1] + dxdt * h

        if Wall:
            if x[i] < -100:
                x[i] = -100
                e = 0.65   #(koeficijent restitucije)
                v[i] = -e * v[i]

    return v, x


def modified_euler_method(f, v0, x0, t, m, k, mu_s, mu_d, g):
    h = t[1] - t[0]
    v = np.zeros_like(t)
    x = np.zeros_like(t)
    v[0] = v0
    x[0] = x0

    for i in range(1, len(t)):
        k1v, k1x = f(v[i-1], x[i-1], t[i-1], m, k, mu_s, mu_d, g)
        v_middle = v[i-1] + 0.5 * k1v * h
        x_middle = x[i-1] + 0.5 * k1x * h
        k2v, k2x = f(v_middle, x_middle, t[i-1] + 0.5 * h, m, k, mu_s, mu_d, g)

        v[i] = v[i-1] + k2v * h
        x[i] = x[i-1] + k2x * h

        if Wall:
            if x[i] < -100:
                x[i] = -100
                e = 0.65   #(koeficijent restitucije)
                v[i] = -e * v[i]

    return v, x    


def runge_kutta_alt_method(f, v0, x0, t, m, k, mu_s, mu_d, g):
    h = t[1] - t[0]
    v = np.zeros_like(t)
    x = np.zeros_like(t)
    v[0] = v0
    x[0] = x0
    counter = 0

    for i in range(1, len(t)): 

        k1v, k1x = f(v[i-1], x[i-1], t[i-1], m, k, mu_s, mu_d, g)
        k1v, k1x = h * k1v, h * k1x

        k2v, k2x = f(v[i-1] + k1v/3, x[i-1] + k1x/3, t[i-1] + h/3, m, k, mu_s, mu_d, g)
        k2v, k2x = h * k2v, h * k2x

        k3v, k3x = f(v[i-1] + k2v/3, x[i-1] + k2x/3, t[i-1] + h/3, m, k, mu_s, mu_d, g)
        k3v, k3x = h * k3v, h * k3x
        
        k4v, k4x = f(v[i-1] + k3v, x[i-1] + k3x, t[i-1] + h, m, k, mu_s, mu_d, g)
        k4v, k4x = h * k4v, h * k4x

        v[i] = v[i-1] + (k1v + 3*k2v + 3*k3v + k4v) / 8
        x[i] = x[i-1] + (k1x + 3*k2x + 3*k3x + k4x) / 8

        if Wall:
            if x[i] < -100:
                x[i] = -100
                e = 0.65   #(koeficijent restitucije)
                v[i] = -e * v[i]

        if v[i] == v[i-1]:
            #print(counter)
            counter += 1
            if counter == 10:
                break

    return v, x


