import algorithms as alg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math


# Parametri
m = 20.0  # masa tela
k = 0.1  # koeficijent opruge
mu_s = 0.5  # koeficijent statičkog trenja
mu_d = 0.1  # koeficijent dinamičkog trenja
g = 9.81  # ubrzanje gravitacije

# Vremenski korak i vremenski opseg
dt = 0.5
t = np.arange(0, 1200, dt)

# Početna brzina i početna pozicija tela
v0 = 100.0
x0 = 1.0

v, x = alg.runge_kutta_method(alg.differential_equation, v0, x0, t, m, k, mu_s, mu_d, g)

v_e, x_e = alg.euler_method(alg.differential_equation, v0, x0, t, m, k, mu_s, mu_d, g)

alg.change_wall_status()

vz, xz = alg.runge_kutta_method(alg.differential_equation, v0, x0, t, m, k, mu_s, mu_d, g)

v_ez, x_ez = alg.euler_method(alg.differential_equation, v0, x0, t, m, k, mu_s, mu_d, g)

kinetic_energy = 0.5 * m * v**2
potential_energy = 0.5 * k * x**2
total_energy = kinetic_energy + potential_energy

kinetic_energy_e = 0.5 * m * v_e**2
potential_energy_e = 0.5 * k * x_e**2
total_energy_e = kinetic_energy_e + potential_energy_e

kinetic_energy_w = 0.5 * m * vz**2
potential_energy_w = 0.5 * k * xz**2
total_energy_w = kinetic_energy_w + potential_energy_w

kinetic_energy_ew = 0.5 * m * v_ez**2
potential_energy_ew = 0.5 * k * x_ez**2
total_energy_ew = kinetic_energy_ew + potential_energy_ew

# Inicijalizacija grafikona
fig, ax = plt.subplots(figsize=(8, 4))

# Inicijalizacija podataka
x_data1, y_data1 = [], []
x_data2, y_data2 = [], []

# Funkcija za inicijalizaciju
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2

# Funkcija za ažuriranje animacije
def update(frame):
    x_data1.append(t[frame])
    y_data1.append(x[frame])
    
    x_data2.append(t[frame])
    y_data2.append(x_e[frame])
    
    line1.set_data(x_data1, y_data1)
    line2.set_data(x_data2, y_data2)
    
    return line1, line2

def update_wall(frame):
    x_data1.append(t[frame])
    y_data1.append(xz[frame])
    
    x_data2.append(t[frame])
    y_data2.append(x_ez[frame])
    
    line1.set_data(x_data1, y_data1)
    line2.set_data(x_data2, y_data2)
    
    return line1, line2

def update_comp(frame):
    x_data1.append(t[frame])
    y_data1.append(x[frame])
    
    x_data2.append(t[frame])
    y_data2.append(xz[frame])
    
    line1.set_data(x_data1, y_data1)
    line2.set_data(x_data2, y_data2)
    
    return line1, line2

def update_energy(frame):
    x_data1.append(t[frame])
    y_data1.append(total_energy[frame])
    
    x_data2.append(t[frame])
    y_data2.append(total_energy_e[frame])
    
    line1.set_data(x_data1, y_data1)
    line2.set_data(x_data2, y_data2)
    
    return line1, line2

def update_energy_wall(frame):
    x_data1.append(t[frame])
    y_data1.append(total_energy_w[frame])
    
    x_data2.append(t[frame])
    y_data2.append(total_energy_ew[frame])
    
    line1.set_data(x_data1, y_data1)
    line2.set_data(x_data2, y_data2)
    
    return line1, line2



# Kreiranje i cuvanje animacija------------------------------------------------------------------------------------
line1, = ax.plot([], [], label='RK4')

ax.set_xlabel('Vreme (s)')
ax.set_ylabel('Pozicija tela')
ax.set_title('Kretanje tela u odnosu na oprugu')
ax.legend()
ax.grid(True)
ax.set_xlim(0, 1200)
ax.set_ylim(-1500, 1500)

line2, = ax.plot([], [], label='Ojler', color='red')
ax.legend()

#ani = FuncAnimation(fig, update, frames=len(t), init_func=init, blit=True, interval=2.0)
#ani.save(r'C:\Users\milan\Desktop\projects\projekat_a4\kretanje_tela_bez_zida.gif',  writer='pillow') #ovaj vec imas

#ani1 = FuncAnimation(fig, update_wall, frames=len(t), init_func=init, blit=True, interval=2.0)
#ani1.save(r'C:\Users\milan\Desktop\projects\projekat_a4\kretanje_tela_sa_zidom.gif',  writer='pillow')

#ani_comp = FuncAnimation(fig, update_comp, frames=len(t), init_func=init, blit=True, interval=2.0)
#ani_comp.save(r'C:\Users\milan\Desktop\projects\projekat_a4\kretanje_comp.gif',  writer='pillow')

ax.set_xlabel('Vreme (s)')
ax.set_ylabel('Ukupna energija tela')
ax.set_title('Energija tela u odnosu na vreme')
ax.set_xlim(0, 1200)
ax.set_ylim(0, 100000)
ax.legend()

#ani2 = FuncAnimation(fig, update_energy, frames=len(t), init_func=init, blit=True, interval=2.0)
#ani2.save(r'C:\Users\milan\Desktop\projects\projekat_a4\energija_tela.gif',  writer='pillow')

#ani3 = FuncAnimation(fig, update_energy_wall, frames=len(t), init_func=init, blit=True, interval=2.0)
#ani3.save(r'C:\Users\milan\Desktop\projects\projekat_a4\energija_tela_zid.gif',  writer='pillow')

# Prikaz animacije
#plt.show()


