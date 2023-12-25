import algorithms as alg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
import matplotlib.lines as lines
from IPython.display import display, clear_output
import time

wall_on = False
choice = input('Wall on(1) or of(0): ')
if input == 1:
    wall_on = True

print(wall_on)

# Parametri
m = 20.0  # masa tela
k = 0.1  # koeficijent opruge
mu_s = 0.02  # koeficijent statičkog trenja
mu_d = 0.01  # koeficijent dinamičkog trenja
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

wall_position = -100.0

if wall_on:    
    wall_position = -105.0


fig, ax = plt.subplots(figsize=(10, 500/80))

#bez zida--------------------------------------------------------------
if not wall_on:
    body = patches.Circle((x[0], 1), 50.0, color='red', label='Telo')

    spring = lines.Line2D([-1350, x[0]], [1, 1], color='blue', linewidth=2, label='Opruga')

    wall = lines.Line2D([wall_position, wall_position], [0, 500], color='black', linestyle='--', label='Zid', linewidth=2)
#-----------------------------------------------------------------------


#zid--------------------------------------------------------------
if wall_on:
    body = patches.Circle((xz[0], 1), 30.0, color='red', label='Telo')

    spring = lines.Line2D([-1350, xz[0]], [1, 1], color='blue', linewidth=2, label='Opruga')

    wall = lines.Line2D([wall_position, wall_position], [-500, 500], color='black', linestyle='--', label='Zid', linewidth=2)
#-----------------------------------------------------------------------

ax.add_patch(body)
ax.add_line(spring) 
ax.add_line(wall)


ax.set_aspect('equal', adjustable='box')

if wall_on:
    ax.set_xlim([-500, 1400])

if not wall_on:
    ax.set_xlim([-1400, 1400])

ax.set_ylim([-300, 800])
ax.legend()

# Funkcija za animaciju
def update(frame):
    x_val = x[frame]
    spring.set_xdata([-1350, x_val])
    body.center = (x_val, 1)
    wall.set_xdata([-1350, -1350])
    
    time.sleep(0.1)
    
    clear_output(wait=True)
    display(fig)
    
    return body, spring, wall

def update_wall(frame):
    x_val = xz[frame]
    spring.set_xdata([-1350, x_val])
    body.center = (x_val, 1)
    wall.set_xdata([-110, -110])
    
    time.sleep(0.1)
    
    clear_output(wait=True)
    display(fig)
    
    return body, spring, wall


# Postavljanje animacije
if not wall_on:
    animation = FuncAnimation(fig, update, frames=len(t), interval=0.01, blit=True, repeat=False)

if wall_on:
    animation = FuncAnimation(fig, update_wall, frames=len(t), interval=0.01, blit=True, repeat=False)

#cuvanje animacije
#print('started...')
#if not wall_on:
#    animation.save(r'C:\Users\milan\Desktop\projects\spring_deformation_differential_modeling\spring_deformation_differential_modeling\animacije_sistema\kretanje_tela.gif',  writer='pillow')

#if wall_on:
#    animation.save(r'C:\Users\milan\Desktop\projects\spring_deformation_differential_modeling\spring_deformation_differential_modeling\animacije_sistema\kretanje_tela_zid.gif',  writer='pillow')
#print('finished')

# Prikaz animacije
display(fig)

# Ručno sinhronizovanje animacije
plt.pause(20)
#plt.pause(5)