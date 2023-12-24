import tkinter as tk
from tkinter import Label, Button, Entry
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from helper import differential_equation, runge_kutta_method, euler_method, modified_euler_method, runge_kutta_alt_method, modified_euler_alt
from tkinter import ttk

def plot_graph(v0, h,m, k, mu_d, mu_s):
    #m = 20.0
    #k = 0.1 
    #mu_s = 0.5  # koeficijent statičkog trenja
    #mu_d = 0.1  # koeficijent dinamičkog trenja
    g = 9.81  # ubrzanje usled gravitacije
    t = np.arange(0, 1100, h)
    x0 = 1.0
    
    wall = True
    if selected_wall.get() == "ne":
        wall = False
    
    # Ocistimo ako je nesto bilo pre
    #ax.clear()
    
    selected_m = selected_method.get()
    if selected_m == "RK4":
        v, x = runge_kutta_method(differential_equation, v0, x0, t, m, k, mu_s, mu_d, g, wall)
    elif selected_m == "Ojler":
        v, x = euler_method(differential_equation, v0, x0, t, m, k, mu_s, mu_d, g, wall)
    elif selected_m == "ModOjler":
        v,x = modified_euler_method(differential_equation, v0, x0, t, m, k, mu_s, mu_d, g, wall)
    elif selected_m == "RK4a":
        v, x = runge_kutta_alt_method(differential_equation, v0, x0, t, m, k, mu_s, mu_d, g, wall)
    else:
        v,x = modified_euler_alt(differential_equation, v0, x0, t, m, k, mu_s, mu_d, g, wall)
        

    selected = selected_option.get()
    #if selected == "Brzina":
        # Plotujemo
     #   ax.plot(t, v, label = 'Brzina')
      #  ax.set_xlabel('Vreme (s)')
       # ax.set_ylabel('Brzina')
        #ax.set_title('Kretanje tela sa oprugom i trenjem')
        #canvas.draw()
    if selected == "Energija":
        kinetic_energy = 0.5 * m * v**2
        potential_energy = 0.5 * k * x**2
        total_energy = kinetic_energy + potential_energy
        ax.plot(t, total_energy, label='Ukupna energija')
        ax.set_xlabel('Vreme (s)')
        ax.set_ylabel('Energija')
        ax.set_title('Energija tela s oprugom i trenjem')
        canvas.draw()
    elif selected == "Pozicija":
        ax.plot(t, x, label='Pozicija (deformacija opruge)')
        ax.set_xlabel('Vreme (s)')
        ax.set_ylabel('Pozicija')
        ax.set_title('Kretanje tela s oprugom i trenjem')
        canvas.draw()

    

def on_button_click():
    v = float(entry_parameter1.get())
    h = float(entry_parameter2.get())
    m = float(entry_parameter3.get())
    k = float(entry_parameter4.get())
    m_d = float(entry_parameter5.get())
    m_s = float(entry_parameter6.get())
    plot_graph(v, h, m, k,m_d,m_s)
    
def on_button_click_clear():
    ax.clear()
    canvas.draw()


# root prozor
root = tk.Tk()
root.geometry("1200x800")
#root.configure(bg="lightgray")
root.title("Spring Deformation Differential Modeling")

#cisto kozmeticki da ih razdvojim
left_frame = tk.Frame(root)
left_frame.pack(side=tk.LEFT, padx=10)
right_frame = tk.Frame(root)
right_frame.pack(side=tk.RIGHT, padx=10)
button_frame = tk.Frame(root)
button_frame.pack(side = tk.BOTTOM, pady=10)



label_parameter1 = Label(left_frame, text="Pocetna brzina:", font=14)
label_parameter1.pack(pady=5)
entry_parameter1 = Entry(left_frame)
entry_parameter1.insert(0,"100.0")
entry_parameter1.pack(pady=5)

label_parameter2 = Label(left_frame, text="Vremenski opseg:", font=14)
label_parameter2.pack(pady=5)
entry_parameter2 = Entry(left_frame)
entry_parameter2.insert(0, "0.5")
entry_parameter2.pack(pady=5)

label_parameter3 = Label(left_frame, text="Masa tela:", font=14)
label_parameter3.pack(pady=5)
entry_parameter3 = Entry(left_frame)
entry_parameter3.insert(0, "20.0")
entry_parameter3.pack(pady=5)

label_parameter4 = Label(left_frame, text="Koeficijent opruge:", font=14)
label_parameter4.pack(pady=5)
entry_parameter4 = Entry(left_frame)
entry_parameter4.insert(0, "0.1")
entry_parameter4.pack(pady=5)

label_parameter5 = Label(left_frame, text="Koef. Dinamickog trenja:", font=14)
label_parameter5.pack(pady=5)
entry_parameter5 = Entry(left_frame)
entry_parameter5.insert(0, "0.1")
entry_parameter5.pack(pady=5)

label_parameter6 = Label(left_frame, text="Koef. Statickog trenja:", font=14)
label_parameter6.pack(pady=5)
entry_parameter6 = Entry(left_frame)
entry_parameter6.insert(0, "0.5")
entry_parameter6.pack(pady=5)
separator = ttk.Separator(right_frame, orient="horizontal")

separator.pack(padx=10, pady=20, fill="x")


label_parameter5 = Label(right_frame, text="Izaberi numericku metodu:")
label_parameter5.configure(font=(17))
label_parameter5.pack(pady=5)
selected_method = tk.StringVar()
#vrsta metode
radio_button5 = tk.Radiobutton(right_frame, text="RK4 metoda", variable=selected_method, value="RK4")
radio_button5.pack(anchor="w",pady=5)
radio_button5.select()
radio_button6 = tk.Radiobutton(right_frame, text="Ojlerova metoda", variable=selected_method, value="Ojler")
radio_button6.pack(anchor="w",pady=5)
radio_button7 = tk.Radiobutton(right_frame, text="Modifikovana Ojlerova metoda", variable=selected_method, value="ModOjler")
radio_button7.pack(anchor="w",pady=5)
radio_button8 = tk.Radiobutton(right_frame, text="RK4 Alternativna Metoda", variable=selected_method, value="RK4a")
radio_button8.pack(anchor="w",pady=5)
radio_button9 = tk.Radiobutton(right_frame, text="Alternativna Modifikovana Ojlerova Metoda", variable=selected_method, value="AltOj")
radio_button9.pack(anchor="w",pady=5)

separator = ttk.Separator(right_frame, orient="horizontal")
separator.pack(padx=10, pady=30, fill="x")

label_parameter6 = Label(right_frame, text="Izaberi grafik:")
label_parameter6.configure(font=(17))
label_parameter6.pack(pady=15)
selected_option = tk.StringVar()
#vrsta plotovanja
#radio_button1 = tk.Radiobutton(right_frame, text="Brzina", variable=selected_option, value="Brzina")
#radio_button1.pack(pady=5)

radio_button2 = tk.Radiobutton(right_frame, text="Energija", variable=selected_option, value="Energija")
radio_button2.pack(pady=5)
radio_button2.select()
radio_button3 = tk.Radiobutton(right_frame, text="Pozicija", variable=selected_option, value="Pozicija")
radio_button3.pack(pady=5)

separator = ttk.Separator(right_frame, orient="horizontal")
separator.pack(padx=10, pady=20, fill="x")

label_parameter7 = Label(right_frame, text="Zid kao ogranicenje:")
label_parameter7.configure(font=(17))
label_parameter7.pack(pady=15)
selected_wall = tk.StringVar()
radio_button_zid_da = tk.Radiobutton(right_frame, text="Da", variable=selected_wall, value="da")
radio_button_zid_da.pack(pady=5)
radio_button_zid_da.select()
radio_button_zid_ne = tk.Radiobutton(right_frame, text="Ne", variable=selected_wall, value="ne")
radio_button_zid_ne.pack(pady=5)

separator = ttk.Separator(right_frame, orient="horizontal")
separator.pack(padx=10, pady=20, fill="x")



# Create a figure and axes for plotting
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# buttoni
button_plot = Button(button_frame, text="Nacrtaj", command=on_button_click, width=20, height=3, font=(14), bg="green", fg="white")
button_plot.pack(side=tk.LEFT,pady=10)

button_clear = Button(button_frame, text="Ocisti", command = on_button_click_clear, width=20, height=3, bg="coral", fg="white", font=(14))
button_clear.pack(pady=10, padx=10)



root.mainloop()
