#!/usr/bin/env python3

from configparser import ConfigParser

import numpy as np

import tkinter as tk
from tkinter.messagebox import askokcancel, showinfo, WARNING

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

config = ConfigParser()
config.read('config.ini')

colorbg = "#f5f5f5"
colorfg = "#dbdbdb"

def clicked(x):
    print("yay! %i"%x)

def killed(event):
    quit("Killed program")

def confirm():
    answer = askokcancel(
        title='Warning',
        message='Are you sure you want to quit the program?',
        icon=WARNING)
    if answer:
        quit()

def clearFolder():
    tmp  = 'Are you sure you want to clear all folders?\n'
    tmp += 'This gets rid of all data that is not necessary '
    tmp += 'to start new trajectories with the HORTENSIA program'
    answer = askokcancel(
        title='Warning',
        message=tmp,
        icon=WARNING)
    if answer:
        showinfo(
            title='Deletion Status',
            message='Data deleted successfully')

def plotData(infile, customy0=False):
    tk.Frame(root, bg=colorbg, width=930, height=600).place(x=260, y=10)

    data = np.loadtxt(infile)
    t = data[:, 0]
    data = data[:, 1:]

    px = 1 / plt.rcParams['figure.dpi']
    figure = Figure(figsize=(926 * px, 596 * px), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, master=root)

    axes = figure.add_subplot()
    for i in range(len(data[0])):
        axes.plot(t, data[:, i])
    axes.set_xlabel("t / fs")
    axes.set_ylabel("some data")
    axes.set_xlim([min(t), max(t)])
    if customy0:
        axes.set_ylim([np.min(data) - abs(np.min(data))*0.1,
                       np.max(data) + abs(np.max(data))*0.1])
    else:
        axes.set_ylim([0.0, np.max(data) + abs(np.max(data))*0.1])
    figure.tight_layout()
    fica = figure_canvas.get_tk_widget()
    fica.place(x=262, y=12)

def plotFreeEl(mode):
    tk.Frame(root, bg=colorbg, width=930, height=600).place(x=260, y=10)

    data = np.loadtxt("freeEl.dat")
    t = data[:, 0]
    ki = data[:, 1:4]
    eke = data[:, 4]
    vde = data[:, 5]

    px = 1 / plt.rcParams['figure.dpi']
    figure = Figure(figsize=(926 * px, 596 * px), dpi=100)
    figure_canvas = FigureCanvasTkAgg(figure, master=root)

    axes = figure.add_subplot()
    if mode == "vde":
        axes.hist(vde, int(nBins.get()), density=True)
        axes.set_xlabel("VDE / eV")
        axes.set_ylabel("hops / percent")
    elif mode == "eke":
        axes.hist(eke*27.211386, int(nBins.get()), density=True)
        axes.set_xlim(0.0, float(config['Continuum']['maxEk']))
        axes.set_xlabel("EKE / eV")
        axes.set_ylabel("hops / percent")
    #for i in range(len(data[0])):
    #    axes.plot(t, data[:, i])
    #axes.set_ylabel("some data")
    #axes.set_xlim([min(t), max(t)])
    #if customy0:
    #    axes.set_ylim([np.min(data) - abs(np.min(data))*0.1,
    #                   np.max(data) + abs(np.max(data))*0.1])
    #else:
    #    axes.set_ylim([0.0, np.max(data) + abs(np.max(data))*0.1])
    figure.tight_layout()
    fica = figure_canvas.get_tk_widget()
    fica.place(x=262, y=12)


root = tk.Tk(className=" HORTENSIA GUI ")
root.geometry("1200x900+50+50")
root.resizable(False,False)
root.configure(bg='white')
root.option_add('*Dialog.msg.font', 'Helvetica 10')

root.bind('<Return>', killed)

### Left box
tk.Frame(root, bg=colorbg, width=240, height=880).place(x=10,y=10)

tk.Label(root,
    text="Ensemble Analysis",
    bg=colorbg,
    font=("Helvetica",12)).place(x=20,y=10,width=220,height=50)
tk.Button(root,
    text="Check Trajectories",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=60,width=220,height=30)
tk.Button(root,
    text="Evaluate Free Electrons",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=95,width=220,height=30)
tk.Button(root,
    text="Check Hopping Angles",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=130,width=220,height=30)
tk.Button(root,
    text="Check Hopping Bond Lengths",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=165,width=220,height=30)

tk.Label(root,
    text="Single Trajectory Analysis",
    bg=colorbg,
    font=("Helvetica",12)).place(x=20,y=195,width=220,height=50)
tk.Button(root,
    text="trajpop.dat",
    bg="white",
    command=lambda: plotData("trajpop.dat")).place(x=20,y=245,width=220,height=30)
tk.Button(root,
    text="energies.dat",
    bg="white",
    command=lambda: plotData("energies.dat", True)).place(x=20,y=280,width=220,height=30)
tk.Button(root,
    text="coefficients.dat",
    bg="white",
    command=lambda: plotData("coefficients.dat")).place(x=20,y=315,width=220,height=30)
tk.Button(root,
    text="couplings.dat",
    bg="white",
    command=lambda: plotData("couplings.dat")).place(x=20,y=350,width=220,height=30)

tk.Button(root,
    text="Hopping EKEs from freeEl.dat",
    bg="white",
    command=lambda: plotFreeEl("eke")).place(x=20,y=405,width=220,height=30)
tk.Button(root,
    text="Hopping VDEs from freeEl.dat",
    bg="white",
    command=lambda: plotFreeEl("vde")).place(x=20,y=440,width=220,height=30)
nBins = tk.StringVar(value=50)
tk.Entry(root,
    bg="white", textvariable=nBins, justify="center", borderwidth=0,
    highlightthickness=0).place(x=20,y=480,width=220,height=30)
"""
tk.Button(root,
    text="Option 7",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=455,width=220,height=30)
tk.Button(root,
    text="Option 8",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=490,width=220,height=30)
tk.Button(root,
    text="Option 9",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=525,width=220,height=30)
tk.Button(root,
    text="Option 10",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=560,width=220,height=30)
"""

tk.Label(root,
    text="Clear All Folders\nfrom Unnecessary Files",
    bg=colorbg,
    font=("Helvetica",12)).place(x=20,y=625,width=220,height=50)
tk.Button(root,
    text="Clear Folders",
    bg="white",
    command=clearFolder).place(x=20,y=680,width=220,height=30)

tk.Label(root,
    text="Pre-Dynamics Options",
    bg=colorbg,
    font=("Helvetica",12)).place(x=20,y=710,width=220,height=50)
tk.Button(root,
    text="Create Wigner Ensemble",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=760,width=220,height=30)
tk.Button(root,
    text="Decontract Basis",
    bg="white",
    command=lambda: clicked(1)).place(x=20,y=795,width=220,height=30)

tk.Button(root,
    text="Quit Program",
    bg="white",
    command=confirm).place(x=20,y=850,width=220,height=30)

tk.Frame(root, bg=colorfg, width=200, height=2).place(x=30,y=614)



### lower box with slider
def get_current_time():
    return "%.1f"%(tReturn.get()*0.2)
def get_current_step():
    return "%i"%tReturn.get()
def slider_changed(event):
    timelabel.configure(text=get_current_time())
    steplabel.configure(text=get_current_step())

tk.Frame(root, bg=colorbg, width=930, height=275).place(x=260,y=620)

tReturn = tk.DoubleVar()

tslide = tk.Scale(
    root,
    from_=0, to=15000,
    bg=colorbg, troughcolor=colorfg,
    showvalue=0,
    orient='horizontal',
    command=slider_changed,
    variable=tReturn)
tslide.place(x=270,y=630,width=910)

tk.Label(root,
    text="Time Step:",
    bg=colorbg,
    anchor='w',
    font=("Helvetica",12)).place(x=565,y=655,width=150,height=20)
steplabel = tk.Label(
    root,
    bg=colorbg,
    anchor='e',
    font=("Helvetica",12),
    text=get_current_step(),)
steplabel.place(x=655,y=655,width=60,height=20)

tk.Label(root,
    text="Time:               fs",
    bg=colorbg,
    anchor='w',
    font=("Helvetica",12)).place(x=725,y=655,width=150,height=20)
timelabel = tk.Label(
    root,
    bg=colorbg,
    anchor='e',
    font=("Helvetica",12),
    text=get_current_time())
timelabel.place(x=775,y=655,width=65,height=20)

tk.Frame(root, bg=colorfg, width=830, height=2).place(x=310,y=684)


### main box
px = 1 / plt.rcParams['figure.dpi']
figure = Figure(figsize=(926 * px, 596 * px), dpi=100)
figure_canvas = FigureCanvasTkAgg(figure, master=root)

fica = figure_canvas.get_tk_widget()
fica.place(x=262, y=12)

root.mainloop()
