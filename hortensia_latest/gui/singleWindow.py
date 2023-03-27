#!/usr/bin/env python3

import os
from configparser import ConfigParser

import numpy as np
from scipy.signal import convolve2d

import tkinter as tk

from matplotlib import cm

import hortensia_latest.misc as misc
from hortensia_latest.gui.baseWindow import AnalysisWin

config = ConfigParser()
config.read('config.ini')

colorbg = "#f5f5f5"
colorab = '#fcfcfc'
colorfg = "#dbdbdb"
coloruw = "#063e79"


class OneWin(AnalysisWin):
    def __init__(self, master, images):
        super().__init__(master, images)
        self.title(" HORTENSIA Single Trajectory Analysis ")

        self.read_config()
        self.make_widgets()

    def read_config(self):
        self.tEns = float(config['Dynamics']['steps']) * \
                    float(config['Dynamics']['timestep'])
        hopdata = np.loadtxt("trajpop.dat")
        self.tTrj = hopdata[-1, 0]
        self.nTraj  = int(hopdata[0, 1])
        self.hopped = int(hopdata[0, 1] - hopdata[-1, 1])
        self.maxEk  = float(config['Continuum']['maxEk'])
        if os.path.exists("freeEl.dat"):
            self.freeEl = np.loadtxt("freeEl.dat")
            self.freeEl[:, 4] *= misc.hartree_to_eV
        else:
            self.freeEl = None

    def make_widgets(self):
        def change_tmax(t):
            self.tmax.set(t)
        def changedConv():
            if self.conv.get():
                self.wt['state'] = 'normal'
                self.we['state'] = 'normal'
            else:
                self.wt['state'] = 'disabled'
                self.we['state'] = 'disabled'

        #tk.Label(self, text="Single Trajectory Analysis", bg=colorbg, 
        #    font="Helvetica 15 bold", fg=coloruw
        #    ).place(x=450, y=15, width=870, height=30, anchor="n")
        #tk.Frame(self, bg=coloruw, width=870, height=2
        #    ).place(x=450, y=43, anchor='n')

        # Option Frames
        tk.Frame(self, bg=colorbg, width=206, height=145
            ).place(x=118, y=655, anchor='n')
        tk.Frame(self, bg=colorbg, width=206, height=145
            ).place(x=339, y=655, anchor='n')
        tk.Frame(self, bg=colorbg, width=206, height=230
            ).place(x=561, y=655, anchor='n')
        tk.Frame(self, bg=colorbg, width=206, height=230
            ).place(x=782, y=655, anchor='n')

        tk.Frame(self, bg=colorbg, width=342, height=70
            ).place(x=100, y=885, anchor='sw')

        # Time management
        tk.Label(self, text="Time plot range", font="Helvetica 12", bg=coloruw,
            fg='white').place(x=118, y=655, width=206, height=30, anchor="n")
        tk.Button(self, text="Max. trajectory time", bg="white", 
            font="Helvetica 9", command=lambda: change_tmax(self.tTrj)
            ).place(x=118, y=695, width=180, height=30, anchor="n")
        tk.Button(self, text="Max. ensemble time", bg="white", 
            font="Helvetica 9", command=lambda: change_tmax(self.tEns)
            ).place(x=118, y=730, width=180, height=30, anchor="n")
        textstr = "From 0 to\t              fs"
        tk.Label(self, text=textstr, anchor="center", font="Helvetica 9", 
            bg=colorbg,
            ).place(x=118, y=769, width=140, height=20, anchor="n")
        self.tmax = tk.DoubleVar(value=self.tEns)
        tk.Entry(self, bg="white", textvariable=self.tmax, justify="center", 
            borderwidth=0, highlightthickness=0, font="Helvetica 9"
            ).place(x=139, y=769, width=50, height=20, anchor="n")

        # Buttons for Histograms
        tk.Label(self, text="Histograms", font="Helvetica 12", bg=coloruw,
            fg='white').place(x=339, y=655, width=206, height=30, anchor="n")
        tk.Button(self, text="Hopping EKEs", bg="white",
            command=lambda: self.plotFreeEl("eke")
            ).place(x=339, y=695, width=180, height=30, anchor="n")
        tk.Button(self, text="Hopping VDEs", bg="white",
            command=lambda: self.plotFreeEl("vde")
            ).place(x=339, y=730, width=180, height=30, anchor="n")
        tk.Label(self, text="Number of bins", font="Helvetica 9", bg=colorbg, 
            anchor='w'
            ).place(x=274, y=769, width=90, height=20, anchor="nw")
        self.nBins = tk.IntVar(value=50)
        tk.Entry(self, bg="white", textvariable=self.nBins, justify="center", 
            borderwidth=0, highlightthickness=0
            ).place(x=404, y=769, width=40, height=20, anchor="ne")

        # Buttons for different files
        tk.Label(self, text="Graph Plots", font="Helvetica 12", bg=coloruw, 
            fg='white').place(x=561, y=655, width=206, height=30, anchor="n")
        tk.Button(self, text="trajpop.dat", bg="white", font="Helvetica 9",
            command=lambda: self.plotData("trajpop.dat", "population", 
                                          renorm=True)
            ).place(x=561, y=695, width=180, height=30, anchor="n")
        tk.Button(self, text="energies.dat", bg="white",font="Helvetica 9",
            command=lambda: self.plotData("energies.dat", "energies / eV", 
                                          customy=True, en=True)
            ).place(x=561, y=730, width=180, height=30, anchor="n")
        tk.Button(self, text="coefficients.dat", bg="white", font="Helvetica 9",
            command=lambda: self.plotData("coefficients.dat", "|c|$^2$ / a.u.")
            ).place(x=561, y=765, width=180, height=30, anchor="n")
        tk.Button(self, text="couplings.dat", bg="white", font="Helvetica 9",
            command=lambda: self.plotData("couplings.dat", "|V| / a.u.", 
            customy=True)
            ).place(x=561, y=800, width=180, height=30, anchor="n")

        # Buttons for 2D plots
        tk.Label(self, text="2D plots", font="Helvetica 12", bg=coloruw, 
            fg='white').place(x=782, y=655, width=206, height=30, anchor="n")
        tk.Button(self, text="Hopping EKEs", bg="white",
            command=lambda: self.plot2D("eke")
            ).place(x=782, y=695, width=180, height=30, anchor="n")
        tk.Button(self, text="Hopping VDEs", bg="white",
            command=lambda: self.plot2D("vde")
            ).place(x=782, y=730, width=180, height=30, anchor="n")

        tk.Label(self, text="Grid size", font="Helvetica 9", bg=colorbg, 
            anchor='w'
            ).place(x=702, y=765, width=103, height=30, anchor="nw")
        tk.Label(self, text="x", font="Helvetica 9", bg=colorbg
            ).place(x=812, y=765, width=100, height=30, anchor="n")
        self.n2D = [tk.IntVar(value=100), tk.IntVar(value=100)]
        tk.Entry(self, bg="white", textvariable=self.n2D[0], justify="center", 
            borderwidth=0, highlightthickness=0
            ).place(x=802, y=765, width=40, height=30, anchor="ne")
        tk.Entry(self, bg="white", textvariable=self.n2D[1], justify="center", 
            borderwidth=0, highlightthickness=0
            ).place(x=862, y=765, width=40, height=30, anchor="ne")
        
        self.conv = tk.BooleanVar(value=False)
        tk.Label(self, text="Gaussian\nconvolution", bg=colorbg,
            font="Helvetica 9"
            ).place(x=735, y=800, width=80, height=30, anchor="n")
        tk.Radiobutton(self, text='Yes', value=True, variable=self.conv,
            bg=colorbg, highlightthickness=0, activebackground=colorab,
            font="Helvetica 9", command=lambda: changedConv()
            ).place(x=820, y=800, width=50, height=30, anchor="ne")
        tk.Radiobutton(self, text='No', value=False, variable=self.conv, 
            bg=colorbg, highlightthickness=0, activebackground=colorab, 
            font="Helvetica 9", command=lambda: changedConv()
            ).place(x=820, y=800, width=50, height=30, anchor="nw")

        tk.Label(self, text="Standard\ndeviation", font="Helvetica 9", 
            bg=colorbg
            ).place(x=735, y=840, width=80, height=30, anchor="n")
        tk.Label(self, text="in t", font="Helvetica 9", bg=colorbg, anchor="w"
            ).place(x=815, y=830, width=60, height=30, anchor="n")
        tk.Label(self, text="in E", font="Helvetica 9", bg=colorbg, anchor="w"
            ).place(x=815, y=850, width=60, height=30, anchor="n")
        self.w2D = [tk.DoubleVar(value=10.0), tk.DoubleVar(value=0.01)]
        self.wt = tk.Entry(self, bg="white", textvariable=self.w2D[0], 
            justify="center", borderwidth=0, highlightthickness=0, 
            state="disabled")
        self.wt.place(x=815, y=833, width=43, height=20, anchor="nw")
        self.we = tk.Entry(self, bg="white", textvariable=self.w2D[1], 
            justify="center", borderwidth=0, highlightthickness=0, 
            state="disabled")
        self.we.place(x=815, y=855, width=43, height=20, anchor="nw")

    def plotData(self, infile, ylabel, renorm=False, customy=False, en=False):
        data = np.loadtxt(infile)
        t = data[:, 0]
        data = data[:, 1:]
        if renorm:
            data /= np.max(data)
        if en:
            data *= misc.hartree_to_eV

        self.ax1.clear()
        for i in range(len(data[0])):
            self.ax1.plot(t, data[:, i])
        self.ax1.set_xlabel("t / fs", fontsize=18.0)
        self.ax1.set_ylabel(ylabel, fontsize=18.0)
        self.ax1.set_xlim([0.0, self.tmax.get()])
        self.ax1.tick_params(labelsize=15.0)
        if customy:
            self.ax1.set_ylim([np.min(data) - abs(np.min(data))*0.1,
                               np.max(data) + abs(np.max(data))*0.1])
        else:
            self.ax1.set_ylim([0.0, 1.0])
        
        self.fica.draw()

    def plotFreeEl(self, mode):
        if self.freeEl is None:
            self.ax1.clear()
            self.ax1.text(0.5, 0.5, "No freeEl.dat file found!", fontsize=18.0, 
                transform=self.ax1.transAxes, ha='center', va='center')
            self.fica.draw()
            return

        self.ax1.clear()
        if mode == "vde":
            vde = self.freeEl[:, 5]
            (counts, bins) = np.histogram(vde, self.nBins.get())
            self.ax1.hist(bins[:-1], bins, weights=counts*100/self.hopped, 
                          ec="k")
            self.ax1.set_xlabel("VDE / eV")
            self.ax1.set_ylabel("hops / percent")
        elif mode == "eke":
            eke = self.freeEl[:, 4]
            (counts, bins) = np.histogram(eke, self.nBins.get(), 
                                          (0.0, self.maxEk))
            self.ax1.hist(bins[:-1], bins, range=(0.0, self.maxEk), 
                          weights=counts*100/self.hopped, ec="k")
            self.ax1.set_xlim(0.0, self.maxEk)
            self.ax1.set_xlabel("EKE / eV")
            self.ax1.set_ylabel("hops / percent")
        
        self.fica.draw()

    def plot2D(self, mode):
        if self.freeEl is None:
            self.ax1.clear()
            self.ax1.text(0.5, 0.5, "No freeEl.dat file found!", fontsize=18.0, 
                transform=self.ax1.transAxes, ha='center', va='center')
            self.fica.draw()
            return
        
        self.ax1.clear()
        if mode == "vde":
            vde   = self.freeEl[:, 5]

            dt = self.tmax.get()/self.n2D[0].get()
            de = (vde.max()-vde.min())/self.n2D[1].get()

            steps = self.freeEl[:, 0]
            steps = steps[steps < self.tmax.get()]

            t   = np.linspace(0, self.tmax.get(), self.n2D[0].get()+1)
            e   = np.linspace(vde.min(), vde.max(), self.n2D[1].get()+1)
            map = np.zeros((self.n2D[0].get()+1, self.n2D[1].get()+1))
            for i,step in enumerate(steps):
                it = int(step/dt)
                ie = int((vde[i]-vde.min())/de)
                map[it, ie] += 1

            if self.conv.get():
                gauss = np.exp(-0.5*(
                    (t[:,None]-self.tmax.get()/2)**2 / self.w2D[0].get()**2 + 
                    (e[None,:]-(vde.max()+vde.min())/2)**2 / 
                        self.w2D[1].get()**2))
                map = convolve2d(map, gauss, mode='same', boundary='symm')

            self.ax1.pcolormesh(t, e, map.T, cmap=cm.viridis)
            self.ax1.set_xlim(0.0, self.tmax.get())
            self.ax1.set_ylim(vde.min(), vde.max())
            self.ax1.set_xlabel("t / fs")
            self.ax1.set_ylabel("VDE / eV")

        elif mode == "eke":
            dt = self.tmax.get()/self.n2D[0].get()
            de = self.maxEk/self.n2D[1].get()

            steps = self.freeEl[:, 0]
            steps = steps[steps < self.tmax.get()]
            eke   = self.freeEl[:, 4]
            eke   = eke[eke < self.maxEk]

            t   = np.linspace(0, self.tmax.get(), self.n2D[0].get()+1)
            e   = np.linspace(0, self.maxEk, self.n2D[1].get()+1)
            map = np.zeros((self.n2D[0].get()+1, self.n2D[1].get()+1))
            for i,step in enumerate(steps):
                it = int(step/dt)
                ie = int(eke[i]/de)
                map[it, ie] += 1

            if self.conv.get():

                gauss = np.exp(-0.5*(
                    (t[:,None]-self.tmax.get()/2)**2 / self.w2D[0].get()**2 + 
                    (e[None,:]-self.maxEk/2)**2 / self.w2D[1].get()**2))
                map = convolve2d(map, gauss, mode='same', boundary='symm')

            self.ax1.pcolormesh(t, e, map.T, cmap=cm.viridis)
            self.ax1.set_xlim(0.0, self.tmax.get())
            self.ax1.set_ylim(0.0, self.maxEk)
            self.ax1.set_xlabel("t / fs")
            self.ax1.set_ylabel("EKE / eV")
        
        self.fica.draw()
