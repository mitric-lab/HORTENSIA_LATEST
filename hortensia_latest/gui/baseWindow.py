#!/usr/bin/env python3

import tkinter as tk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

colorbg = "#f5f5f5"
colorab = '#fcfcfc'
colorfg = "#dbdbdb"
coloruw = "#063e79"


class BaseWin(tk.Toplevel):
    def __init__(self, master, images):
        super().__init__(master)
        self.images = images

        self.protocol("WM_DELETE_WINDOW", self.on_close)
        self.title(" HORTENSIA Analysis Tool ")
        self.geometry("900x900+450+50")
        self.resizable(False, False)
        self.configure(bg="white")
        self.iconphoto(False, self.images[-1])
        self.option_add('*Dialog.msg.font', 'Helvetica 10')

    def on_close(self):
        self.destroy()
        self.master.quit()


class AnalysisWin(BaseWin):
    def __init__(self, master, images):
        super().__init__(master, images)
        self.logo = self.images[0].zoom(2,2)
        self.make_basic_appearance()

    def make_basic_appearance(self):
        # Background Würzburg logo
        bg = tk.Label(self, image=self.logo, bg='white')
        bg.place(x=0, y=0, width=1200, height=900)

        # Quit button
        #tk.Frame(self, bg=colorbg, width=100, height=100
        #    ).place(x=15, y=885, anchor='sw')
        quitbutton = tk.Button(self, borderwidth=0, highlightthickness=0,
            background='white', command=self.confirm)
        quitbutton['image'] = self.images[2]
        quitbutton.place(x=15, y=885, width=70, height=70, anchor="sw")


        ### Main box
        self.mainframe = tk.Frame(self, bg=colorbg, width=870, height=580)
        self.mainframe.place(x=15, y=60)
        px = 1 / plt.rcParams['figure.dpi']
        self.figure = Figure(figsize=(850 * px, 560 * px), dpi=100)
        self.fica = FigureCanvasTkAgg(self.figure, master=self)

        frame = tk.Frame(self, bg=colorbg, width=850, height=30)
        frame.place(x=25, y=595, anchor="nw")
        toolbar = NavigationToolbar2Tk(self.fica, frame)
        toolbar.config(bg='white')
        [i.config(bg='white') for i in toolbar.winfo_children()]
        toolbar.update()
        
        self.ax1 = self.figure.add_subplot()
        self.ax1.set_xlabel("", fontsize=18.0)
        self.ax1.set_ylabel("", fontsize=18.0)
        self.ax1.set_xlim([0.0, 1.0])
        self.ax1.tick_params(labelsize=15.0)
        self.ax1.set_ylim([0.0, 1.0])
        
        self.ax1.set_position([0.11, 0.11, 0.85, 0.85])

        self.fica.get_tk_widget().place(x=25, y=70)
        self.fica.draw()

    def confirm(self):
        def confirmclick(event):
            if event == 0:
                self.quit()
            elif event == 1:
                answer.destroy()

        answer = BaseWin(self, self.images)
        answer.title("Quitting Program")
        answer.geometry("300x120+500+440")
        tk.Label(answer, image=self.images[1], bg='white'
            ).place(x=0, y=0, width=300, height=120)

        tk.Label(answer, text='Do you really want to quit?',
            bg=colorbg, fg=coloruw, font=("Helvetica",12)
            ).place(x=15, y=15, width=270, height=50)
        tk.Frame(answer, bg=coloruw, width=270, height=2).place(x=15, y=63)

        tk.Button(answer,
            text="Quit program", bg="white",
            command=lambda: confirmclick(0)).place(
                x=135, y=80, width=85, height=30, anchor='ne')
        tk.Button(answer,
            text="Cancel", bg="white",
            command=lambda: confirmclick(1)).place(
                x=165, y=80, width=85, height=30, anchor='nw')