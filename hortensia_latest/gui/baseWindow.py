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
        self.geometry("1200x900+50+50")
        self.resizable(False, False)
        self.configure(bg="white")
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

        ### Left box
        tk.Frame(self, bg=colorbg, width=280, height=575).place(x=15, y=15)

        ### Lower box with slider
        def get_current_time():
            return "%.1f"%(tReturn.get()*0.2)
        def get_current_step():
            return "%i"%tReturn.get()
        def slider_changed(event):
            timelabel.configure(text=get_current_time())
            steplabel.configure(text=get_current_step())

        tk.Frame(self, bg=colorbg, width=1170, height=275).place(x=15,y=610)

        tReturn = tk.DoubleVar()
        tslide = tk.Scale(self,
            from_=0, to=15000, bg=colorbg, troughcolor=colorfg,
            showvalue=0, orient='horizontal',
            command=slider_changed, variable=tReturn)
        tslide.place(x=30,y=630,width=1140)

        tk.Label(self, 
            text="Time Step:", bg=colorbg, anchor='w',font=("Helvetica",12)
            ).place(x=565,y=655,width=150,height=20)
        steplabel = tk.Label(self,
            bg=colorbg, anchor='e', font=("Helvetica",12), 
            text=get_current_step(),)
        steplabel.place(x=655,y=655,width=60,height=20)

        tk.Label(self,
            text="Time:               fs", bg=colorbg, anchor='w',
            font=("Helvetica",12)).place(x=725,y=655,width=150,height=20)
        timelabel = tk.Label(self,
            bg=colorbg, anchor='e', font=("Helvetica",12),
            text=get_current_time())
        timelabel.place(x=775,y=655,width=65,height=20)

        # Quit button
        quitbutton = tk.Button(self, borderwidth=0, highlightthickness=0,
            background='white', command=self.confirm)
        quitbutton['image'] = self.images[2]
        quitbutton.place(x=30, y=870, width=70, height=70, anchor="sw")

        ### Main box
        self.mainframe = tk.Frame(self, bg=colorbg, width=875, height=575)
        self.mainframe.place(x=310, y=15)
        px = 1 / plt.rcParams['figure.dpi']
        self.figure = Figure(figsize=(865 * px, 565 * px), dpi=100)
        self.fica = FigureCanvasTkAgg(self.figure, master=self)

        frame = tk.Frame(self, bg=colorbg, width=875, height=30)
        frame.place(x=315, y=585, anchor="sw")
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

        self.fica.get_tk_widget().place(x=315, y=20)
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