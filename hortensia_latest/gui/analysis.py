#!/usr/bin/env python3

import os

import tkinter as tk

import hortensia_latest.gui.infos as gI
from hortensia_latest.gui.baseWindow import BaseWin
from hortensia_latest.gui.singleWindow import OneWin
from hortensia_latest.gui.ensembleWindow import AllWin

colorbg = "#f5f5f5"
colorab = '#fcfcfc'
colorfg = "#dbdbdb"
coloruw = "#063e79"
#coloruw = "#f78c00"
#colorki = "#6e6d6d"

modpath = os.path.dirname(gI.__file__)
hortcwd = "%s/images/"%modpath

class AnalysisRoot(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('')
        self.geometry('600x450+50+50')
        self.resizable(False, False)
        self.withdraw()  # Hide default root Tk window.

        logo       = tk.PhotoImage(file=hortcwd+'logo_old.png')
        logoS      = tk.PhotoImage(file=hortcwd+'logo_old_small.png')
        quitButton = tk.PhotoImage(file=hortcwd+'quitbutton_blue.png')
        question   = tk.PhotoImage(file=hortcwd+'question_white.png')
        questionB  = tk.PhotoImage(file=hortcwd+'question_black.png')

        images = [logo, logoS, quitButton, question, questionB]

        TitleScreen(self.master, images)
        self.mainloop()


class TitleScreen(BaseWin):
    def __init__(self, master, images):
        super().__init__(master, images)
        self.geometry('400x300+450+350')
        self.iconphoto(False, tk.PhotoImage(file=hortcwd+'logo_new.png'))

        self.logo   = self.images[0].zoom(2,2)
        self.logo   = self.logo.subsample(3,3)

        self.make_widgets()

    def make_widgets(self):
        # Background Würzburg logo
        bg = tk.Label(self, image=self.logo, bg='white')
        bg.place(x=0, y=0, width=400, height=300)

        header = tk.Label(self,
            text="Do you want to analyse a \nsingle trajectory or an ensemble?", 
            fg=coloruw, bg=colorbg, font="Helvetica 12 bold")
        header.place(x=15, y=15, width=370, height=50)
        tk.Frame(self, bg=coloruw, width=370, height=2).place(x=15, y=63)

        # Button for single trajectory evaluation
        one = tk.Button(self, text="Single Trajectory Evaluation", bg="white",
                        font="Helvetica 12", command=self.make_one)
        one.place(x=200, y=96, width=250, height=40, anchor="n")
        stext  = "Opens a window for the visual representation of\n"
        stext += "output data of a single trajectory\n"
        stext += "One needs to be directly in the output folder"
        one.bind("<Enter>", lambda event: self.hover(stext))
        one.bind("<Leave>", lambda event: self.hover(""))

        # Button for ensemble evaluation
        all = tk.Button(self, text="Ensemble Evaluation", bg="white",
                        font="Helvetica 12", command=self.make_all)
        all.place(x=200, y=146, width=250, height=40, anchor="n")
        etext  = "Opens a window for the visual representation of\n"
        etext += "output data of the whole ensemble of trajectories\n"
        etext += "One needs to be in the folder containing the\n"
        etext += "TRAJ_X.out folders"
        all.bind("<Enter>", lambda event: self.hover(etext))
        all.bind("<Leave>", lambda event: self.hover(""))

        # quit button
        quitbutton = tk.Button(self, borderwidth=0, highlightthickness=0,
            background='white', command=self.confirm)
        quitbutton['image'] = self.images[2]
        quitbutton.place(x=15, y=215, width=70, height=70)
        qtext = "\nQuits the program, leaving everything unchanged"
        quitbutton.bind("<Enter>", lambda event: self.hover(qtext))
        quitbutton.bind("<Leave>", lambda event: self.hover(""))

        # tooltips frame
        tk.Frame(self, bg=colorbg).place(x=95, y=215, width=290, height=70)
        self.tooltip = tk.Label(self, font=("Helvetica", 9), text='', 
            bg=colorbg, fg=coloruw, anchor='nw', justify='left')
        self.tooltip.place(x=100, y=220, width=280, height=60)

    def hover(self, text):
        self.tooltip["text"] = text

    def confirm(self):
        def confirmclick(event):
            if event == 0:
                answer.destroy()

                saved = BaseWin(self, self.images)
                saved.title("Information")
                saved.geometry("300x120+500+440")
                tk.Label(saved, image=self.images[1], bg='white'
                    ).place(x=0, y=0, width=300, height=120)
                tk.Label(saved, text='config.ini\nwritten successfully!',
                    bg=colorbg, fg=coloruw, font=("Helvetica",12)
                    ).place(x=150, y=60, width=270, height=60, anchor='c')
                tk.Frame(saved, bg=coloruw, width=270, height=2
                    ).place(x=15, y=88)

                saved.after(1000, lambda: self.quit())

            elif event == 1:
                self.quit()
            elif event == 2:
                answer.destroy()

        answer = BaseWin(self, self.images)
        answer.title("Quiting Program")
        answer.geometry("300x120+500+440")
        tk.Label(answer, image=self.images[1], bg='white'
            ).place(x=0, y=0, width=300, height=120)

        tk.Label(answer, text='Do you want to quit, \n"+\
            "leaving everything unchanged?', bg=colorbg, fg=coloruw, 
            font=("Helvetica",12)
            ).place(x=15, y=15, width=270, height=50)
        tk.Frame(answer, bg=coloruw, width=270, height=2).place(x=15, y=63)

        tk.Button(answer,
            text="Quit program", bg="white",
            command=lambda: confirmclick(1)).place(
                x=135, y=80, width=85, height=30, anchor='ne')
        tk.Button(answer,
            text="Cancel", bg="white",
            command=lambda: confirmclick(2)).place(
                x=165, y=80, width=85, height=30, anchor='nw')

    def make_one(self):
        self.destroy()
        self.app = OneWin(self.master, self.images)

    def make_all(self):
        self.destroy()
        self.app = AllWin(self.master, self.images)


if __name__ == '__main__':
    AnalysisRoot()
