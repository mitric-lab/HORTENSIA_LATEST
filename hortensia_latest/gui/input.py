#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

import tkinter as tk
from tkinter import ttk

import hortensia_latest.gui.infos as gI

colorbg = "#f5f5f5"
colorab = '#fcfcfc'
colorfg = "#dbdbdb"
coloruw = "#063e79"

modpath = os.path.dirname(gI.__file__)
hortensiacwd  = "%s/images/"%modpath

def makeFolders():
    if createFolders.get():
        trajs = nrTraj.get()
        for i in range(trajs):
            if not os.path.isdir("TRAJ_%i"%i):
                os.system("mkdir TRAJ_%i"%i)

def callWigner():
    if doWigner.get():
        if os.path.isdir("INITSTRUCT"):
            os.system("rm -f INITSTRUCT/structure*in")
        else:
            os.system("mkdir INITSTRUCT")
        os.system("cp %s %s INITSTRUCT"%(input1.get(), input2.get()))

        cwd = os.getcwd()
        os.chdir(cwd+"/INITSTRUCT")

        from hortensia_latest.wigner.WignerBase import WignerEnsemble
        qc = qcWigner.get()
        if qc == "g09":
            from hortensia_latest.wigner.HarmonicGaussian import \
                HarmonicGaussian as HarmonicInterface
        elif qc == "QChem":
            from hortensia_latest.wigner.HarmonicQChem import \
                HarmonicQChem as HarmonicInterface

        if distType.get() == "w":
            DistClass = WignerEnsemble(HarmonicInterface(linWigner.get(),
                input1.get(), input2.get()), float(modes.get()), nrInCon.get(),
                "structure", distType.get().upper(), [])
        else:
            if len(modes.get().split(",")) == 1:
                modesV1H = np.asarray([modes.get().split(",")], dtype=int)
            else:
                modesV1H = np.asarray(modes.get().split(","), dtype=int)
            DistClass = WignerEnsemble(HarmonicInterface(linWigner.get(),
                input1.get(), input2.get()), 0.1, nrInCon.get(), "structure",
                distType.get().upper(), modesV1H)
        DistClass.getWignerEnsemble()

        os.chdir(cwd)

def copyToFolders():
    if cpFiles.get():
        cpcommand = "cp config.ini "
        if os.path.exists("basis"):
            cpcommand += "basis "
        else:
            print("basis file missing!\nfile can't be copied to folders")

        if os.path.isdir("INITSTRUCT"):
            struct = True
        else:
            struct = False

        for i in range(nrTraj.get()):
            if os.path.isdir("TRAJ_%i"%i):
                os.system("%s TRAJ_%i"%(cpcommand, i))
                if struct:
                    cmd  = "cp INITSTRUCT/structure_%i.in "%i
                    cmd += "TRAJ_%i/structure.in"%i
                    os.system(cmd)

def confirm():
    def confirmclick(event):
        if event == 0:
            data = packData()

            import hortensia_latest.gui.writeInput as writeInput
            writeInput.write(data)
            if subScript.get():
                writeInput.writeSub(Nproc.get(), steps.get()*timestep.get(),
                                       qcmethod.get())

            answer.destroy()

            saved = tk.Toplevel(root)
            saved.title("Information")
            saved.geometry("300x120+200+215")
            saved.resizable(False,False)
            saved.configure(bg='white')
            saved.iconphoto(False, favicon)
            tk.Label(saved, image=logoS, bg='white'
                ).place(x=0, y=0, width=300, height=120)
            tk.Label(saved, text='config.ini\nwritten successfully!',
                bg=colorbg, fg=coloruw, font=("Helvetica",12)
                ).place(x=150, y=60, width=270, height=60, anchor='c')
            tk.Frame(saved, bg=coloruw, width=270, height=2).place(x=15, y=88)

            makeFolders()
            callWigner()
            copyToFolders()

            saved.after(1000, lambda: root.destroy())

        elif event == 1:
            quit()
        elif event == 2:
            answer.destroy()

    answer = tk.Toplevel(root)
    answer.title("Quiting Program")
    answer.geometry("300x120+200+215")
    answer.resizable(False, False)
    answer.configure(bg='white')
    answer.iconphoto(False, favicon)
    tk.Label(answer, image=logoS, bg='white'
        ).place(x=0, y=0, width=300, height=120)

    tk.Label(answer, text='Do you want\nto save your settings?',
        bg=colorbg, fg=coloruw, font=("Helvetica",12)
        ).place(x=15, y=15, width=270, height=50)
    tk.Frame(answer, bg=coloruw, width=270, height=2).place(x=15, y=63)

    tk.Button(answer,
        text="Save & Quit", bg="white", fg='black', activeforeground='black',
        command=lambda: confirmclick(0)).place(
            x=15, y=80, width=85, height=30)
    tk.Button(answer,
        text="Don't Save", bg="white", fg='black', activeforeground='black',
        command=lambda: confirmclick(1)).place(
            x=150, y=80, width=85, height=30, anchor='n')
    tk.Button(answer,
        text="Cancel", bg="white", fg='black', activeforeground='black',
        command=lambda: confirmclick(2)).place(
            x=285, y=80, width=85, height=30, anchor='ne')

def packData():
    data = []
    if len(gridConf) == 1:
        tmp0, tmp1, tmp2, tmp3, tmp4, tmp5 = '', '', '', '', '', ''
        tmp0 = "%s"%gridConf[0][0]
        if gridConf[0][0] != 'cubic':
            tmp1 = "%s"%gridConf[0][1]
            tmp2 = "%s"%gridConf[0][2]
            if gridConf[0][0] == 'fib':
                tmp3 = "%s"%gridConf[0][3]
        else:
            tmp4 = gridConf[0][1].split(',')
            tmp5 = gridConf[0][2].split(',')
        tmp = [tmp0, tmp1, tmp2, tmp3, tmp4, tmp5]
    else:
        tmp0, tmp1, tmp2, tmp3 = '[', '[', '[', '['
        tmp4, tmp5 = ['[','[','['], ['[','[','[']
        for i in range(len(gridConf)):
            tmp0 += "'%s',"%gridConf[i][0]

            if gridConf[i][0] == 'cubic':
                kx,ky,kz = gridConf[i][1].split(',')
                nx,ny,nz = gridConf[i][2].split(',')
                tmp4[0] += "%s,"%kx
                tmp4[1] += "%s,"%ky
                tmp4[2] += "%s,"%kz
                tmp5[0] += "%s,"%nx
                tmp5[1] += "%s,"%ny
                tmp5[2] += "%s,"%nz
            else:
                tmp1 += "%s,"%gridConf[i][1]
                tmp2 += "%s,"%gridConf[i][2]
                if gridConf[i][0] == 'fib':
                    tmp3 += "%s,"%gridConf[i][3]

        tmp = []
        for i in [tmp0,tmp1,tmp2,tmp3]:
            if i != '[':
                i = i[:-1] + ']'
                tmp.append(i)
            else:
                tmp.append('')
        if tmp4[0] != '[':
            tmp4[0] = tmp4[0][:-1] + ']'
            tmp4[1] = tmp4[1][:-1] + ']'
            tmp4[2] = tmp4[2][:-1] + ']'
            tmp5[0] = tmp5[0][:-1] + ']'
            tmp5[1] = tmp5[1][:-1] + ']'
            tmp5[2] = tmp5[2][:-1] + ']'
        else:
            tmp4 = ''
            tmp5 = ''

    data.append(tmp[0])
    data.append(tmp[1])
    data.append(tmp[2])
    data.append(tmp[3])
    data.append(tmp4)
    data.append(tmp5)
    data.append(exact2elInt.get())
    data.append(intNk.get())
    data.append(intkSkip.get())
    data.append(anion_states.get())
    data.append(starting_state.get())
    data.append(precision.get())
    data.append(max_micro.get())
    data.append(Eshift.get())
    data.append(orthotype.get())
    data.append(timestep.get())
    data.append(steps.get())
    data.append(popsteps.get())
    data.append(nrpertraj.get())
    data.append(restart.get())
    data.append(printlevel.get())
    data.append(skipstep.get())
    data.append([coupprint.get(), coefprint.get(), nadiaprint.get(),
                 probprint.get(), enerprint.get()])
    data.append(qcmethod.get())
    data.append(charge.get())
    data.append(mult.get())
    data.append(func.get())
    data.append(Nproc.get())
    data.append(convergence.get())
    data.append(maxiter.get())
    data.append(zpeInclude.get())
    data.append(zpeDiff.get())
    data.append(vibExEn.get())
    data.append(adiaDecay.get())

    return data

### root window
root = tk.Tk(className=" HORTENSIA Input Generator ")
root.geometry("600x450+50+50")
root.resizable(False,False)
root.configure(bg='white')
root.option_add('*Dialog.msg.font', 'Helvetica 10')

### Images
logo       = tk.PhotoImage(file=hortensiacwd+'logo_normal.png')
logoS      = tk.PhotoImage(file=hortensiacwd+'logo_small.png')
favicon    = tk.PhotoImage(file=hortensiacwd+'favicon.png')
quitButton = tk.PhotoImage(file=hortensiacwd+'quitbutton.png')
question   = tk.PhotoImage(file=hortensiacwd+'question_white.png')
questionB  = tk.PhotoImage(file=hortensiacwd+'question_black.png')

### Set favicon
root.iconphoto(False, favicon)

### create a notebook
notebook = ttk.Notebook(root)
notebook.place(width=600, height=450)
s = ttk.Style()
s.theme_use('default')
s.configure('TNotebook', background=colorbg)
s.configure('TNotebook.Tab', background=colorbg, foreground='black')
s.map("TNotebook.Tab", background=[("selected", coloruw)],
                       foreground=[("selected", 'white')])

### create frames
frame1 = tk.Frame(notebook, width=600, height=450, bg='white')
frame2 = tk.Frame(notebook, width=600, height=450, bg='white')
frame3 = tk.Frame(notebook, width=600, height=450, bg='white')
frame4 = tk.Frame(notebook, width=600, height=450, bg='white')
frame5 = tk.Frame(notebook, width=600, height=450, bg='white')
frame6 = tk.Frame(notebook, width=600, height=450, bg='white')
frame7 = tk.Frame(notebook, width=600, height=450, bg='white')
frame1.place()
frame2.place()
frame3.place()
frame4.place()
frame5.place()
frame6.place()
frame7.place()

### add frames to notebook
notebook.add(frame1, text='Dynamics')
notebook.add(frame2, text='Quantum Chemistry')
notebook.add(frame3, text='States')
notebook.add(frame4, text='2e-Integrals')
notebook.add(frame5, text='Free Electrons')
notebook.add(frame6, text='Folders')
notebook.add(frame7, text='Wigner')

### For testing
#def killed(event):
#    quit("Killed program")
#root.bind('<Return>', killed)
notebook.select(frame7) #!!!


################################################################################
#                                  Frame 1                                     #
################################################################################
tk.Label(frame1, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover1(event, text):
    tooltip1["text"] = text

def disableSkip(event, state):
    if state.get() == -1:
        skipPrint['state'] = 'normal'
    else:
        skipPrint['state'] = 'disabled'

opt  = ["full", "summed", "none"]

# headline
tk.Label(frame1,
    text="General options", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=20, y=10, width=270, height=35)
tk.Frame(frame1, bg=coloruw, width=270, height=2).place(x=20, y=43)

# nuclear steps
steps = tk.IntVar(value=5000)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=20, y=60)
tk.Label(frame1, text="Number of nuclear dynamics steps", bg=coloruw,
    fg='white', font=("Helvetica", 10)).place(x=20, y=60, width=270, height=20)
tk.Entry(frame1, bg=colorbg, fg='black', textvariable=steps,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=80, width=270, height=25)

# timestep
timestep = tk.DoubleVar(value=0.2)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=20, y=115)
tk.Label(frame1, text="Nuclear dynamics time step (fs)", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=115, width=270, height=20)
tk.Entry(frame1, bg=colorbg, fg='black', textvariable=timestep,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=135, width=270, height=25)

# number of popsteps per nuclear step
popsteps = tk.IntVar(value=100)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=20, y=170)
tk.Label(frame1, text="Number of integration steps", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=170, width=270, height=20)
tk.Entry(frame1, bg=colorbg, fg='black', textvariable=popsteps,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=190, width=270, height=25)
# tooltip
q0 = tk.Label(frame1, image=question, bg=coloruw)
q0.place(x=245, y=175, width=10, height=10)
q0.bind("<Enter>", lambda event, arg=gI.q0: hover1(event, arg))
q0.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))


# number of subtrajs per traj
nrpertraj = tk.IntVar(value=100)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=20, y=225)
tk.Label(frame1, text="Number of subtrajectories", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=225, width=270, height=20)
tk.Entry(frame1, bg=colorbg, fg='black', textvariable=nrpertraj,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=245, width=270, height=25)
# tooltip
q1 = tk.Label(frame1, image=question, bg=coloruw)
q1.place(x=237, y=230, width=10, height=10)
q1.bind("<Enter>", lambda event, arg=gI.q1: hover1(event, arg))
q1.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

# restart option
restart = tk.BooleanVar(value=False)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=20, y=280)
tk.Label(frame1, text="Trajectory restart", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=280, width=270, height=20)
tk.Radiobutton(frame1, text='New Trajectory', value=False, variable=restart,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=20, y=300, width=135, height=25)
tk.Radiobutton(frame1, text='Restart', value=True, variable=restart, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=155, y=300, width=135, height=25)

# headline
tk.Label(frame1, text="Printing options", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=310, y=10, width=270, height=35)
tk.Frame(frame1, bg=coloruw, width=270, height=2).place(x=310, y=43)

# name of output file
output = tk.StringVar(value='out.out')
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=310, y=60)
tk.Label(frame1, text="Output file name", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=60, width=270, height=20)
tk.Entry(frame1, bg=colorbg, fg='black', textvariable=output,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=80, width=270, height=25)

# couplings print options
coupprint  = tk.StringVar(value='full')
nadiaprint = tk.StringVar(value='full')
coefprint  = tk.StringVar(value='full')
probprint  = tk.StringVar(value='full')
enerprint  = tk.StringVar(value='full')

tk.Frame(frame1, bg=colorbg, width=270, height=155).place(x=310, y=60)
tk.Label(frame1, text="Print level of calculated variables", bg=coloruw,
    fg='white', font=("Helvetica", 10)).place(x=310, y=60, width=270, height=20)

tk.Label(frame1, text="total couplings", bg=colorbg, fg='black',
    font=("Helvetica", 9)).place(x=310, y=80, width=180, height=25)
coupOpt = tk.OptionMenu(frame1, coupprint, *['full', 'summed', 'none'])
coupOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
coupOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
coupOpt.place(x=580, y=80, width=90, height=25, anchor='ne')
# tooltip
q2 = tk.Label(frame1, image=questionB, bg=colorbg)
q2.place(x=445, y=87, width=10, height=10)
q2.bind("<Enter>", lambda event, arg=gI.q2: hover1(event, arg))
q2.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

tk.Label(frame1, text="non-ad./diab. couplings", bg=colorbg, fg='black',
    font=("Helvetica", 9)).place(x=310, y=107, width=180, height=25)
nadiOpt = tk.OptionMenu(frame1, nadiaprint, *['full', 'summed', 'none'])
nadiOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
nadiOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
nadiOpt.place(x=580, y=107, width=90, height=25, anchor='ne')
# tooltip
q3 = tk.Label(frame1, image=questionB, bg=colorbg)
q3.place(x=471, y=114, width=10, height=10)
q3.bind("<Enter>", lambda event, arg=gI.q2: hover1(event, arg))
q3.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

tk.Label(frame1, text="coefficients", bg=colorbg, fg='black',
    font=("Helvetica", 9)).place(x=310, y=134, width=180, height=25)
coefOpt = tk.OptionMenu(frame1, coefprint, *['full', 'summed', 'none'])
coefOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
coefOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
coefOpt.place(x=580, y=134, width=90, height=25, anchor='ne')
# tooltip
q4 = tk.Label(frame1, image=questionB, bg=colorbg)
q4.place(x=435, y=141, width=10, height=10)
q4.bind("<Enter>", lambda event, arg=gI.q4: hover1(event, arg))
q4.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

tk.Label(frame1, text="probabilities", bg=colorbg, fg='black',
    font=("Helvetica", 9)).place(x=310, y=161, width=180, height=25)
probOpt = tk.OptionMenu(frame1, probprint, *['full', 'summed', 'none'])
probOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
probOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
probOpt.place(x=580, y=161, width=90, height=25, anchor='ne')
# tooltip
q5 = tk.Label(frame1, image=questionB, bg=colorbg)
q5.place(x=437, y=168, width=10, height=10)
q5.bind("<Enter>", lambda event, arg=gI.q5: hover1(event, arg))
q5.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

tk.Label(frame1, text="energies", bg=colorbg, fg='black', font=("Helvetica", 9)
    ).place(x=310, y=188, width=180, height=25)
enerOpt = tk.OptionMenu(frame1, enerprint, *['full', 'none'])
enerOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
enerOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
enerOpt.place(x=580, y=188, width=90, height=25, anchor='ne')
# tooltip
q6 = tk.Label(frame1, image=questionB, bg=colorbg)
q6.place(x=428, y=195, width=10, height=10)
q6.bind("<Enter>", lambda event, arg=gI.q6: hover1(event, arg))
q6.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

# printlevel option
printlevel = tk.IntVar(value=0)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=310, y=225)
tk.Label(frame1, text="Print level", bg=coloruw, fg='white',
    font=("Helvetica",10)).place(x=310, y=225, width=270, height=20)
optMenu1 = tk.OptionMenu(frame1, printlevel, *[-1, 0, 1],
    command=lambda event, arg=printlevel: disableSkip(event, arg))
optMenu1.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
optMenu1['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
optMenu1.place(x=310, y=245, width=270, height=25)
# tooltip
q7 = tk.Label(frame1, image=question, bg=coloruw)
q7.place(x=481, y=230, width=10, height=10)
q7.bind("<Enter>", lambda event, arg=gI.q7: hover1(event, arg))
q7.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))

# number of steps to skip for print
skipstep = tk.IntVar(value=0)
tk.Frame(frame1, bg=colorbg, width=270, height=45).place(x=310, y=280)
tk.Label(frame1, text="Skipped steps", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=280, width=270, height=20)
skipPrint = tk.Entry(frame1, bg=colorbg, fg='black', textvariable=skipstep,
    font=("Helvetica", 9), justify='center', borderwidth=0,
    highlightthickness=0, state='disabled')
skipPrint.place(x=310, y=300, width=270, height=25)

# tooltips frame
tk.Frame(frame1, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip1 = tk.Label(frame1, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip1.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton1 = tk.Button(frame1, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton1['image'] = quitButton
quitbutton1.place(x=20, y=340, width=70, height=70)
quitbutton1.bind("<Enter>", lambda event, arg=gI.qb: hover1(event, arg))
quitbutton1.bind("<Leave>", lambda event, arg=""   : hover1(event, arg))


################################################################################
#                                  Frame 2                                     #
################################################################################
tk.Label(frame2, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

# headline
tk.Label(frame2, text="Options for molecular properties", fg=coloruw,
    bg=colorbg, font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame2, bg=coloruw, width=270, height=2).place(x=165, y=43)

# molecular charge
charge = tk.IntVar(value=-1)
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=20, y=60)
tk.Label(frame2, text="Molecular charge", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=60, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=charge,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=80, width=270, height=25)

# multiplicity
mult = tk.IntVar(value=2)
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=20, y=115)
tk.Label(frame2, text="Multiplicity of system", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=115, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=mult,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=135, width=270, height=25)

# QC Method
qcmethod = tk.StringVar(value='g09')
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=310, y=60)
tk.Label(frame2, text="External quantum chemistry program", bg=coloruw,
    fg='white', font=("Helvetica", 10)).place(x=310, y=60, width=270, height=20)
tk.Radiobutton(frame2, text='Gaussian09', value='g09', variable=qcmethod,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=310, y=80, width=90, height=25)
tk.Radiobutton(frame2, text='Gaussian16', value='g16', variable=qcmethod,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=400, y=80, width=90, height=25)
tk.Radiobutton(frame2, text='QChem', value='qchem', variable=qcmethod,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=490, y=80, width=90, height=25)

# functional
func = tk.StringVar(value='wB97XD')
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=310, y=115)
tk.Label(frame2, text="Functional for DFT calculation", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=115, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=func,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=135, width=270, height=25)

# Number of threads
Nproc = tk.IntVar(value=2)
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=310, y=170)
tk.Label(frame2, text="Number of threads for DFT calculation", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=310, y=170, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=Nproc,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=190, width=270, height=25)

# SCF convergence
convergence = tk.IntVar(value=8)
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=310, y=225)
tk.Label(frame2, text="SCF convergence criterion (10\u207B\u207F)", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=310, y=225, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=convergence,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=245, width=270, height=25)

# maximum number of SCF iterations
maxiter = tk.IntVar(value=100)
tk.Frame(frame2, bg=colorbg, width=270, height=45).place(x=310, y=280)
tk.Label(frame2, text="Maximum number of SCF iterations", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=310, y=280, width=270, height=20)
tk.Entry(frame2, bg=colorbg, fg='black', textvariable=maxiter,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=300, width=270, height=25)

# tooltips frame
tk.Frame(frame2, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip2 = tk.Label(frame2, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip2.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton2 = tk.Button(frame2, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton2['image'] = quitButton
quitbutton2.place(x=20, y=340, width=70, height=70)


################################################################################
#                                  Frame 3                                     #
################################################################################
tk.Label(frame3, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover3(event, text):
    tooltip3["text"] = text

# headline
tk.Label(frame3,
    text="Options for state configuration", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame3, bg=coloruw, width=270, height=2).place(x=165, y=43)

# anion States
anion_states  = tk.StringVar(value='0')
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=20, y=60)
tk.Label(frame3, text="Anion states (comma separated)", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=60, width=270, height=20)
tk.Entry(frame3, bg=colorbg, fg='black', textvariable=anion_states,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=80, width=270, height=25)
# tooltip
q8 = tk.Label(frame3, image=question, bg=coloruw)
q8.place(x=261, y=65, width=10, height=10)
q8.bind("<Enter>", lambda event, arg=gI.q8: hover3(event, arg))
q8.bind("<Leave>", lambda event, arg=""   : hover3(event, arg))

# starting state
starting_state = tk.IntVar(value=0)
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=310, y=60)
tk.Label(frame3, text="Starting state (index of state list)", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=310, y=60, width=270, height=20)
tk.Entry(frame3, bg=colorbg, fg='black', textvariable=starting_state,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=80, width=270, height=25)

# energy shift
Eshift = tk.DoubleVar(value=0.0)
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=165, y=115)
tk.Label(frame3, text="Energy shift anion/neutral (in eV)", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=165, y=115, width=270, height=20)
tk.Entry(frame3, bg=colorbg, fg='black', textvariable=Eshift,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=165, y=135, width=270, height=25)
# tooltip
q11 = tk.Label(frame3, image=question, bg=coloruw)
q11.place(x=408, y=120, width=10, height=10)
q11.bind("<Enter>", lambda event, arg=gI.q11: hover3(event, arg))
q11.bind("<Leave>", lambda event, arg=""    : hover3(event, arg))

# precision
precision = tk.DoubleVar(value=99.5)
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=20, y=170)
tk.Label(frame3, text="Precision for excited states in %", bg=coloruw,
    fg='white', font=("Helvetica", 10)).place(x=20, y=170, width=270, 
    height=20)
tk.Entry(frame3, bg=colorbg, fg='black', textvariable=precision,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=190, width=270, height=25)
# tooltip
q9 = tk.Label(frame3, image=question, bg=coloruw)
q9.place(x=403, y=175, width=10, height=10)
q9.bind("<Enter>", lambda event, arg=gI.q9: hover3(event, arg))
q9.bind("<Leave>", lambda event, arg=""   : hover3(event, arg))

# max microstates
max_micro = tk.IntVar(value=50)
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=310, y=170)
tk.Label(frame3, text="Maximum number of microstates", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=170, width=270, height=20)
tk.Entry(frame3, bg=colorbg, fg='black', textvariable=max_micro,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=310, y=190, width=270, height=25)
# tooltip
q10 = tk.Label(frame3, image=question, bg=coloruw)
q10.place(x=550, y=175, width=10, height=10)
q10.bind("<Enter>", lambda event, arg=gI.q10: hover3(event, arg))
q10.bind("<Leave>", lambda event, arg=""    : hover3(event, arg))

# Whether to include adiabatic decay in case of negative VDE
adiaDecay = tk.BooleanVar(value=False)
tk.Frame(frame3, bg=colorbg, width=270, height=45).place(x=165, y=225)
tk.Label(frame3, text="Include adiabatic decay for neg. VDE", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=165, y=225, width=270, height=20)
rbAdia1 = tk.Radiobutton(frame3, text='True', value=True,
    variable=adiaDecay, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9))
rbAdia1.place(x=165, y=245, width=135, height=25)
rbAdia2 = tk.Radiobutton(frame3, text='False', value=False,
    variable=adiaDecay, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9))
rbAdia2.place(x=300, y=245, width=135, height=25)

# tooltips frame
tk.Frame(frame3, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip3 = tk.Label(frame3, font=("Helvetica", 9),
    text='', bg=colorbg, fg=coloruw, anchor='nw', justify='left')
tooltip3.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton3 = tk.Button(frame3, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton3['image'] = quitButton
quitbutton3.place(x=20, y=340, width=70, height=70)


################################################################################
#                                  Frame 4                                     #
################################################################################
tk.Label(frame4, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover4(event, text):
    tooltip4["text"] = text

def disableExact():
    if orthotype.get() == "mo":
        rbExact1['state'] = 'normal'
        rbExact2['state'] = 'normal'
        disableInt()
    else:
        rbExact1['state'] = 'disabled'
        rbExact2['state'] = 'disabled'
        intNkOpt['state'] = 'disabled'
        intkSkipEntry['state'] = 'disabled'

def disableInt():
    if exact2elInt.get() == 'True':
        intNkOpt['state']      = 'normal'
        intkSkipEntry['state'] = 'normal'
    else:
        intNkOpt['state']      = 'disabled'
        intkSkipEntry['state'] = 'disabled'

# headline
tk.Label(frame4,
    text="Options for two-electron integrals", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame4, bg=coloruw, width=270, height=2).place(x=165, y=43)

# type of orthogonalization
orthotype = tk.StringVar(value="mo")
tk.Frame(frame4, bg=colorbg, width=270, height=45).place(x=165, y=60)
tk.Label(frame4, text="Type of orthogonalization", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=165, y=60, width=270, height=20)
rbOrthoM = tk.Radiobutton(frame4, text='MO', value='mo', variable=orthotype,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableExact())
rbOrthoM.place(x=165, y=80, width=67, height=25)
rbOrthoS = tk.Radiobutton(frame4, text='State', value='state',
    variable=orthotype, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableExact())
rbOrthoS.place(x=232, y=80, width=68, height=25)
rbOrthoD = tk.Radiobutton(frame4, text='Dyson', value='dyson',
    variable=orthotype, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableExact())
rbOrthoD.place(x=300, y=80, width=68, height=25)
rbOrthoN = tk.Radiobutton(frame4, text='None', value='none', variable=orthotype,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableExact())
rbOrthoN.place(x=368, y=80, width=67, height=25)
# tooltip
q24 = tk.Label(frame4, image=question, bg=coloruw)
q24.place(x=384, y=65, width=10, height=10)
q24.bind("<Enter>", lambda event, arg=gI.q24: hover4(event, arg))
q24.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))
rbOrthoM.bind("<Enter>", lambda event, arg=gI.q25: hover4(event, arg))
rbOrthoM.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))
rbOrthoS.bind("<Enter>", lambda event, arg=gI.q26: hover4(event, arg))
rbOrthoS.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))
rbOrthoD.bind("<Enter>", lambda event, arg=gI.q27: hover4(event, arg))
rbOrthoD.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))
rbOrthoN.bind("<Enter>", lambda event, arg=gI.q28: hover4(event, arg))
rbOrthoN.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))

# exact integrals
exact2elInt = tk.StringVar(value='False')
tk.Frame(frame4, bg=colorbg, width=270, height=45).place(x=165, y=115)
tk.Label(frame4, text="Calculation for 4-center integrals", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=165, y=115, width=270, height=20)
rbExact1 = tk.Radiobutton(frame4, text='Exact', value='True',
    variable=exact2elInt, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableInt())
rbExact1.place(x=165, y=135, width=135, height=25)
rbExact2 = tk.Radiobutton(frame4, text='Approximate', value='False',
    variable=exact2elInt, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9),
    command=lambda: disableInt())
rbExact2.place(x=300, y=135, width=135, height=25)
# tooltip
q13 = tk.Label(frame4, image=question, bg=coloruw)
q13.place(x=407, y=120, width=10, height=10)
q13.bind("<Enter>", lambda event, arg=gI.q13: hover4(event, arg))
q13.bind("<Leave>", lambda event, arg=""    : hover4(event, arg))

# number of k vectors
intNk = tk.IntVar(value=12)
tk.Frame(frame4, bg=colorbg, width=270, height=45).place(x=165, y=170)
tk.Label(frame4, text="Number of k vectors per energy", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=165, y=170, width=270, height=20)
intNkOpt = tk.OptionMenu(frame4, intNk, *[2, 6, 8, 12])
intNkOpt.configure(state='disabled', bg=colorbg, fg='black',
    activebackground=colorab, activeforeground='black', borderwidth=0,
    font=("Helvetica", 9), highlightthickness=0)
intNkOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
intNkOpt.place(x=165, y=190, width=270, height=25)

# number of k energies
intkSkip = tk.IntVar(value=20)
tk.Frame(frame4, bg=colorbg, width=270, height=45).place(x=165, y=225)
tk.Label(frame4,
    text="Number of k energies", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=165, y=225, width=270, height=20)
intkSkipEntry = tk.Entry(frame4, bg=colorbg, fg='black', textvariable=intkSkip,
    font=("Helvetica", 9), justify='center', borderwidth=0,
    highlightthickness=0, state='disabled')
intkSkipEntry.place(x=165, y=245, width=270, height=25)

# tooltips frame
tk.Frame(frame4, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip4 = tk.Label(frame4, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip4.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton4 = tk.Button(frame4, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton4['image'] = quitButton
quitbutton4.place(x=20, y=340, width=70, height=70)


################################################################################
#                                  Frame 5                                     #
################################################################################
tk.Label(frame5, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover5(event, text):
    tooltip5["text"] = text

def refreshSelection(nGrids):
    selected.set(1)
    selectOpt['menu'].delete(0,'end')
    try:
        nGrids = int(nGrids.get())
    except:
        return
    for i in range(1,nGrids+1):
        selectOpt['menu'].add_command(label=i,
                                      command=lambda name=i: selection(name))
    if len(gridConf) < nGrids:
        for _ in range(nGrids-len(gridConf)):
            gridConf.append(['fib','1.5','1000','96', False, 0.0, 0.0])

def selection(name):
    selected.set(name)
    PWtype.set(gridConf[int(name)-1][0])
    maxEn.set(gridConf[int(name)-1][1])
    nEnergy.set(gridConf[int(name)-1][2])
    nkPerEn.set(gridConf[int(name)-1][3])
    zpeInclude.set(gridConf[int(name)-1][4])
    zpeDiff.set(gridConf[int(name)-1][5])
    vibExEn.set(gridConf[int(name)-1][6])
    if PWtype.get() != 'fib':
        fibEntry['state'] = 'disabled'
    else:
        fibEntry['state'] = 'normal'

def changeType(pwtype, selected):
    pwtype   = PWtype.get()
    selected = selected.get()
    gridConf[selected-1][0] = pwtype
    if pwtype == 'cubic':
        gridConf[selected-1][1] = '0.1,0.1,0.1'
        gridConf[selected-1][2] = '200,200,200'
        gridConf[selected-1][3] = ''
        maxEn.set('0.1,0.1,0.1')
        nEnergy.set('200,200,200')
        nkPerEn.set('')
        klabel['text']    = 'Maximum values of kx,ky,kz (in a.u.)'
        nElabel['text']   = 'Number of values for kx,ky,kz'
        fibEntry['state'] = 'disabled'
        zpeI1['state']    = 'disabled'
        zpeI2['state']    = 'disabled'
        zpeInclude.set(False)
    else:
        gridConf[selected-1][1] = '1.5'
        gridConf[selected-1][2] = '200,200,200'
        maxEn.set('1.5')
        nEnergy.set('1000')
        klabel['text']  = 'Maximum plane wave energy (in eV)'
        nElabel['text'] = 'Number of plane wave energies'
        zpeI1['state']  = 'normal'
        zpeI2['state']  = 'normal'
        if pwtype == 'snub':
            gridConf[selected-1][3] = '24'
            nkPerEn.set('24')
            fibEntry['state'] = 'disabled'
        else:
            gridConf[selected-1][3] = '96'
            nkPerEn.set('96')
            fibEntry['state'] = 'normal'

def changeMaxEn(maxEn, pwtype, selected):
    maxEn   = maxEn.get()
    pwtype   = pwtype.get()
    selected = selected.get()
    if pwtype == 'fib' or pwtype == 'snub':
        try:
            float(maxEn)
        except:
            return
    elif pwtype == 'cubic':
        try:
            [float(i) for i in maxEn.split(',')]
            if len(maxEn.split(',')) != 3:
                return
        except:
            return
    gridConf[selected-1][1] = maxEn

def changeNEnergy(nEnergy, pwtype, selected):
    nEnergy  = nEnergy.get()
    pwtype   = pwtype.get()
    selected = selected.get()
    if pwtype == 'fib' or pwtype == 'snub':
        try:
            int(nEnergy)
        except:
            return
    elif pwtype == 'cubic':
        try:
            [float(i) for i in nEnergy.split(',')]
            if len(nEnergy.split(',')) != 3:
                return
        except:
            return
    gridConf[selected-1][2] = nEnergy

def changeNkPerEn(nkPerEn, pwtype, selected):
    nkPerEn  = nkPerEn.get()
    pwtype   = pwtype.get()
    selected = selected.get()
    try:
        int(nkPerEn)
    except:
        return
    gridConf[selected-1][3] = nkPerEn

def changeZpeInclude(zpeInclude, selected):
    zpeInclude = zpeInclude.get()
    selected   = selected.get()
    gridConf[selected-1][4] = zpeInclude
    if zpeInclude:
        maxEntry['state'] = 'disabled'
        zpeEntry['state'] = 'normal'
        vibEntry['state'] = 'normal'
        maxEn.set(zpeDiff.get() + vibExEn.get())
    else:
        maxEntry['state'] = 'normal'
        zpeEntry['state'] = 'disabled'
        vibEntry['state'] = 'disabled'

def changeZpeDiff(zpeDiff, vibExEn, selected):
    zpeDiff  = zpeDiff.get()
    vibExEn  = vibExEn.get()
    selected = selected.get()
    gridConf[selected-1][5] = zpeDiff
    maxEn.set(zpeDiff + vibExEn)


def changeVibExEn(zpeDiff, vibExEn, selected):
    zpeDiff  = zpeDiff.get()
    vibExEn  = vibExEn.get()
    selected = selected.get()
    gridConf[selected-1][6] = vibExEn
    maxEn.set(zpeDiff + vibExEn)

gridConf = [['fib', '1.5', '1000', '96', False, 0.0, 0.0]]

# headline
tk.Label(frame5,
    text="Options for free electron grid", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame5, bg=coloruw, width=270, height=2).place(x=165, y=43)

# number of different settings
nGrids = tk.IntVar(value='1')
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=20, y=60)
tk.Label(frame5, text="Number of different grids", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=60, width=270, height=20)
tk.Entry(frame5, bg=colorbg, fg='black', justify='center', textvariable=nGrids,
    font=("Helvetica", 9), borderwidth=0, highlightthickness=0).place(
    x=20, y=80, width=270, height=25)
nGrids.trace_add('write', lambda a,b,c,d=nGrids: refreshSelection(d))
# tooltip
q14 = tk.Label(frame5, image=question, bg=coloruw)
q14.place(x=237, y=65, width=10, height=10)
q14.bind("<Enter>", lambda event, arg=gI.q14: hover5(event, arg))
q14.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# selected grid
selected = tk.IntVar(value=1)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=310, y=60)
tk.Label(frame5, text="Selected grid", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=60, width=270, height=20)
selectOpt = tk.OptionMenu(frame5, selected, *[1])
selectOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
selectOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
selectOpt.place(x=310, y=80, width=270, height=25)

# type of discretization
PWtype = tk.StringVar(value='fib')
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=165, y=115)
tk.Label(frame5, text="Type of discretization", bg=coloruw, fg='white',
    font=("Helvetica",10)).place(x=165, y=115, width=270, height=20)
tk.Radiobutton(frame5, text="Fibonacci", value="fib", font=("Helvetica", 9),
    variable=PWtype, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black',
    command=lambda b=PWtype, c=selected: changeType(PWtype, selected)
    ).place(x=165, y=135, width=90, height=25)
tk.Radiobutton(frame5, text="Snub", value="snub", variable=PWtype, bg=colorbg,
    fg='black', activebackground=colorab, activeforeground='black',
    command=lambda b=PWtype, c=selected: changeType(PWtype, selected),
    highlightthickness=0, font=("Helvetica", 9)
    ).place(x=300, y=135, width=90, height=25, anchor='n')
tk.Radiobutton(frame5, text="Cubic", value="cubic", variable=PWtype, bg=colorbg,
    fg='black', activebackground=colorab, activeforeground='black',
    command=lambda b=PWtype, c=selected: changeType(PWtype, selected),
    highlightthickness=0, font=("Helvetica", 9)
    ).place(x=435, y=135, width=90, height=25, anchor='ne')
# tooltip
q15 = tk.Label(frame5, image=question, bg=coloruw)
q15.place(x=370, y=120, width=10, height=10)
q15.bind("<Enter>", lambda event, arg=gI.q15: hover5(event, arg))
q15.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# number of different pw energies ('fib','snub') or number of k_i
nEnergy = tk.StringVar(value=1000)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=20, y=170)
nElabel = tk.Label(frame5, text="Number of plane wave energies", bg=coloruw,
    fg='white', font=("Helvetica", 10))
nElabel.place(x=20, y=170, width=270, height=20)
tk.Entry(frame5, bg=colorbg, fg='black', textvariable=nEnergy,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=20, y=190, width=270, height=25)
nEnergy.trace_add('write',
    lambda a,b,c,d=nEnergy: changeNEnergy(d, PWtype, selected))

# number of different k per energy ('fib')
nkPerEn = tk.StringVar(value=96)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=310, y=170)
tk.Label(frame5, text="Number of plane waves per energy", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=310, y=170, width=270, height=20)
fibEntry = tk.Entry(frame5, bg=colorbg, fg='black', textvariable=nkPerEn,
    justify='center', font=("Helvetica", 9), borderwidth=0,
    highlightthickness=0)
fibEntry.place(x=310, y=190, width=270, height=25)
nkPerEn.trace_add('write',
    lambda a,b,c,d=nkPerEn: changeNkPerEn(d, PWtype, selected))

# whether to include zero-point energy
zpeInclude = tk.BooleanVar(value=False)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=20, y=225)
tk.Label(frame5, text="Include zero-point energy of states", bg=coloruw, 
    fg='white', font=("Helvetica",10)).place(x=20, y=225, width=270, height=20)
zpeI1 = tk.Radiobutton(frame5, text="Include", value=True,
    font=("Helvetica", 9), variable=zpeInclude, bg=colorbg, fg='black',
    highlightthickness=0, activebackground=colorab, activeforeground='black')
zpeI1.place(x=20, y=245, width=135, height=25)
zpeI2 = tk.Radiobutton(frame5, text="Don't Include", value=False, 
    variable=zpeInclude, bg=colorbg, fg='black', activebackground=colorab, 
    activeforeground='black', highlightthickness=0, font=("Helvetica", 9))
zpeI2.place(x=155, y=245, width=135, height=25)
zpeInclude.trace_add('write', 
    lambda a,b,c,d=zpeInclude: changeZpeInclude(d, selected))
# tooltip
q24 = tk.Label(frame5, image=question, bg=coloruw)
q24.place(x=268, y=230, width=10, height=10)
q24.bind("<Enter>", lambda event, arg=gI.q29: hover5(event, arg))
q24.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# maximum pw energy ('fib', 'snub') or maximum k_i ('cubic')
maxEn = tk.StringVar(value=1.5)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=310, y=225)
klabel = tk.Label(frame5, text="Maximum plane wave energy (in eV)", bg=coloruw,
    fg='white', font=("Helvetica", 10))
klabel.place(x=310, y=225, width=270, height=20)
maxEntry = tk.Entry(frame5, bg=colorbg, fg='black', textvariable=maxEn,
    font=("Helvetica", 9), justify='center', borderwidth=0,
    highlightthickness=0)
maxEntry.place(x=310, y=245, width=270, height=25)
maxEn.trace_add('write', lambda a,b,c,d=maxEn: changeMaxEn(d, PWtype, selected))
# tooltip
q16 = tk.Label(frame5, image=question, bg=coloruw)
q16.place(x=562, y=230, width=10, height=10)
q16.bind("<Enter>", lambda event, arg=gI.q16: hover5(event, arg))
q16.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# Difference of ZPE between initial and neutral ground state in eV
zpeDiff = tk.DoubleVar(value=0.0)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=20, y=280)
zlabel = tk.Label(frame5, text="ZPE difference (in eV)", 
    bg=coloruw, fg='white', font=("Helvetica", 10))
zlabel.place(x=20, y=280, width=270, height=20)
zpeEntry = tk.Entry(frame5, bg=colorbg, fg='black', textvariable=zpeDiff,
    font=("Helvetica", 9), justify='center', borderwidth=0, state='disabled',
    highlightthickness=0)
zpeEntry.place(x=20, y=300, width=270, height=25)
zpeDiff.trace_add('write', 
    lambda a,b,c,d=zpeDiff: changeZpeDiff(d, vibExEn, selected))
# tooltip
q25 = tk.Label(frame5, image=question, bg=coloruw)
q25.place(x=228, y=285, width=10, height=10)
q25.bind("<Enter>", lambda event, arg=gI.q30: hover5(event, arg))
q25.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# energy of excited vibrational mode in eV
vibExEn = tk.DoubleVar(value=0.0)
tk.Frame(frame5, bg=colorbg, width=270, height=45).place(x=310, y=280)
vlabel = tk.Label(frame5, text="Vibrational excitation energy (in eV)", 
    bg=coloruw, fg='white', font=("Helvetica", 10))
vlabel.place(x=310, y=280, width=270, height=20)
vibEntry = tk.Entry(frame5, bg=colorbg, fg='black', textvariable=vibExEn,
    font=("Helvetica", 9), justify='center', borderwidth=0, state='disabled',
    highlightthickness=0)
vibEntry.place(x=310, y=300, width=270, height=25)
vibExEn.trace_add('write', 
    lambda a,b,c,d=zpeDiff: changeVibExEn(d, vibExEn, selected))
# tooltip
q26 = tk.Label(frame5, image=question, bg=coloruw)
q26.place(x=560, y=285, width=10, height=10)
q26.bind("<Enter>", lambda event, arg=gI.q31: hover5(event, arg))
q26.bind("<Leave>", lambda event, arg=""    : hover5(event, arg))

# tooltips frame
tk.Frame(frame5, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip5 = tk.Label(frame5, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip5.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton5 = tk.Button(frame5, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton5['image'] = quitButton
quitbutton5.place(x=20, y=340, width=70, height=70)


################################################################################
#                                  Frame 6                                     #
################################################################################
tk.Label(frame6, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover6(event, text):
    tooltip6["text"] = text

# headline
tk.Label(frame6, text="Options for folder management", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame6, bg=coloruw, width=270, height=2).place(x=165, y=43)

# number of trajectories
nrTraj = tk.IntVar(value=100)
tk.Frame(frame6, bg=colorbg, width=270, height=45).place(x=165, y=60)
tk.Label(frame6, text="Number of trajectories", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=165, y=60, width=270, height=20)
tk.Entry(frame6, bg=colorbg, fg='black', textvariable=nrTraj,
    font=("Helvetica", 9), justify='center', borderwidth=0, highlightthickness=0
    ).place(x=165, y=80, width=270, height=25)
# tooltip
q17 = tk.Label(frame6, image=question, bg=coloruw)
q17.place(x=375, y=65, width=10, height=10)
q17.bind("<Enter>", lambda event, arg=gI.q17: hover6(event, arg))
q17.bind("<Leave>", lambda event, arg=""    : hover6(event, arg))

# create folder structure
createFolders = tk.BooleanVar(value=False)
tk.Frame(frame6, bg=colorbg, width=270, height=45).place(x=165, y=115)
tk.Label(frame6, text="Create dynamics folder structure", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=165, y=115, width=270, height=20)
tk.Radiobutton(frame6, text='Yes', value=True, variable=createFolders,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=165, y=135, width=135, height=25)
tk.Radiobutton(frame6, text='No', value=False, variable=createFolders,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=300, y=135, width=135, height=25)
# tooltip
q18 = tk.Label(frame6, image=question, bg=coloruw)
q18.place(x=407, y=120, width=10, height=10)
q18.bind("<Enter>", lambda event, arg=gI.q18: hover6(event, arg))
q18.bind("<Leave>", lambda event, arg=""    : hover6(event, arg))

# copy everything to folders
cpFiles = tk.BooleanVar(value=False)
tk.Frame(frame6, bg=colorbg, width=270, height=45).place(x=165, y=170)
tk.Label(frame6, text="Copy input to folders", bg=coloruw, fg='white',
    font=("Helvetica",10)).place(x=165, y=170, width=270, height=20)
tk.Radiobutton(frame6, text='Yes', value=True, variable=cpFiles, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=165, y=190, width=135, height=25)
tk.Radiobutton(frame6, text='No', value=False, variable=cpFiles, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=300, y=190, width=135, height=25)
# tooltip
q19 = tk.Label(frame6, image=question, bg=coloruw)
q19.place(x=371, y=175, width=10, height=10)
q19.bind("<Enter>", lambda event, arg=gI.q19: hover6(event, arg))
q19.bind("<Leave>", lambda event, arg=""    : hover6(event, arg))

# write submit script
subScript = tk.BooleanVar(value=False)
tk.Frame(frame6, bg=colorbg, width=270, height=45).place(x=165, y=225)
tk.Label(frame6, text="Write submit script for dynamics", bg=coloruw,
    fg='white', font=("Helvetica", 10)
    ).place(x=165, y=225, width=270, height=20)
tk.Radiobutton(frame6, text='Yes', value=True, variable=subScript, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=165, y=245, width=135, height=25)
tk.Radiobutton(frame6, text='No', value=False, variable=subScript, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9)
    ).place(x=300, y=245, width=135, height=25)
# tooltip
q20 = tk.Label(frame6, image=question, bg=coloruw)
q20.place(x=404,y=230,width=10,height=10)
q20.bind("<Enter>", lambda event, arg=gI.q20: hover6(event, arg))
q20.bind("<Leave>", lambda event, arg=""    : hover6(event, arg))

# tooltips frame
tk.Frame(frame6, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip6 = tk.Label(frame6, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip6.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton6 = tk.Button(frame6, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton6['image'] = quitButton
quitbutton6.place(x=20, y=340, width=70, height=70)


################################################################################
#                                  Frame 7                                     #
################################################################################
tk.Label(frame7, image=logo, bg='white').place(x=0, y=0, width=600, height=450)

def hover7(event, text):
    tooltip7["text"] = text

def activateDist(doWigner):
    doWigner = doWigner.get()
    if doWigner:
        linWignerRB1["state"]  = "normal"
        linWignerRB2["state"]  = "normal"
        rbDist1["state"]       = "normal"
        rbDist3["state"]       = "normal"
        inConEntry["state"]    = "normal"
        modeTempEntry["state"] = "normal"
        WigQCOpt["state"]      = "normal"
        inpEntry1["state"]     = "normal"
        inpEntry2["state"]     = "normal"
    else:
        linWignerRB1["state"]  = "disabled"
        linWignerRB2["state"]  = "disabled"
        rbDist1["state"]       = "disabled"
        rbDist3["state"]       = "disabled"
        inConEntry["state"]    = "disabled"
        modeTempEntry["state"] = "disabled"
        WigQCOpt["state"]      = "disabled"
        inpEntry1["state"]     = "disabled"
        inpEntry2["state"]     = "disabled"

def changeDist(distType):
    distType = distType.get()
    if distType == "v1":
        modeTemp["text"] = "Excited Modes (comma-separated)"
        modes.set("")
    elif distType == "w":
        modeTemp["text"] = "Temperature (in K)"
        modes.set("100")

def changeQC(event, qc):
    qc = qc.get()
    if qc == 'g09':
        input1.set("freq.fchk")
        input2.set("freq.log")
        inpLabel1["text"] = "Formatted g09 checkpoint file"
        inpLabel2["text"] = "g09 standard output file"
        inpEntry2["state"] = "normal"
    elif qc == "QChem":
        input1.set("freq.log")
        input2.set("HESS")
        inpLabel1["text"] = "QChem standard ouput file"
        inpLabel2["text"] = "QChem HESS file (generated with -save)"
        inpEntry2["state"] = "normal"

# headline
tk.Label(frame7,
    text="Options for Wigner distribution", fg=coloruw, bg=colorbg,
    font="Helvetica 12").place(x=165, y=10, width=270, height=35)
tk.Frame(frame7, bg=coloruw, width=270, height=2).place(x=165, y=43)

# Wigner
doWigner = tk.BooleanVar(value=False)
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=20, y=60)
tk.Label(frame7, text="Generate initial conditions", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=60, width=270, height=20)
tk.Radiobutton(frame7, text='Yes', value=True, variable=doWigner, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9),
    command=lambda c=doWigner: activateDist(c)).place(
    x=20, y=80, width=135, height=25)
tk.Radiobutton(frame7, text='No', value=False, variable=doWigner, bg=colorbg,
    fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 9),
    command=lambda c=doWigner: activateDist(c)).place(
    x=155, y=80, width=135, height=25)
# tooltip
q21 = tk.Label(frame7, image=question, bg=coloruw)
q21.place(x=240, y=65, width=10, height=10)
q21.bind("<Enter>", lambda event, arg=gI.q21: hover7(event, arg))
q21.bind("<Leave>", lambda event, arg=""    : hover7(event, arg))

# linear molecule
linWigner = tk.BooleanVar(value=False)
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=310, y=60)
tk.Label(frame7, text="Molecule is linear", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=60, width=270, height=20)
linWignerRB1 = tk.Radiobutton(frame7, state="disabled", text='Yes', value=True,
    variable=linWigner, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9))
linWignerRB1.place(x=310, y=80, width=135, height=25)
linWignerRB2 = tk.Radiobutton(frame7, state="disabled", text='No', value=False,
    variable=linWigner, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica",9))
linWignerRB2.place(x=445, y=80, width=135, height=25)

# type of distribution
distType = tk.StringVar(value='v1')
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=20, y=115)
tk.Label(frame7,
    text="Type of distribution", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=20, y=115, width=270, height=20)
rbDist1 = tk.Radiobutton(frame7, state="disabled",
    text='|\u03a8(x)|\u00B2|\u03a8(q)|\u00B2', value='v1', variable=distType,
    bg=colorbg, fg='black', highlightthickness=0, activebackground=colorab,
    activeforeground='black', font=("Helvetica", 10),
    command=lambda c=distType: changeDist(c))
rbDist1.place(x=20, y=135, width=135, height=25)
rbDist3 = tk.Radiobutton(frame7, state="disabled", text='Wigner', value='w',
    variable=distType, bg=colorbg, fg='black', highlightthickness=0,
    activebackground=colorab, activeforeground='black', font=("Helvetica", 9),
    command=lambda c=distType: changeDist(c))
rbDist3.place(x=155, y=135, width=135, height=25)
# tooltip
q22 = tk.Label(frame7, image=question, bg=coloruw)
q22.place(x=218, y=120, width=10, height=10)
q22.bind("<Enter>", lambda event, arg=gI.q22: hover7(event, arg))
q22.bind("<Leave>", lambda event, arg=""    : hover7(event, arg))

# number of initial conditions
nrInCon = tk.IntVar(value=100)
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=310, y=115)
tk.Label(frame7, text="Number of initial conditions", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=310, y=115, width=270, height=20)
inConEntry = tk.Entry(frame7, bg=colorbg, fg='black', textvariable=nrInCon,
    borderwidth=0, font=("Helvetica", 9), justify='center',
    highlightthickness=0, state="disabled")
inConEntry.place(x=310, y=135, width=270, height=25)

# excited modes (h or v1) or temperature (w)
modes = tk.StringVar(value="")
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=165, y=170)
modeTemp = tk.Label(frame7, text="Excited Modes (comma-separated)", bg=coloruw,
    fg='white', font=("Helvetica", 10))
modeTemp.place(x=165, y=170, width=270, height=20)
modeTempEntry = tk.Entry(frame7, bg=colorbg, fg='black', textvariable=modes,
    borderwidth=0, font=("Helvetica", 9), justify='center',
    highlightthickness=0, state="disabled")
modeTempEntry.place(x=165, y=190, width=270, height=25)
# tooltip
q23 = tk.Label(frame7, image=question, bg=coloruw)
q23.place(x=414, y=175, width=10, height=10)
q23.bind("<Enter>", lambda event, arg=gI.q23: hover7(event, arg))
q23.bind("<Leave>", lambda event, arg=""    : hover7(event, arg))

# qc program for Wigner
qcWigner = tk.StringVar(value='g09')
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=165, y=225)
tk.Label(frame7, text="QC program for distribution", bg=coloruw, fg='white',
    font=("Helvetica", 10)).place(x=165, y=225, width=270, height=20)
WigQCOpt = tk.OptionMenu(frame7, qcWigner, *["g09", "QChem"], 
    command=lambda event, arg=qcWigner: changeQC(event, arg))
WigQCOpt.configure(bg=colorbg, fg='black', activebackground=colorab,
    activeforeground='black', borderwidth=0, highlightthickness=0,
    font=("Helvetica", 9))
WigQCOpt['menu'].configure(bg='white', fg='black', activebackground=colorab,
    activeforeground='black')
WigQCOpt["state"] = "disabled"
WigQCOpt.place(x=165, y=245, width=270, height=25)

# input1
input1 = tk.StringVar(value="freq.fchk")
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=20, y=280)
inpLabel1 = tk.Label(frame7, text="Formatted g09 checkpoint file", bg=coloruw,
    fg='white', font=("Helvetica", 10))
inpLabel1.place(x=20, y=280, width=270, height=20)
inpEntry1 = tk.Entry(frame7, bg=colorbg, fg='black', textvariable=input1,
    borderwidth=0, font=("Helvetica", 9), justify='center',
    highlightthickness=0, state="disabled")
inpEntry1.place(x=20, y=300, width=270, height=25)

# input2
input2 = tk.StringVar(value="freq.log")
tk.Frame(frame7, bg=colorbg, width=270, height=45).place(x=310, y=280)
inpLabel2 = tk.Label(frame7, text="g09 standard output file", bg=coloruw,
    fg='white', font=("Helvetica", 10))
inpLabel2.place(x=310, y=280, width=270, height=20)
inpEntry2 = tk.Entry(frame7, bg=colorbg, fg='black', textvariable=input2,
    borderwidth=0, font=("Helvetica", 9), justify='center',
    highlightthickness=0, state="disabled")
inpEntry2.place(x=310, y=300, width=270, height=25)

# tooltips frame
tk.Frame(frame7, bg=colorbg).place(x=110, y=340, width=470, height=70)
tooltip7 = tk.Label(frame7, font=("Helvetica", 9), text='', bg=colorbg,
    fg=coloruw, anchor='nw', justify='left')
tooltip7.place(x=115, y=345, width=460, height=60)

# quit button
quitbutton7 = tk.Button(frame7, borderwidth=0, highlightthickness=0,
    background='white', command=confirm)
quitbutton7['image'] = quitButton
quitbutton7.place(x=20, y=340, width=70, height=70)

root.mainloop()
