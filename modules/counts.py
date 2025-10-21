from tkinter import Toplevel
from tkinter import messagebox
import tkinter as tk
from tkinter import *

import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler
from math import sin, pi, exp, sqrt
import numpy as np
import csv
from modules.data import Spectra, list_of_spectra, Counts, list_of_counts, list_of_spectra_c

def show_list_for_counts(self):
    self.list_for_counts.delete(0, tk.END)
    for i in list_of_spectra:
        self.list_for_counts.insert(tk.END, i.name)

def show_plot_counts(self):
    plt.clf()
    list_of_counts.clear()
    selected_items = self.list_for_counts.curselection()
    data = [self.list_for_counts.get(i) for i in selected_items]
    a_min = float(self.emin_as.get())
    b_max = float(self.emax_as.get())
    i1 = 1
    i1_error = 0
    N = []
    C = []
    Error = []
    if int(self.SC.get()) == 0:
        for i in list_of_spectra:
            ev_data = str(data)
            ev_a = str(i.name)
            if ev_a in ev_data:
                as_counts = 0
                err = 0
                for o in range(0, len(i.energy[:])):
                    if i.energy[o] > a_min and i.energy[o] < b_max:
                        as_counts = as_counts+float(i.final_yield[o])
                        err = sqrt(float(i.counts[o]) + float(i.background[o]))
                N.append(i.name)
                C.append(as_counts)
                Error.append(err)
                list_of_counts.append(Counts(i.name, as_counts, err))
    else:
        for i in list_of_spectra:
            ev_data = str(data)
            ev_a = str(i.name)
            if ev_a in ev_data:
                as_counts = 0
                as_y = 0
                as_b = 0
                err = 0
                err_2 = 0
                for o in range(0, len(i.energy[:])):
                    if i.energy[o] > a_min and i.energy[o] < b_max:
                        as_counts = as_counts + float(i.final_yield[o])
                        as_y = as_y + float(i.counts[o])
                        as_b = as_b + float(i.background[o])
                err = sqrt(as_y + as_b)
                if i1 == 1:
                    i1 = as_counts
                    i1_error = err
                err_2 = (as_counts/i1) * ((err/as_counts) + (i1_error/i1))
                N.append(i.name)
                C.append(as_counts/i1)
                Error.append(err_2)
                list_of_counts.append(Counts(i.name, as_counts/i1, err_2))
    if int(self.EB.get()) == 1:
        plt.bar(N, C, color='black')
        plt.errorbar(N, C, yerr=Error, fmt="o", color="r")
    else:
        plt.bar(N, C, color='orange')
    plt.ylabel('Counts')
    self.canvas_1.draw()

def show_plot_counts_2(self):
    plt.clf()
    list_of_counts.clear()
    a_min = float(self.emin_as.get())
    b_max = float(self.emax_as.get())
    i1 = 1
    i1_error = 0
    N = []
    C = []
    Error = []
    if int(self.SC.get()) == 0:
        for i in list_of_spectra_c:
            as_counts = 0
            err = 0
            for o in range(0, len(i.energy[:])):
                if i.energy[o] > a_min and i.energy[o] < b_max:
                    as_counts = as_counts+float(i.final_yield[o])
                    err = sqrt(float(i.counts[o]) + float(i.background[o]))
            N.append(i.name)
            C.append(as_counts)
            Error.append(err)
            list_of_counts.append(Counts(i.name, as_counts, err))
    else:
        for i in list_of_spectra_c:
            as_counts = 0
            as_y = 0
            as_b = 0
            err = 0
            err_2 = 0
            for o in range(0, len(i.energy[:])):
                if i.energy[o] > a_min and i.energy[o] < b_max:
                    as_counts = as_counts + float(i.final_yield[o])
                    as_y = as_y + float(i.counts[o])
                    as_b = as_b + float(i.background[o])
            err = sqrt(as_y + as_b)
            if i1 == 1:
                i1 = as_counts
                i1_error = err
            err_2 = (as_counts/i1) * ((err/as_counts) + (i1_error/i1))
            N.append(i.name)
            C.append(as_counts/i1)
            Error.append(err_2)
            list_of_counts.append(Counts(i.name, as_counts/i1, err_2))
    if int(self.EB.get()) == 1:
        plt.bar(N, C, color='black')
        plt.errorbar(N, C, yerr=Error, fmt="o", color="r")
    else:
        plt.bar(N, C, color='orange')
    plt.ylabel('Counts')
    self.canvas_1.draw()

def save_counts(self, file_name):
    with open(file_name, "w") as file:
        file.write("Name Counts Error \n")
        for i in list_of_counts:
            file.write(str(i.name) + " ")
            file.write(str(i.counts) + " ")
            file.write(str(i.error) + "\n")
    messagebox.showinfo('Info', 'Done!')

