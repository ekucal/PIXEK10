from tkinter import Toplevel
from tkinter import messagebox
import tkinter as tk
from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler
from modules.data import Spectra, list_of_spectra, list_of_spectra_c 

def add_to_list(self, name_):
    list_of_spectra.append(Spectra(name_, self.spectrum.counts, self.spectrum.channel, self.spectrum.energy, self.spectrum.material, self.spectrum.crystal_orientation, self.spectrum.axial_channel, self.spectrum.tilt_angle, self.spectrum.filter, self.spectrum.background, self.spectrum.final_yield))
    messagebox.showinfo('Info', 'Done!')

def add_to_list_c(self):
    list_of_spectra_c.append(Spectra(self.spectrum.tilt_angle, self.spectrum.counts, self.spectrum.channel, self.spectrum.energy, self.spectrum.material, self.spectrum.crystal_orientation, self.spectrum.axial_channel, self.spectrum.tilt_angle, self.spectrum.filter, self.spectrum.background, self.spectrum.final_yield))
    messagebox.showinfo('Info', 'Done!')

def show_list(self):
    self.list_for_plot.delete(0, tk.END)
    for i in list_of_spectra:
        self.list_for_plot.insert(tk.END, i.name)

def show_plot_list(self):
    plt.clf()
    selected_items = self.list_for_plot.curselection()
    data = [self.list_for_plot.get(i) for i in selected_items]
    for i in list_of_spectra:
        ev_data = str(data)
        ev_a = str(i.name)
        if ev_a in ev_data:
            plt.plot(i.energy, i.final_yield, '-', label=i.name, markersize=0.5) 
    plt.yscale('log')
    plt.xlabel('Energy')
    plt.ylabel('Yield')
    plt.legend()
    self.canvas_1.draw()

def show_plot_list_c(self):
    plt.clf()
    for i in list_of_spectra_c:
        plt.plot(i.energy, i.counts, '.', label=i.tilt_angle, markersize=0.8) 
    plt.yscale('log')
    plt.xlabel('Energy')
    plt.ylabel('Yield')
    plt.legend()
    self.canvas_1.draw()

def show_plot_list_c_b(self):
    plt.clf()
    for i in list_of_spectra_c:
        plt.plot(i.energy, i.final_yield, '-', label=i.tilt_angle, markersize=0.5) 
    plt.yscale('log')
    plt.xlabel('Energy')
    plt.ylabel('Yield')
    plt.legend()
    self.canvas_1.draw()