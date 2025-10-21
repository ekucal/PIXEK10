
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler
import tkinter as tk
from tkinter import *
from tkinter import Toplevel
from tkinter import messagebox
import threading
from modules.data import Spectra, list_of_spectra, list_of_spectra_c
import numpy as np
from math import ceil

def plot_spectra(self):
    plt.clf()
    self.new_window_plot = 0
    plt.plot(self.spectrum.channel, self.spectrum.counts, '.', markersize=2) 
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    plt.minorticks_on()
    plt.yscale('log')
    plt.xlabel('Channels')
    plt.ylabel('Yield')
    self.canvas_1.draw()


def plot_background(self, rb):
    self.new_window_plot = 1
    plt.plot(self.spectrum.energy, self.spectrum.background, '--', label=rb)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    plt.minorticks_on()
    plt.yscale('log')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Yield')
    plt.legend()
    self.canvas_1.draw() 

def plot_background_list(self, root):
    new_window_thread = threading.Thread(target=plot_background_list2(self, root))
    new_window_thread.start() 

def plot_background_list2(self, root):
    chart_window = Toplevel(root)
    chart_window.title("Plot")
    chart_window.geometry("1600x1000")
    dl = len(list_of_spectra_c)
    ncols_dl = int(ceil(dl/4))
    # sprawdzic pierwiastek - do ciecia do okna
    nrows=4  
    #ncols=int(np.ceil(sqrt(dl)))
    #nrows = int(np.floor(dl/ncols))

    fig2 = plt.figure(figsize=(16, 10))
    axX = 1
    axY = 1
    axN = 1
    for sp in list_of_spectra_c:
        plt.subplot(ncols_dl, nrows, axN) 
        plt.plot(sp.energy, sp.counts, '.', markersize=1, color='black')
        plt.plot(sp.energy, sp.background, '--', color='red')
        plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
        plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
        plt.minorticks_on()
        plt.yscale('log')
        plt.xlim(2.5, 8)
        plt.xlabel('Energy')
        if axN == 1 or axN == 5 or axN == 9 or axN == 13:
            plt.ylabel('Yield')
        else: 
            #plt.yticks([])
            plt.gca().set_yticklabels([])
            plt.grid(True, axis='y')
        plt.title(sp.tilt_angle, loc='center')
        axN = axN + 1
    plt.tight_layout()
    canvas_2 = tk.Canvas()
    canvas_2 = FigureCanvasTkAgg(fig2, master=chart_window) 
    canvas_2.get_tk_widget().place(x=10, y=10)  
    canvas_2.draw()
    toolbar_2 = NavigationToolbar2Tk(canvas_2, chart_window)
    toolbar_2.update()
    toolbar_2.place(x=0, y=0)

def new_window2(self, root):
    chart_window = Toplevel(root)
    chart_window.title("Plot")
    chart_window.geometry("525x475")
    fig_2 = plt.figure() 
    plt.clf()
    if self.new_window_plot == 0:
        plt.plot(self.spectrum.channel, self.spectrum.counts, '.', markersize=2)
    if self.new_window_plot == 1:
        plt.plot(self.spectrum.energy, self.spectrum.counts, '.', markersize=2)
        plt.plot(self.spectrum.energy, self.spectrum.background, '--')
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    plt.minorticks_on()
    plt.yscale('log')
    plt.xlabel('Energy')
    plt.ylabel('Yield')
    canvas_2 = tk.Canvas()
    canvas_2 = FigureCanvasTkAgg(fig_2, master=chart_window) 
    canvas_2.get_tk_widget().place(x=10, y=10)  
    canvas_2.draw()
    toolbar_2 = NavigationToolbar2Tk(canvas_2, chart_window)
    toolbar_2.update()
    toolbar_2.place(x=10, y=430)

def plot_line(self, xx):
    Y_u = []
    X_u = []
    Y_u.append(0)
    Y_u.append(1000000)
    X_u.append(xx)
    X_u.append(xx)
    plt.plot(X_u, Y_u, linestyle='-', color='red')
    self.canvas_1.draw()

def new_window(self, root):
    new_window_thread = threading.Thread(target=new_window2(self, root))
    new_window_thread.start()


    