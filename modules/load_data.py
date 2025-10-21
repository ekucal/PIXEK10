import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler
import os
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfile
from tkinter import messagebox
import csv
import json
from modules.data import Spectra, list_of_spectra, list_of_spectra_c
from modules.config import InitialParameters
from modules.plot_list import show_plot_list_c, show_plot_list_c_b

def open_folder(self, root):
    filetypes_ = (('pixie files','*.spe'),('data files', '*.dat'),('All files', '*.*'))
    file = inputFileName = askopenfile(parent=root, title='Open a file', initialdir=self.open_spectra_file_loc, filetypes=filetypes_)
    self.open_spectra_file.set(file.name) 
    self.open_spectra_file_loc = self.open_spectra_file

def open_file(self):
    data_storage = []
    data_storage.append(self.spectra_name.get())
    data_storage.append(self.angle.get())
    data_storage.append(self.dose.get())
    A = self.angle.get()
    N = self.spectra_name.get()
    D = []
    K = []
    B = []
    with open(self.open_spectra_file.get()) as f:
        lines = f.readlines()[8:828]
        e_1 = 0
        for line in lines:
            e = line.strip().split()
            e_2 = 0
            for element in e:
                D.append(float(e[e_2])*10/float(self.dose.get()))
                K.append(e_1)
                B.append(0)
                e_1 = e_1 + 1
                e_2 = e_2 + 1
    f.close()
    self.spectrum = Spectra(N, D, K, B, ' ', ' ', ' ', A, B, B, D)


def save_to_folder(self, file_name):
    with open(file_name, "w") as file:
        #file.write('# Config parameters:')
        #file.write(self.open_spectra_file.get(), self.name.get(), self.angle.get(), self.dose.get(), self.cal_a.get(), self.cal_b.get(), self.ul.get(), self.um.get(), self.SGf.get(), self.sg_min.get(), self.sg_max.get(), self.rb_entry.get(), self.wavelet_family_names_entry.get(), self.wavelet_names.get(), self.level_number.get(), self.maxiter.get(), self.oM.get(), self.a0.get(), self.a1.get(), self.a2.get(), self.a3.get(), self.a4.get(),  self.x_1.get(),  self.y_1.get(),  self.x_2.get(),  self.y_2.get())
        #file.write('# Spectrum:')
        matrix = [self.spectrum.channel, self.spectrum.energy,  self.spectrum.counts, self.spectrum.background, self.spectrum.final_yield]
        dane = [list(row) for row in zip(*matrix)]
        writer = csv.writer(file)
        writer.writerows(dane)
        messagebox.showinfo('Info', 'Done!')

def spectra_to_dict(spectra):
    return spectra.__dict__

def save_list_c(self):
    file_name = self.file_list.get()
    with open(file_name, "w") as file:
        json.dump([spectra_to_dict(spectra) for spectra in list_of_spectra_c], file, indent=4)
    messagebox.showinfo('Info', 'Done!')

def dict_to_spectra(d):
    return Spectra(**d)


def read_list_c(self):
    file_name = self.open_spectra_list_file.get()
    #list_of_spectra_c.clear()
    with open(file_name, "r") as file:
        list_of_spectra_c_from_file = json.load(file)
        for spectra in list_of_spectra_c_from_file:
            list_of_spectra_c.append(Spectra(spectra.get("name"), spectra.get("counts"), spectra.get("channel"), spectra.get("energy"), spectra.get("material"), spectra.get("crystal_orientation"), spectra.get("axial_channel"), spectra.get("tilt_angle"), spectra.get("filter"), spectra.get("background"), spectra.get("final_yield")))
    show_plot_list_c(self)
    messagebox.showinfo('Info', 'Done!')
    

def open_folder_2(self, root):
    filetypes_ = (('json files','*.json'), ('All files', '*.*'))
    file = inputFileName = askopenfile(parent=root, title='Open a file', initialdir=self.open_spectra_list_file_loc, filetypes=filetypes_)
    self.open_spectra_list_file.set(file.name) 
    self.open_spectra_list_file_loc = self.open_spectra_list_file

