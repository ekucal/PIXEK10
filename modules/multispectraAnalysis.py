import numpy as np
from math import sin, pi, exp, sqrt, log

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler

import tkinter as tk
from tkinter import IntVar,StringVar
from tkinter import filedialog as fd
from tkinter import Toplevel
from tkinter import messagebox

from scipy.signal import find_peaks
from skued import gaussian, spectrum_colors
from scipy.signal import savgol_filter
from skued import baseline_dwt

import os,re
import pywt
import threading
import csv


# import  modules.multispectraAnalysis as openFolder, removeNoise, removeBackground, saveData
import  modules.inisettings as iniSet


def openFolder(self):
    new_spectrum = []
    currDir=''
    idir=iniSet.getValue('default','inipath')
    if idir == '' : 
        idir=self.currDir
    fdRes = fd.askopenfile( title='Open a file', initialdir=idir, filetypes=[("Savad spectra", "*.txt")])
    if not fdRes :
        self.spectrum.clear()
        print('WARNING: no data loaded')    
    E = []
    Y = []
    Yr = []
    with open(fdRes.name, 'r', encoding="utf-8") as file:
        infoRead = file.read()
        file.seek(0)
        for line in file:
            line = line.strip()
            if not line.startswith("#"): 
                x,y, yr = line.split()
                E.append(float(x))
                Y.append(float(y))
                Yr.append(float(yr))
        material = re.findall('#    Material: (\S*)', infoRead)
        tiltAngle = re.findall('#    Tilt angle: (\S*)', infoRead)
        axis = re.findall('#    Axis: (\S*)', infoRead)
        plane = re.findall('#    Plane: (\S*)', infoRead)
        dose = re.search('#    Dose: (\S*)', infoRead)
        if dose:
            dose = dose.group(1)
        file.close()
    number = len(self.listSpectraO) + 1
    Info = {'material': material[0], 'tiltAngle': tiltAngle[0], 'axis': axis[0], 'plane': plane[0], 'dose': dose}
    new_spectrum.append(E)
    new_spectrum.append(Y)
    new_spectrum.append(Info)
    new_spectrum.append(number)
    new_spectrum.append(yr)
    self.listSpectraO.append(new_spectrum)
    newTree = self.tree.insert('', 'end', text=str(number), values=(str(number), str(tiltAngle[0]), str(plane[0]), str(axis[0]), str(material[0])))
    return E,Y,Info

def addCalibratedpectrum(self):
    new_spectrum = []
    number = len(self.listSpectraO) + 1
    Info = {'material': self.spectrum.Material, 'tiltAngle': self.spectrum.TiltAngle, 'axis': self.spectrum.AxialChannel, 'plane': self.spectrum.Plane, 'dose': self.spectrum.Dose}
    new_spectrum.append(self.spectrum.Xc)
    new_spectrum.append(self.spectrum.Yi)
    new_spectrum.append(Info)
    new_spectrum.append(number)
    self.listSpectraO.append(new_spectrum)
    newTree = self.tree.insert('', 'end', text=str(number), values=(str(number), str(self.spectrum.TiltAngle), str(self.spectrum.Plane), str(self.spectrum.AxialChannel), str(self.spectrum.Material)))
    print('done')






    


    

