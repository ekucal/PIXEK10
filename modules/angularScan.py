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

import  modules.inisettings as iniSet


# self.selectedScan, self.selectedScans

def addDataAS(self):
    print('addDataAS')

'''
def showplotAS(self):
    fig=self.__fig__
    fig.clear()
    self.__plt__=fig.add_subplot(111)
    for item in self.selectedScan:
        spectrum=self.treeAS.item(item)['values']
        idSpectrum = int(spectrum[0])
        for s in self.PIXEcount:
            if s[0] == idSpectrum:
                self.__plt__.plot(s[2], s[1], 's', markersize=10,  linestyle='None') #, legend=s[4])
    #self.__plt__.legend(loc='lower right', fontsize=15)
    #self.__plt__.set_ylabel('Normalised Yield', fontsize=20, fontweight='bold')
    #self.__plt__.set_xlabel('Tilt angle [degrees]', fontsize=20, fontweight='bold')
    #self.__plt__.tick_params(axis='y', labelsize=15)
    #self.__plt__.tick_params(axis='x', labelsize=15)
    #self.__plt__.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    #self.__plt__.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    #self.__plt__.minorticks_on()
    fig.canvas.draw()
'''

def showplotAS(self):
    fig=self.__fig__
    fig.clear()
    plt=fig.add_subplot(111)
    for item in self.selectedScan:
        spectrum=self.treeAS.item(item)['values']
        idSpectrum = int(spectrum[0])
        for s in self.PIXEcount:
            if s[0] == idSpectrum:
                plt.plot(s[2], s[1], 's', markersize=10,  linestyle='None') #, legend=s[4])
    #plt.legend(loc='lower right', fontsize=15)
    plt.set_ylabel('Normalised Yield', fontsize=20, fontweight='bold')
    plt.set_xlabel('Tilt angle [degrees]', fontsize=20, fontweight='bold')
    plt.tick_params(axis='y', labelsize=15)
    plt.tick_params(axis='x', labelsize=15)
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    plt.minorticks_on()
    fig.canvas.draw()

def plotAS(self): 
    fig=self.__fig__
    fig.clear()
    self.__plt__=fig.add_subplot(111)
    for item in self.PIXEcount:
        spectrum=self.treeAS.item(item)['values']
        idSpectrum = int(spectrum[0])
        if item[0] == idSpectrum:
            self.__plt__.plot(item[0], item[1], '.', markersize=0.5, color='red', linestyle='None')
    self.__plt__.legend(loc='lower right', fontsize=15)
    self.__plt__.ylabel('Normalised Yield', fontsize=20, fontweight='bold')
    self.__plt__.xlabel('Tilt angle [degrees]', fontsize=20, fontweight='bold')
    self.__plt__.tick_params(axis='y', labelsize=15)
    self.__plt__.tick_params(axis='x', labelsize=15)
    self.__plt__.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
    self.__plt__.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
    self.__plt__.minorticks_on()
    fig.canvas.draw()

# lepszy zapis do angular scanu