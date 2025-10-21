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


from math import sin, pi
from sys import exc_info
import matplotlib.pyplot as plt
import re
import numpy as np
from scipy import stats

from tkinter import filedialog

def extract_values(file_path, value):
    regex_string = regex[value]
    value = re.findall(regex_string, open(file_path,  encoding="utf8", errors='ignore').read(), re.DOTALL)
    float_values = []
    for element in value: 
        new_element = float(element)
        float_values.append(new_element)
    return float_values


regex = {
    'c1' : f"(\d*)*\d*",
    'c2' : f"\d*\s*(\d*)\n"
}

def openRandom(self):
    fdRes = fd.askopenfile(title='Open a file')
    fileName=fdRes.name
    if re.search('\.dat$',fileName):
        self.RBSrandom = extract_values(fileName, 'c2')
        fig=self.__fig__
        fig.clear()
        self.__plt__=fig.add_subplot(111)
        self.__plt__.plot(self.RBSrandom, '.', label='Random', markersize=0.5, linestyle='none')
        self.__plt__.set_xlabel('Channel')
        self.__plt__.set_ylabel('Yield')
        fig.canvas.draw()

def openRBSC(self):
    fdRes = fd.askopenfile(title='Open a file')
    fileName=fdRes.name
    Channel = []
    RBS_yield = []
    if re.search('\.dat$',fileName):
        self.RBSscan = extract_values(fileName, 'c2')
        dose = float(self.doseRBS.get())/float(self.doseRBSa.get())
        for a in range(0, len(self.RBSscan)-1):
            Channel.append(a)
            RBS_yield.append(float(self.RBSscan[a])*dose)
        angle = self.angleRBS.get()
        fig=self.__fig__
        #fig.clear()
        #self.__plt__=fig.add_subplot(111)
        self.RBSlist.append((Channel, RBS_yield, angle))
        self.__plt__.plot(RBS_yield, '.', label=angle, markersize=0.5, linestyle='none')
        self.__plt__.set_xlabel('Channel')
        self.__plt__.set_ylabel('Yield')
        fig.canvas.draw()

def showRBSlim(self): # tu dodano przekazywanie Info do zapisu
    fig=self.__fig__
    X1 = []
    X2 = []
    Y = []
    x1 = float(self.channelRBSmin.get())
    x2 = float(self.channelRBSmax.get())
    Y.append(0)
    Y.append(2e4)
    X1.append(x1)
    X1.append(x1)
    X2.append(x2)
    X2.append(x2)
    self.__plt__.plot(X1, Y, label='min', color='yellow')
    self.__plt__.plot(X2, Y, label='max', color='yellow')
    self.__plt__.fill_betweenx(Y, X1, X2, color='yellow', alpha=0.2)
    fig.canvas.draw()

    
def RBSas(self):
    x1 = float(self.channelRBSmin.get())
    x2 = float(self.channelRBSmax.get())
    nameRBS = self.nameRBS.get()
    self.RBScounts.clear()
    cr = 0
    for i in range(len(self.RBSrandom)):
        if x1 <= i <= x2:
            cr = cr + float(self.RBSrandom[i])
    for item in self.RBSlist:
        c = 0
        for i in range(len(item[0])):
            if x1 <= i <= x2:
                c = c + float(item[1][i])
        c = c /cr
        self.RBScounts.append([c, float(item[2])])
    #self.RBScountL[nameRBS] = self.RBScounts
    le = len(self.treeAS.get_children()) + 1
    self.treeAS.insert('', 'end', text=str(le), values=(str(nameRBS)))
    Ang = []
    co = []
    for rbs in self.RBScounts:
        Ang.append(str(rbs[1]))
        co.append(rbs[0])
    fig=self.__fig__
    fig.clear()
    self.__plt__=fig.add_subplot(111)
    #print(Ang)
    #print(co)
    self.__plt__.bar(Ang, co)
    self.__plt__.set_xlabel('Angle')
    self.__plt__.set_ylabel('Counts')
    fig.canvas.draw()

def RBSsave(self):
    newTree = self.treeAS.insert('', 'end', values=(str(PeakCounts.id), str(PeakCounts.peak)))
    dRes = fd.asksaveasfile(title='Save a file')
    if dRes is None:
        print('no file selected')
        return
    Counts = []
    Angle = []
    for rbs in self.RBScounts:
        Angle.append(rbs[1])
        Counts.append(rbs[0])
    fid=open(dRes.name,'w')
    fid.write('#    Angle [deg]    Counts   \n')
    for i in range(0,len(Counts[:])):
        ftok=f'{Angle[i]:9.6f} {Counts[i]:12.6f}\n'
        fid.write(ftok)                                             
    fid.close()
