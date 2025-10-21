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

# loadData,
# self.selectedSpectraC, listSpectraCountS, showPlotCounts, showRange

class PeakCounts:
    def __init__(self, _id_, _peak_, _material_, _axis_, _plane_, _tiltangle_, _counts_, _error_):
        self.id = _id_
        self.peak = _peak_
        self.material = _material_
        self.axis = _axis_
        self.plane = _plane_
        self.tiltangle = _tiltangle_
        self.counts = _counts_
        self.error = _error_
    def __str__(self):
        return {
            "id": self.id,
            "peak": self.peak,
            "material": self.material,
            "axis": self.axis,
            "plane": self.plane,
            "tiltangle": self.tiltangle,
            "counts": self.counts,
            "error": self.error,
        }




def loadData(self):
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
    B = []
    with open(fdRes.name, 'r', encoding="utf-8") as file:
        infoRead = file.read()
        file.seek(0)
        for line in file:
            line = line.strip()
            if not line.startswith("#"): 
                x,yr = line.split()
                E.append(float(x))
                Yr.append(float(yr))
                B.append(float(yr))
                Y.append(float(yr))
        material = re.findall('#    Material: (\S*)', infoRead)
        tiltAngle = re.findall('#    Tilt angle: (\S*)', infoRead)
        axis = re.findall('#    Axis: (\S*)', infoRead)
        plane = re.findall('#    Plane: (\S*)', infoRead)
        dose = re.search('#    Dose: (\S*)', infoRead)
        if dose:
            dose = dose.group(1)
        file.close()
    number = len(self.listSpectraCount) + 1
    Info = {'material': material[0], 'tiltAngle': tiltAngle[0], 'axis': axis[0], 'plane': plane[0], 'dose': dose}
    new_spectrum.append(E)
    new_spectrum.append(Y)
    new_spectrum.append(Info)
    new_spectrum.append(number)
    self.listSpectraCount.append(new_spectrum)
    newTree = self.treeC.insert('', 'end', text=str(number), values=(str(number), str(tiltAngle[0]), str(plane[0]), str(axis[0]), str(material[0])))
    return E,Y,Info

def addSpectra(self, x, y, spectrumInfo):
    new_spectrum = []
    number = len(self.listSpectraCount) + 1
    Info = spectrumInfo
    new_spectrum.append(x)
    new_spectrum.append(y)
    new_spectrum.append(Info)
    new_spectrum.append(number)
    self.listSpectraCount.append(new_spectrum)
    newTree = self.treeC.insert('', 'end', text=str(number), values=(str(number), str(Info['tiltAngle']), str(Info['plane']), str(Info['axis']), str(Info['material'])))
    print('done')

def plotSpectraCount(self):
    fig=self.__fig__
    fig.clear()
    self.listSpectraCountS.clear()
    self.__plt__=fig.add_subplot(111)
    for item in self.selectedSpectraC:
        spectrum=self.treeC.item(item)['values']
        idSpectrum = int(spectrum[0])
        for s in self.listSpectraCount:
            if s[3] == idSpectrum:
                self.listSpectraCountS.append(s)
                self.__plt__.set_yscale('log')
                #label=s[2]['tiltAngle']
                
                self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
                self.legend=self.__plt__.legend()
    fig.canvas.draw()

def showRange(self):
    fig=self.__fig__
    a_min = float(self.emin_as.get())
    b_max = float(self.emax_as.get())
    #self.__plt__=fig.add_subplot(111)
    A = []
    B = []
    A.append(a_min)
    B.append(b_max)
    A.append(a_min)
    B.append(b_max)
    Ym = []
    Ym.append(0)
    Ym.append(1e6)
    self.__plt__.plot(A, Ym, '-', color='orange')
    self.__plt__.plot(B, Ym, '-', color='orange')
    self.__plt__.fill_betweenx(Ym, A, B, color='orange', alpha=0.5)
    fig.canvas.draw()


def showPlotCounts(self):
    fig=self.__fig__
    fig.clear()
    self.__plt__=fig.add_subplot(111)
    self.Counts.clear()
    a_min = float(self.emin_as.get())
    b_max = float(self.emax_as.get())
    i1 = 0
    i1_error = 0
    li = []
    ang = []
    for s in self.listSpectraCountS:
        i1 = i1 + 1
        as_counts = 0
        err = 0
        for o in range(0, len(s[0][:])):
            if float(s[0][o]) > a_min and float(s[0][o]) < b_max:
                as_counts = as_counts+float(s[1][o])
                #err = sqrt(float(s[1][o]) + float(s[1][o]))
        PeakCounts.id = i1
        PeakCounts.peak = a_min
        PeakCounts.material = s[2]['material']
        PeakCounts.axis = s[2]['axis']
        PeakCounts.plane = s[2]['plane']
        PeakCounts.tiltangle = s[2]['tiltAngle']
        PeakCounts.counts = as_counts
        PeakCounts.error = err
        self.Counts.append(PeakCounts(PeakCounts.id, PeakCounts.peak, PeakCounts.material, PeakCounts.axis, PeakCounts.plane, PeakCounts.tiltangle, PeakCounts.counts, PeakCounts.error))
        self.N.append(i1)
        self.C.append(as_counts)
        self.Error.append(err)
        ang.append(float(s[2]['tiltAngle']))
    numberScan = len(self.PIXEcount) + 1
    nameScan = str(str(a_min) + ' - ' + str(b_max))
    li.append(numberScan)
    li.append(self.N)
    li.append(ang)
    li.append(self.C)
    li.append(nameScan)
    self.PIXEcount.append(li)
    le = len(self.treeAS.get_children()) + 1
    self.treeAS.insert('', 'end', text=str(le), values=(str(numberScan)))
    self.__plt__.bar(ang, self.C, yerr=self.Error)
    fig.canvas.draw()
    return self.N, self.C, self.Error

def showPlotCountsN(self):
    fig=self.__fig__
    fig.clear()
    self.__plt__=fig.add_subplot(111)
    self.Counts.clear()
    a_min = float(self.emin_as.get())
    b_max = float(self.emax_as.get())
    i1 = 1
    i1_error = 0
    for s in self.listSpectraCountS:
        R_counts = 0
        if s[2]['tiltAngle'] == 'R':
            for o in range(0, len(s[0][:])):
                if float(s[0][o]) > a_min and float(s[0][o]) < b_max:
                    R_counts = R_counts+float(s[1][o])
        else:
            R_counts = 1
    for s in self.listSpectraCountS:
        i1 = i1 + 1
        as_counts = 0
        err = 0
        for o in range(0, len(s[0][:])):
            if float(s[0][o]) > a_min and float(s[0][o]) < b_max:
                as_counts = as_counts+float(s[1][o])
                #err = sqrt(float(s[1][o]) + float(s[1][o]))
        as_counts = as_counts / R_counts
        PeakCounts.id = i1
        PeakCounts.peak = a_min
        PeakCounts.material = s[2]['material']
        PeakCounts.axis = s[2]['axis']
        PeakCounts.plane = s[2]['plane']
        PeakCounts.tiltangle = s[2]['tiltAngle']
        PeakCounts.counts = as_counts
        PeakCounts.error = err
        self.Counts.append(PeakCounts(PeakCounts.id, PeakCounts.peak, PeakCounts.material, PeakCounts.axis, PeakCounts.plane, PeakCounts.tiltangle, PeakCounts.counts, PeakCounts.error))
        self.N.append(i1)
        self.C.append(as_counts)
        self.Error.append(err)
    nameScan = str(str(a_min) + ' - ' + str(b_max))
    li.append(numberScan)
    li.append(self.N)
    li.append(ang)
    li.append(self.C)
    li.append(nameScan)
    self.PIXEcount.append(li)
    le = len(self.treeAS.get_children()) + 1
    self.treeAS.insert('', 'end', text=str(le), values=(str(numberScan)))
    self.__plt__.bar(ang, self.C, yerr=self.Error)
    fig.canvas.draw()
    return self.N, self.C, self.Error

def saveCount(self):
    from datetime import datetime
    idir=iniSet.getValue('default','inipath')
    if idir == '' : idir=self.currDir
    dRes = fd.asksaveasfile(title='Save a file', initialdir=idir,                                                                
                            filetypes=[("Background free TXT","*.txt"),
                                ("All files","*.*")])
    if dRes is None:
        print('no file selected')
        return
    n=self.N
    c=self.C
    e=self.Error
    #ang = self.Counts[0].tiltangle
    newTree = self.treeAS.insert('', 'end', values=(str(PeakCounts.id), str(PeakCounts.peak)))
    fid=open(dRes.name,'w')
    #fid.write('#    DateTime: ' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')+'\n')
    #fid.write('#    Material: ' + self.spectrum.Material + '\n')
    #fid.write('#    Axis: '  + self.spectrum.Axis + '\n')
    #fid.write('#    Plane: ' + self.spectrum.Plane + '\n')
    #fid.write('#    Tilt angle: ' + self.spectrum.TiltAngle + '\n')
    #fid.write()
    for i in range(0,len(self.N[:])):
            #ftok=f'{n[i]} {self.Counts[i].tiltangle} {c[i]:12.6f} {e[i]:12.6f}\n'
            ftok=f'{self.Counts[i].tiltangle} {c[i]:12.6f}\n'
            fid.write(ftok)                                             
    fid.close()
