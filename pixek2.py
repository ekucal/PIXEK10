#!/home/skrobask/bin/pyt3

#########################################
###      PIXE.spe data analysis       ###
#########################################

import numpy as np
from math import sin, pi, exp, sqrt, log
from math import ceil
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

import ttkbootstrap as ttk
from ttkbootstrap.style import Bootstyle
from tkinter.filedialog import askdirectory
from ttkbootstrap.dialogs import Messagebox
from ttkbootstrap.constants import *
from tkinter.scrolledtext import ScrolledText
from tooltip import ToolTip

import os,re
import pywt
import threading
import csv

#########################################

from modules.inputData import SpectrumData
from modules.calibrationData import XaxisMapping
from modules.filteringData import Filter
from modules.labeledInput import LabeledInput
from modules.clipPlotWindow import clipPlotWindow
from modules.showCWT import showCWT
import  modules.inisettings as iniSet
from  modules.multispectraAnalysis import openFolder, addCalibratedpectrum
from modules.calCounts import loadData, plotSpectraCount,  showPlotCounts, showRange, PeakCounts, saveCount, addSpectra, showPlotCountsN
from modules.angularScan import plotAS, addDataAS, showplotAS
from modules.RBSC import openRandom, openRBSC, showRBSlim, RBSas, RBSsave


# NUMERICAL VALUES SEPRATOR
SEP=';'
DEFRANGE='0 '+SEP+' -1'
# REAL NUMBER EXPRESSION
RN="[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
# TEST FOR INT/DOUBLE NUMBER CONDITIONS
CNI='^\s*[0-9]+\s*'+SEP+'\s*[-]?[0-9]+\s*$'
CND='^\s*'+RN+'\s*;\s*'+RN+'\s*$'

class AppGUI():

    spectrum=SpectrumData()
    xmapping=XaxisMapping()
    filtr=Filter()

    # List of spectra in diffrents steps of analysis

    listSpectraO = []
    listSpectraA = []
    listSpectraC = []
    listSpectraR = []
    listSpectraF = []
    listSpectraB = []

    listSpectraCount = []
    listSpectraCountS = []

    selectedSpectrum = None
    selectedSpectra = []
    selectedSpectraC = []

    selectedScan = []
    selectedScans = []

    Counts = []
    PIXEcount = []
    listCounts = []

    RBSrandom = []
    RBSscan = []
    RBSlist = []
    RBScounts = []
    RBScountL = {}

    list_of_spectra = []
    list_of_spectra_c = []
    list_of_spectra_f = []
    list_of_spectra_b = []
   
    a_calibration = 1
    b_calibration = 0

    a_calibration_M = 1
    b_calibration_M = 0
    
    N = []
    C = []
    Error = []

    currDir=''
    buttons={}
    edits={}
    comboBoxes={}
    
    actParams={}
    
    __fig__=None
    __plt__=None
    
    grid=True
    legendSwitch=True

    dataInfo={}

    ### Open file (self.spectrum.Xi, self.spectrum.Yi)

    def openFile(self):
        idir=iniSet.getValue('default','inipath')
        if idir == '' : idir=self.currDir
        fdRes = fd.askopenfile( title='Open a file', initialdir=idir, filetypes=self.spectrum.fileTypes)
        if not fdRes :
            self.spectrum.clear()
            print('WARNING: no data loaded')
            self.open_spectra_file.set("NO DATA")
            return
        iniSet.setValue('default','inipath',os.path.dirname(fdRes.name))
        if self.spectrum.loadSpectrum(fdRes.name):
            self.open_spectra_file.set(self.spectrum.ifileName.split('/')[-1])            
            self.plotData([self.spectrum.Xi],[self.spectrum.Yi],[self.spectrum.plotScheme],
                          dlabels=['Raw data'],
                          plt_kwargs=[{}],
                          xlabel="Channels",ylabel=self.spectrum.ylabel)
            dataSize=str(self.spectrum.dataSize())
            self.spectra_size.set("Number of points: "+dataSize)
            self.spectra_range.set("0  "+SEP+' '+dataSize)
            
            self.spectrum.Material = self.material.get()
            self.spectrum.CrystalOrientation = self.material.get()
            self.spectrum.AxialChannel = self.axis.get()
            self.spectrum.Plane = self.plane.get()
            self.spectrum.TiltAngle = self.angle.get()
            self.spectrum.Dose = self.dose_e.get()
            self.spectrum.DosePlot = self.dose.get()
        else:
            self.open_spectra_file.set("ERROR")            

    def plotData(self,X:list[np.ndarray],Y:list[np.ndarray],
                 plt_args=list[str],dlabels=list[str],
                 plt_kwargs=None,xlabel=None,ylabel=None):
        #plt_args=plt_args or {}
        #plt_kwargs=plt_kwargs or {}
        fig=self.__fig__
        fig.clear()
        self.__plt__=fig.add_subplot(111)
        for i in range(0,len(X)):  
            #plt_kw=(plt_kwargs[i])
            if self.Scale.get()==0:
                self.__plt__.plot(X[i],Y[i]) #,plt_args[i],**plt_kw)
            if self.Scale.get()==1:
                self.__plt__.set_yscale('log')
                self.__plt__.legend(loc='upper right')
                self.__plt__.plot(X[i],Y[i]) #,plt_args[i],**plt_kw)
        self.__plt__.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
        self.__plt__.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
        self.__plt__.minorticks_on()
        self.__plt__.set_ylim(1,1.5*self.spectrum.Yi.max())
        if dlabels is not None:
            self.legend=self.__plt__.legend(dlabels)
        if xlabel is not None:
            self.__plt__.set_xlabel(xlabel)
        if ylabel is not None:
            self.__plt__.set_ylabel(ylabel)
        fig.canvas.draw()
    
    ### Calibration of X axis (self.spectrum.Xc)

    def dataMapped(self):      
        #if self.dataReduce()==False:
         #   return False       
        if re.search(CND,self.mapping_axb.get()) is None:
           print('ERROR: calibration a, b prms wrong format')
           tk.messagebox.showerror('ERROR',message='ERROR: calibration a, b prms wrong format')
           self.mapping_axb.set("1 ; 0")
           return False

        method=self.spectra_mapping.get()
        a,b=1,0
        if method=='Ignore':
            a,b=1,0
            self.a_calibration = a
            self.b_calibration = b
        if method=='E=ax+b':
            a,b=self.mapping_axb.get().split(SEP)
            self.a_calibration = a
            self.b_calibration = b
        if method=='Peaks':
            P1e = float(self.peak1energy.get())
            P2e = float(self.peak2energy.get())
            P1c = float(self.peak1channel.get())
            P2c = float(self.peak2channel.get())
            a = (P2e-P1e)/(P2c-P1c)
            b = P2e - a * P2c
            self.a_calibration = a
            self.b_calibration = b
        self.spectrum.Xc=self.xmapping(self.spectrum.Xi,float(a),float(b),method)
        self.xmapping.xlabel='Channels' if method=='Ignore' else 'Energy [keV]' 
        self.dataInfo['cMethod']=method
        self.dataInfo['axb']=self.mapping_axb.get()
        return True
        
    def replotMappedData(self):
        method=self.spectra_mapping.get()
        if self.xmapping.isIgnore(method):
            return True
        if self.dataMapped():            
            self.plotData([self.spectrum.Xc],[self.spectrum.Yi],['.'],dlabels=['Raw spectrum'],                          
                          plt_kwargs=[{'markersize' : 1, 'color' : 'black'}],
                          xlabel=self.xmapping.xlabel,ylabel=self.spectrum.ylabel)
        addCalibratedpectrum(self)


    ### Reduction of data (self.spectrum.Xr, self.spectrum.Yr)

    def replotInputData(self):
        if self.dataReduce():            
            self.plotData([self.spectrum.Xr],[self.spectrum.Yr],['.'],dlabels=['Raw spectrum'],plt_kwargs=[{'markersize' : 1, 'color' : 'black'}],
                          xlabel="Energy [keV]",ylabel=self.spectrum.ylabel)     

    def dataReduce(self):
        if self.spectrum.empty():
            print('WARNING: spectrum data no loaded')
            tk.messagebox.showwarning('WARNING',message='spectrum data no loaded')
            return False

        spectra_range=self.spectra_range.get()
        if re.search(CNI,spectra_range) is None:
            print('ERROR: spectra range wrong format')
            tk.messagebox.showerror('ERROR',message='spectra range wrong format')
            self.spectra_range.set(DEFRANGE)
            return False

        fromTo=spectra_range.split(SEP)
        # calibration method
        from_=int((int(fromTo[0])-self.b_calibration)//self.a_calibration)
        to_  =int((int(fromTo[1])-self.b_calibration)//self.a_calibration)
        spSize=self.spectrum.dataSize()

        if to_<0:
            to_+=spSize

        if from_>to_:
            print('ERROR: spectra range wrong values: from>to')
            tk.messagebox.showerror('ERROR',message='ERROR: spectra range wrong values: from>to')
            self.spectra_range.set(DEFRANGE)
            return False

        if from_>spSize or to_>spSize:
            print('ERROR: value(s) out of range')
            tk.messagebox.showerror('ERROR',message='ERROR: values out of range')
            self.spectra_range.set(DEFRANGE)
            return False

        method=self.spectra_reduction.get()            
        self.spectrum.Xr,self.spectrum.Yr=self.getReducedData(from_,to_,method)
        self.spectrum.Y = self.spectrum.Yr
        self.dataInfo['dataRange']=spectra_range
        self.dataInfo['reduceMethod']=method
        return True
    
    def getReducedData(self,from_,to_,rtype):
        self._from=from_
        self._to=to_
        self.rtype=rtype
        if rtype=='none':
            self.spectrum.ylabel='Yield'
            return self.spectrum.Xc[from_:to_],self.spectrum.Yi[from_:to_]
        if rtype=='sqrt':
            nYi=self.spectrum.Yi[from_:to_]
            ignoreInd=np.where(nYi<0)
            nYi[ignoreInd]=0
            self.ylabel='sqrt(counts)'
            return self.spectrum.Xc[from_:to_],nYi**0.5
        if rtype=='log10':
            nYi=self.spectrum.Yi[from_:to_]
            ignoreInd=np.where(nYi<1)
            nYi[ignoreInd]=1
            self.ylabel='log10(counts)'
            return self.spectrum.Xc[from_:to_],np.log10(nYi)

        print('ERROR: not recognized reduction type')
        return None,None

    ### Filtering of data (self.spectrum.Xf, self.spectrum.Yf)

    

    def filterData(self):                  
        if self.dataMapped()==False:
            return False
                                
        filterKernelSize=self.filterKernelSize.get()     
        
        if filterKernelSize<1 or  filterKernelSize>self.spectrum.dataSize():
            print('ERROR: kernel size out of range')
            tk.messagebox.showerror('ERROR',message='ERROR: kernel size out of range')
            self.filterKernelSize.set(1)
            return False
        
        method=self.filterMethod.get()
        
        if method==self.filtr.methods[1]:
            filterPolyOrder=0
            filterDerivative=0
        else:
            filterPolyOrder=self.filterPolyOrder.get()
            filterDerivative=self.filterDeriv.get()
                                  
        if self.spectrum.Yr is None:
            self.spectrum.Yf=self.filtr(self.spectrum.Yi,filterKernelSize,method,filterPolyOrder,filterDerivative)  
        else:
            self.spectrum.Yf=self.filtr(self.spectrum.Yr,filterKernelSize,method,filterPolyOrder,filterDerivative)     
        self.spectrum.Y = self.spectrum.Yf
        self.dataInfo['filterMethod']=method
        self.dataInfo['filterSize']=filterKernelSize
        self.dataInfo['filterPoly']=filterPolyOrder
        self.dataInfo['filterDeriv']=filterDerivative
        
        return True



    def replotFilteredData(self):        
        if self.filtr.isIgnore(self.filterMethod.get()):
            return  False
        
        if self.filterData():            
            coeff= self.spectrum.Yr.max()/self.spectrum.Yf.max() if  self.filterDeriv.get()>0 else 1
            if  self.filterDeriv.get()==2:
                coeff*=-1
            
            if self.filterShowRawData.get()==1:                                                      
                self.plotData([self.spectrum.Xr,self.spectrum.Xr],[self.spectrum.Yr,self.spectrum.Yr*coeff],[self.spectrum.plotScheme,self.filtr.plotScheme],
                              dlabels=['raw','filtered'],plt_kwargs=[{'alpha' : 0.5,'markersize' : 3}, {}],
                              xlabel=self.xmapping.xlabel,ylabel="values")
            else:
                self.plotData([self.spectrum.Xr],[self.spectrum.Yf],[self.filtr.plotScheme],dlabels=['Filtered'],                              
                              plt_kwargs=[{}],xlabel=self.xmapping.xlabel,ylabel=self.spectrum.ylabel)


    ### Background removal (self.spectrum.Xb, self.spectrum.Yb)

    def plotRemoveBackground(self):
        #if self.filterData()==False:
         #   return
        wlevel=self.waveletLevels.get()
        wname=self.waveletName.get()
        self.baseline=baseline_dwt(self.spectrum.Y,level=wlevel,wavelet=wname,max_iter=100)
        if self.rtype=='none':
            self.spectrum.Yb=self.spectrum.Y-self.baseline
            y = self.spectrum.Yr
            b = self.baseline
            if self.spectrum.Xr is None:
                x = self.spectrum.Xc
            else:
                x=self.spectrum.Xr
        if self.rtype=='sqrt':
            self.spectrum.Yb=(self.spectrum.Y-self.baseline)**2
            y = self.spectrum.Yr**2
            b = self.baseline**2
            if self.spectrum.Xr is None:
                x = self.spectrum.Xc
            else:
                x=self.spectrum.Xr
        if self.rtype=='log10':
            self.spectrum.Yb=10**(self.spectrum.Y-self.baseline)
            y = 10**self.spectrum.Yr
            b = 10**self.baseline
            if self.spectrum.Xr is None:
                x = self.spectrum.Xc
            else:
                x=self.spectrum.Xr
        self.spectrum.Xb=x
        self.spectrum.Yb=y
        self.plotData([x,x,x],
                     [y,b,self.spectrum.Yb],
                     [self.filtr.plotScheme,'-k','-m'],
                     dlabels=['Raw data','Background','Background free'],
                     plt_kwargs=[{'alpha' : 0.5,'linewidth' : 0.75},{'alpha' : 0.5,'linewidth' : 0.75},{}],
                     xlabel=self.xmapping.xlabel,ylabel=self.spectrum.ylabel)
        self.dataInfo['waveletname']=wname
        self.dataInfo['waveletlevel']=wlevel   

    #### Save

    ######################################################################################################
    ######################################################################################################
    ######## Multispectra analysis #######################################################################
    ######################################################################################################

    # listSpectraO = []
    # listSpectraA = [] -> tylko wybrane z listy TreeC
    # listSpectraC = [] -> kalibracja, tylko z parametrami a, b, w innym przypadku trzeba przejsc kalibracje w add spectra
    # listSpectraR = []
    # listSpectraF = []
    # listSpectraB = []

    ### Open

    def showSpectra(self):
        fig=self.__fig__
        fig.clear()
        self.listSpectraA.clear()
        self.__plt__=fig.add_subplot(111)
        for item in self.selectedSpectra:
            spectrum=self.tree.item(item)['values']
            idSpectrum = int(spectrum[0])
            for s in self.listSpectraO:
                if s[3] == idSpectrum:
                    self.listSpectraA.append(s)
                    if self.Scale.get()==0:
                        self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
                    if self.Scale.get()==1:
                        self.__plt__.set_yscale('log')
                        self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
        self.__plt__.set_xlabel('Energy  [keV]')
        self.__plt__.set_ylabel('Yield')
        #self.legend=self.__plt__.legend()
        fig.canvas.draw()

    ### Calbration

    def replotMappedDataM(self):
        fig=self.__fig__
        fig.clear()
        self.__plt__=fig.add_subplot(111)
        a,b=self.mapping_axbM.get().split(SEP)
        a_calibration_M = a
        b_calibration_M = b
        for s in self.listSpectraA:
            s[0]=self.xmapping(s[0],float(a),float(b),'E=ax+b')
            if self.Scale.get()==0:
                self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
            if self.Scale.get()==1:
                self.__plt__.set_yscale('log')
                self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
            self.listSpectraC.append(s)
        fig.canvas.draw()

    ### Reduction

    def replotM(self):
        fig=self.__fig__
        fig.clear()
        self.dataReduceM()
        self.__plt__=fig.add_subplot(111)
        for s in self.listSpectraR:
            if self.Scale.get()==0:
                self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
            if self.Scale.get()==1:
                self.__plt__.set_yscale('log')
                self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)
        self.__plt__.set_xlabel('Energy  [keV]')
        self.__plt__.set_ylabel('Yield')
        fig.canvas.draw()

    def dataReduceM(self):
        spectra_ran=self.spectra_range_M.get()
        fromTo=spectra_ran.split(SEP)
        method=self.spectra_reduction_M.get()   
        if self.listSpectraC:
            for s in self.listSpectraC: 
                sx = s[0]
                sy = s[1]
                from__ = int(fromTo[0])
                to__ = int(fromTo[1])
                sx = np.array(s[0])
                sy = np.array(s[1])
                for i in range(len(sx)):
                    if sx[i] < from__:
                        from_ = i
                    if sx[i] < to__:
                        to_ = i
                s[0],s[1]=self.getReducedDataM(from_,to_,method, s[0], s[1])
                self.listSpectraR.append(s)
        else:
            for s in self.listSpectraA:          
                sx = s[0]
                sy = s[1]
                from__ = int(fromTo[0])
                to__ = int(fromTo[1])
                sx = np.array(s[0])
                sy = np.array(s[1])
                for i in range(len(sx)):
                    if sx[i] < from__:
                        from_ = i
                    if sx[i] < to__:
                        to_ = i
                s[0],s[1]=self.getReducedDataM(from_,to_,method, s[0], s[1])
                self.listSpectraR.append(s)
             

    
    def getReducedDataM(self,from_,to_,rtype, sx, sy):
        self._from=from_
        self._to=to_
        self.rtypeM=rtype
        if rtype=='none':
            return sx[from_:to_],sy[from_:to_]
        if rtype=='sqrt':
            nYi=np.array(sy[from_:to_])
            ignoreInd=np.where(nYi<0)
            nYi[ignoreInd]=0
            return sx[from_:to_],nYi**0.5
        if rtype=='log10':
            nYi=np.array(sy[from_:to_])
            ignoreInd=np.where(nYi<1)
            nYi[ignoreInd]=1
            return sx[from_:to_],np.log10(nYi)
        print('ERROR: not recognized reduction type')
        return None,None

    
        
    ### Filtering


    ### Background removal

    def plotRemoveBackgroundM(self):
        fig=self.__fig__
        fig.clear()
        wlevel=self.waveletLevelsM.get()
        wname=self.waveletNameM.get()
        self.__plt__=fig.add_subplot(111)
        for s in self.listSpectraR:
            baseline=baseline_dwt(s[1],level=wlevel,wavelet=wname,max_iter=100)
            if self.rtypeM=='none':
                Yf=s[1]-baseline
                y = s[1]
                b = baseline
            if self.rtypeM=='sqrt':
                Yf=(s[1]-baseline)**2
                y = s[1]**2
                b = baseline**2
            if self.rtypeM=='log10':
                Yf=10**(s[1]-baseline)
                y = 10**s[1]
                b = 10**baseline
            x=s[0]
            s.append(baseline)
            if self.Scale.get()==1:
                self.__plt__.set_yscale('log')
            labelB = str(s[2]['tiltAngle'])
            self.__plt__.plot(s[0], Yf, '-', label=labelB, markersize=0.5)
            s.append(Yf)
            self.listSpectraB.append(s)
            addSpectra(self, x, Yf, y, s[2])
            #self.listSpectraB.append(Yf)
        self.__plt__.set_xlabel('Energy  [keV]')
        self.__plt__.set_ylabel('Yield')
        self.legend=self.__plt__.legend()
        fig.canvas.draw()
        self.plot_background_list(self.root)
        
        #self.dataInfo['waveletname']=wname
        #self.dataInfo['waveletlevel']=wlevel      


        #{} self.__plt__.plot(s[0], s[1], '-', label=s[2]['tiltAngle'], markersize=0.5)

    def plot_background_list(self, root):
        new_window_thread = threading.Thread(target=self.plot_background_list2(root))
        new_window_thread.start() 

    def plot_background_list2(self, root):
        chart_window = Toplevel(root)
        chart_window.title("Plot")
        chart_window.geometry("1600x1000")
        fig2 = plt.figure(figsize=(16, 10))

        dl = len(self.listSpectraB)
        ncols_dl = int(ceil(dl/4))
        nrows=4  
        axX = 1
        axY = 1
        axN = 1

        for sp in self.listSpectraB:
            plt.subplot(ncols_dl, nrows, axN) 
            plt.plot(sp[0], sp[1], '.', markersize=1, color='black')
            plt.plot(sp[0], sp[4], '--', color='red')
            plt.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
            plt.grid(which='minor', linestyle=':', linewidth='0.5', color='lightgray')
            plt.minorticks_on()
            plt.yscale('log')
            plt.xlim(2.5, 8)
            plt.xlabel('Energy')
            if axN == 1 or axN == 5 or axN == 9 or axN == 13:
                plt.ylabel('Yield')
            else: 
                plt.gca().set_yticklabels([])
                plt.grid(True, axis='y')
            titleB = str(sp[2]['tiltAngle'])
            plt.title(titleB, loc='center')
            axN = axN + 1
        plt.tight_layout()
        canvas_2 = tk.Canvas()
        canvas_2 = FigureCanvasTkAgg(fig2, master=chart_window) 
        canvas_2.get_tk_widget().place(x=10, y=10)  
        canvas_2.draw()
        toolbar_2 = NavigationToolbar2Tk(canvas_2, chart_window)
        toolbar_2.update()
        toolbar_2.place(x=0, y=0)    

##############################################
    def saveMany(self):
        # wybierz folder
        from datetime import datetime
        idir= filedialog.askdirectory(title="Select folder")
        for sp in self.listSpectraB:
            name = sp[2]['tiltAngle'] + '.txt'
            with open(f"{idir}/{name}", 'w') as fid:
                fid.write('#    DateTime: '+datetime.today().strftime('%Y-%m-%d %H:%M:%S')+'\n')
                fid.write('#    Material: ' + sp[2]['material'] + '\n')
                fid.write('#    Axis: ' + sp[2]['axis']  + '\n')
                fid.write('#    Plane: ' + sp[2]['plane']  + '\n')
                fid.write('#    Tilt angle: ' + sp[2]['tiltAngle']  + '\n')
                fid.write('#    Dose: ' + sp[2]['dose']  + '\n')
                fid.write('#    Energy [keV]    Raw Yield       Yield   \n')
            for i in range(0,len(sp[0][:])):
                ftok=f'{sp[0][i]:9.6f} {sp[1][i]:12.6f} {sp[5][i]:12.6f}\n'
                fid.write(ftok)
            fid.close()
   
    


########################################

    def setXLim(self):
        if self.__plt__ is None:
            return

        xlim=self.__plt__.get_xlim()
        xmin,xmax=xlim[0],xlim[1]

        X=self.spectrum.Xi

        ind=np.argwhere( (X>xmin) & (X<xmax) )
        imin,imax=ind[0,0],ind[-1,0]

        self.spectra_range.set(str(imin)+' '+SEP+' '+str(imax))
        
        
        
 
    def comboChangedMapping(self,event):
        evalue=event.widget.get()        

        state0='disabled' if evalue=='Ignore' else 'enabled'                         
        self.edits['eab']['state']=state0
        self.buttons['mappPlot']['state']=state0
        


    def comboChangedFiterMethods(self,event):
        evalue=event.widget.get()                    
        
        state0='disabled' if evalue=='Ignore' else 'enabled'   
        state1='disabled' if evalue!=self.filtr.methods[2] else 'enabled'
        
        self.buttons['BfilterPot']['state']=state0
        self.edits['ekersize']['state']=state0
        self.comboBoxes['filterPoly']['state']=state1
        self.comboBoxes['filterDeriv']['state']=state1

        
        
    def comboChangedWaveletFamilyM(self,event):
        evalue=event.widget.get().split()        
        wlcollectionM=pywt.wavelist(evalue[-1])
        self.waveletNameM.set(wlcollectionM[3])
        self.comboBoxes['waveletName']['values']=wlcollectionM
            
        maxLevelsM=10 #pywt.dwt_max_level(to_-from_,self.waveletNameM.get()) 
        self.waveletLevelsM.set(maxLevelsM)
        self.comboBoxes['waveletLevels']['values']=' '.join(map(str,np.arange(1,maxLevelsM+1,1)))
             
    def comboChangedWaveletFamily(self,event):
        evalue=event.widget.get().split()        
        wlcollection=pywt.wavelist(evalue[-1])
        self.waveletName.set(wlcollection[3])
        self.comboBoxes['waveletName']['values']=wlcollection
        
        #spectra_range=self.spectra_range.get()
        #fromTo=spectra_range.split(SEP)
        #from_=int(fromTo[0])
        #to_  =int(fromTo[1])
        
        #if to_<0:
         #   spSize=self.spectrum.dataSize()
          #  to_+=spSize
            
        maxLevels=12 #pywt.dwt_max_level(to_-from_,self.waveletName.get()) 
        self.waveletLevels.set(maxLevels)
        self.comboBoxes['waveletLevels']['values']=' '.join(map(str,np.arange(1,maxLevels+1,1)))
    
    def gridOnOff(self):
        if self.__plt__ is None:
            return
        
        self.__plt__.grid()
        self.__fig__.canvas.draw()
        
        
        
    def legendOnOff(self):  
        if self.__plt__ is None:
            return
        
        self.legendSwitch = not self.legendSwitch
        self.legend.set_visible(self.legendSwitch)                
        self.__fig__.canvas.draw()
    
    
    
    def showClipPlot(self):        
        import io
        b=io.BytesIO()
        self.__fig__.savefig(b,format='jpeg')
        
        child=clipPlotWindow(self.root)
        child.drawImage(b)
        
    
        
            
    def cwtPlot(self):        
        if self.dataReduce():                              
            cwt=showCWT(self.root)
            cwt.X,cwt.Y=self.X,self.Y
            cwt.plot()
        
        





    def saveData(self):        
        from datetime import datetime
        idir=iniSet.getValue('default','inipath')
        if idir == '' : idir=self.currDir
        dRes = fd.asksaveasfile(title='Save a file', initialdir=idir,                                                                
                                filetypes=[("Background free TXT","*.txt"),
                                           #("Background free BIN","*.bfbin"),
                                           ("Background free NPY","*.npy"),
                                           ("All files","*.*")])
        if dRes is None:
            print('no file selected')
            return
        
        x=self.spectrum.Xb
        print('x',x)
        y=self.spectrum.Yb
        inf=self.spectrum.Info
        #yfr=self.Yfiltr
        #b=self.baseline
        #ynb=self.Yfree  
        ext=dRes.name[-3:]
        
        #if ext=='txt':
        fid=open(dRes.name,'w')
        
             #fid.write('#ver: 0 \n')
        fid.write('#    DateTime: '+datetime.today().strftime('%Y-%m-%d %H:%M:%S')+'\n')
        fid.write('#    Spectrum file: '+ self.spectrum.ifileName +'\n')
        fid.write('#    Material: ' + self.spectrum.Material + '\n')
        fid.write('#    Axis: ' + self.spectrum.AxialChannel + '\n')
        fid.write('#    Plane: ' + self.spectrum.Plane + '\n')
        fid.write('#    Tilt angle: ' + self.spectrum.TiltAngle + '\n')
        fid.write('#    Dose: ' + self.spectrum.Dose + '\n')

        
        #fid.write('#    Data Range: '+self.dataInfo['dataRange']+'\n')        
        #fid.write('#    Reduce method: '+self.dataInfo['reduceMethod']+'\n')
        #fid.write('#    Calibration method: '+self.dataInfo['cMethod']+'\n')
        #fid.write('#cParams: '+self.dataInfo['axb']+'\n')            
        #fid.write('#    Filter method: '+self.dataInfo['filterMethod']+'\n')
        #fid.write('#    Filter size: '+str(self.dataInfo['filterSize'])+'\n')
        #fid.write('#    Filter poly: '+str(self.dataInfo['filterPoly'])+'\n')
        #fid.write('#    Filter deriv: '+str(self.dataInfo['filterDeriv'])+'\n')
        #fid.write('#    Wavelet name: '+str(self.dataInfo['waveletname'])+'\n')
        #fid.write('#    Wavelet level: '+str(self.dataInfo['waveletlevel'])+'\n')
        #fid.write('#    Data size RC: '+str(self.X.shape[0])+' 4\n')
        #fid.write('#--- x ------ input ------ baseline ------ output ----\n')
        fid.write('#    Energy [keV]    Yield   \n')
        #str(inf[1]), str(inf[2]), str(inf[3])+'\n')                                    
        for i in range(0,len(x[:])):
            ftok=f'{x[i]:9.6f} {y[i]:12.6f}\n'
            fid.write(ftok)                                             
        fid.close()
        #else:
         #   if ext =='bfbin':
         #       return
         #   else:
         #       ext=dRes.name[-3:]
         #       if ext=='npy':                    
         #           np.save(dRes.name,np.array([x,yfr,b,ynb]))
         #       else:
         #           tk.messagebox.showerror('ERROR',message='data not saved: wrong file type')
                    
        
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
       
        


    def __init__(self, root):
        super().__init__()

        self.currDir=os.getcwd()

        self.root = root
        self.root.title("pixek")
        self.root.geometry("1200x700+100+100")
        #self.root.bind('<Destroy>',self.Destroy)
        self.root.protocol("WM_DELETE_WINDOW",self.Destroy)
        #self.root.bind('<Close>',self.Destroy)
        ttk.Style("cyborg") #superhero")
        #ttk.Style("darkly")

        ######################

        labelColor = "danger"
        entryColor = "light"
        buttonColor = "succes"


        ### Loading of initial parameters from file "config.json"

        #initial_parameters = InitialParameters.load_config("config.json")

        ### Buttonbar in the top of app frame

        buttonbar = ttk.Frame(self.root, bootstyle="primary")
        buttonbar.pack(fill=X, pady=5, side=TOP)

        btn1 = ttk.Button(buttonbar, text="Add spectrum", command=self.show_frame_spectra, bootstyle="primary")
        btn1.pack(side=tk.LEFT, padx=5)
        btn2 = ttk.Button(buttonbar, text="Multispectra analysis", command=self.show_frame_plot, bootstyle="primary")
        btn2.pack(side=tk.LEFT, padx=5)
        btn3 = ttk.Button(buttonbar, text="Counts", command=self.show_frame_count, bootstyle="primary")
        btn3.pack(side=tk.LEFT, padx=5)
        btn4 = ttk.Button(buttonbar, text="Angular scan", command=self.show_frame_angular_scan, bootstyle="primary")
        btn4.pack(side=tk.LEFT, padx=5)
        btn5 = ttk.Button(buttonbar, text="Info", command=self.show_info, bootstyle="primary")
        btn5.pack(side=tk.LEFT, padx=5)

        ### Main app frame

        self.notebook = ttk.Frame(self.root)
        self.notebook.pack(fill="both", expand=True) 

        self.notebook_1 = ttk.Frame(self.notebook)
        self.notebook_1.pack(fill="both", side=LEFT) 

        self.notebook_2 = ttk.Frame(self.notebook) 
        self.notebook_3 = ttk.Frame(self.notebook) 
        self.notebook_4 = ttk.Frame(self.notebook)
        self.notebook_5 = ttk.Frame(self.notebook)

        # Info frame! (no place for plot)

        self.frame_6 = ttk.Frame(self.notebook, style='primary.TFrame')
        self.frame_6.pack(fill="both", expand=True, side=RIGHT, pady=10, padx=2)

        ### Plotting frame

        self.__fig__ = plt.figure(figsize=(4,3), dpi = 150)
        

        self.canvas_1 = tk.Canvas()
        self.canvas_1= FigureCanvasTkAgg(self.__fig__, master=self.frame_6)        
        self.canvas_1.get_tk_widget().pack(side=TOP,fill=BOTH,expand=True)

        self.toolbar_1 = NavigationToolbar2Tk(self.canvas_1,self.frame_6,pack_toolbar=False)
        self.toolbar_1.update()
        #self.toolbar_1.place(x=25, y=840)
        self.toolbar_1.pack(side=TOP,fill=BOTH)


        #self.spacePng=tk.PhotoImage(file='space.png')
        #,image=self.spacePng
        
        s=ttk.Style()
        s.configure('Bfromto.TButton',
                    width=24,height=32,
                    borderwidth=0)
        
        self.gridPng=tk.PhotoImage(file='grid.png')        
        BgraphLimits=ttk.Button(self.toolbar_1,image=self.gridPng,command=self.gridOnOff,style='Bfromto.TButton',padding=0)
        BgraphLimits.pack(side=tk.LEFT,padx=0)
        
        self.legendPng=tk.PhotoImage(file='legend.png')        
        BgraphLimits=ttk.Button(self.toolbar_1,image=self.legendPng,command=self.legendOnOff,style='Bfromto.TButton',padding=0)
        BgraphLimits.pack(side=tk.LEFT,padx=0)
        
        BgraphLimits=ttk.Button(self.toolbar_1,text="Accept range",style='Bfromto.TButton',command=self.setXLim)
        BgraphLimits.pack(side=tk.LEFT,padx=5,pady=5)
        self.tooltip=ToolTip(BgraphLimits) #,msg="accept data range",delay=2.0)
        
        BclipPlot=ttk.Button(self.toolbar_1,text="CLONE",style='Bfromto.TButton',command=self.showClipPlot)
        BclipPlot.pack(side=tk.LEFT,padx=5,pady=5)
        
        #self.ScalePlot=ToolTip(BclipPlot)
        self.Scale = tk.IntVar()
        self.Scale.set(1)
        LabelScale = ttk.Label(self.toolbar_1, text="Scale:").pack(side=tk.LEFT,padx=5,pady=5)
        self.Scale0 = ttk.Radiobutton(self.toolbar_1, text='none', variable=self.Scale, value=0).pack(side=tk.LEFT,padx=5,pady=5)
        self.Scale1 = ttk.Radiobutton(self.toolbar_1, text='LOG', variable=self.Scale, value=1).pack(side=tk.LEFT,padx=5,pady=5)

        #self.new_window_plot = 0
        #ttk.Button(self.frame_6, text="New window", style="TButton", command=lambda: new_window(self,self.root)).place(x=28, y=28)

        ### Add spectra options (notebook_1)

        self.add_spectra_steps = ["Load data", "Calibration", "Data reduction", "Remove noise", "Remove background", "Save"]
        
        self.frameLoadData = ttk.Frame(self.notebook_1)
        self.frameCalibration = ttk.Frame(self.notebook_1)
        self.frameDataReduce = ttk.Frame(self.notebook_1)
        self.frameFiltering = ttk.Frame(self.notebook_1)
        self.frameWavelets = ttk.Frame(self.notebook_1)
        self.frame_5 = ttk.Frame(self.notebook_1)

        for step in self.add_spectra_steps:
            frame_step = ttk.Frame(self.notebook_1)
            frame_step.pack(pady=5, fill=X)
            
            if step == "Load data":
                frame = self.frameLoadData
            if step == "Calibration":
                frame = self.frameCalibration
            if step == "Data reduction":
                frame = self.frameDataReduce
            if step == "Remove noise":
                frame = self.frameFiltering
            if step == "Remove background":
                frame = self.frameWavelets
            if step == "Save":
                frame = self.frame_5

            var_step = IntVar() # sprawdza ktory przycisk jest klikniety
            var_step.set(0)

            self.bt_1 = ttk.Checkbutton(frame_step, width=40, text=step, variable = var_step, style='Toolbutton')
            self.bt_1.config(command=lambda var=var_step, btn=self.bt_1, sub=frame, frame=frame_step: self.toggle(var, btn, sub, frame))
            self.bt_1.pack(pady=5, padx=2)

        ####

        self.list_of_spectra = []

        ### Frame 1 - notebook 1 (Add spectra) #######
        ################ Load data ###################

        ttk.Label(self.frameLoadData, text="Material").place(x=15, y=20)
        self.material = StringVar()
        self.material.set(iniSet.getValue('info','material'))
        self.material_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.material, bootstyle=entryColor).place(x=160, y=15)
        
        ttk.Label(self.frameLoadData, text="Axis").place(x=15, y=50)
        self.axis = StringVar()
        self.axis.set(iniSet.getValue('info','axis'))
        self.axis_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.axis, bootstyle=entryColor).place(x=160, y=45)
        
        ttk.Label(self.frameLoadData, text="Plane").place(x=15, y=80)
        self.plane = StringVar()
        self.plane.set(iniSet.getValue('info','plane'))
        self.plane_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.plane, bootstyle=entryColor).place(x=160, y=75)
        
        ttk.Label(self.frameLoadData, text="Tilt angle [deg]").place(x=15, y=110)
        self.angle = StringVar()
        self.angle.set(iniSet.getValue('info','angle'))
        self.angle_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.angle, bootstyle=entryColor).place(x=160, y=105)
        
        ttk.Label(self.frameLoadData, text="Dose [μC]").place(x=15, y=140)
        self.dose_e = StringVar()
        self.dose_e.set(iniSet.getValue('info','dose_e'))
        self.dose_e_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.dose_e, bootstyle=entryColor).place(x=160, y=135)
        
        ttk.Label(self.frameLoadData, text="Plot with dose [μC] as...").place(x=15, y=170)
        self.dose = StringVar()
        self.dose.set(iniSet.getValue('info','dose'))
        self.dose_entry = ttk.Entry(self.frameLoadData, width=8, textvariable=self.dose, bootstyle=entryColor).place(x=160, y=165)
        
        #  Load data from file

        ttk.Button(self.frameLoadData, text="Open file", command=self.openFile).place(x=10, y=210)
        
        self.open_spectra_file = StringVar()
        #self.open_spectra_file_loc = initial_parameters.init_file_loc
        inpFile=iniSet.getValue('default','inpfile')

        if inpFile=='':
            self.open_spectra_file.set('NO FILE SELECTED')            
        else:
            self.open_spectra_file.set(inpFile)            
            self.spectrum.loadSpectrum(iniSet.getValue('default','inipath')+'/'+self.open_spectra_file.get())

        self.open_spectra_file_entry = ttk.Entry(self.frameLoadData, width=22, textvariable=self.open_spectra_file, bootstyle=entryColor).place(x=10, y=250)

        self.spectra_size=StringVar()
        if self.spectrum.empty():
            self.spectra_size.set("Number of points:")
        else:
            self.spectra_size.set("Number of points: "+str(self.spectrum.dataSize()))

        ttk.Label(self.frameLoadData, textvariable=self.spectra_size).place(x=10,y=280)
        
        if inpFile=='':
            srange="NO DATA"
        else:
            inirange= iniSet.getValue('reduction','fromto')
            srange= inirange if inirange!='' else '0 ; -1'



        ### Data Reduction ####

        self.spectra_range=StringVar()
        self.spectra_range.set(srange)
        editArgs={'width': 15, 'justify' : 'center', 'text': self.spectra_range}
        labelArgs={'justify': 'right','width' :15}
        LabeledInput(self.frameDataReduce,label='Range [from;to]',input_class=ttk.Entry,input_args=editArgs,label_args=labelArgs).place(x=15,y=20)

        sreduction='none'

        self.spectra_reduction=tk.StringVar()
        self.spectra_reduction.set(sreduction)
        comboArgs={'width': 13, 'justify': 'center','values': self.spectrum.reducTypes}
        LabeledInput(self.frameDataReduce,label='Reduction ',input_class=ttk.Combobox,input_args=comboArgs,label_args=labelArgs, input_var=self.spectra_reduction).place(x=15,y=60)

        ttk.Button(self.frameDataReduce, text="Reduce", command=self.replotInputData).place(x=15, y=100)

        ### Frame 2 - notebook 1 (Add spectra) #######
        ################ Calibration #################

        #onoff=iniSet.getValue('calibration','onoff')
        cMethod=iniSet.getValue('calibration','method')
        
        self.spectra_mapping=tk.StringVar()
        self.spectra_mapping.set(cMethod)
        comboArgs={'width': 13, 'justify': 'center','values': self.xmapping.methods}
        CmBMapping=LabeledInput(self.frameCalibration,label='Method ',input_class=ttk.Combobox,
                                input_args=comboArgs,label_args=labelArgs, input_var=self.spectra_mapping)
        CmBMapping.place(x=20,y=40)
        CmBMapping.input.bind("<<ComboboxSelected>>",self.comboChangedMapping )        
        
        
        bstate='disabled' if self.spectra_mapping.get()=="Ignore" else 'enabled'
        
        axb=iniSet.getValue('calibration','ax+b')

        self.mapping_axb=tk.StringVar()
        self.mapping_axb.set(axb)
        editArgs={'width': 15, 'justify' : 'center', 'text': self.mapping_axb, 'state': bstate}
        labelArgs={'justify': 'right','width' :15}
        Eab=LabeledInput(self.frameCalibration,label='[ a ; b ]',input_class=ttk.Entry,input_args=editArgs,label_args=labelArgs)
        Eab.place(x=20,y=70)
        self.edits['eab']=Eab.input

        ttk.Button(self.frameCalibration, text="Show spectrum maxima", style="TButton", command=self.findPeaks).place(x=15, y=110)

        ttk.Label(self.frameCalibration, text="Element").place(x=25, y=150)
        self.element = StringVar()
        self.element.set("U")
        self.element_entry = ttk.Entry(self.frameCalibration, width=5, textvariable=self.element, bootstyle=entryColor).place(x=80, y=145)

        ttk.Label(self.frameCalibration, text="Energy").place(x=25, y=180)
        #U_L_alfa = 13.5970 # keV
        #U_M_alfa = 3.1703  # keV

        self.peak1energy = StringVar()
        self.peak1energy.set("3.1703")
        self.peak1e_entry = ttk.Entry(self.frameCalibration, width=8, textvariable=self.peak1energy, bootstyle=entryColor).place(x=15, y=210)

        self.peak2energy = StringVar()
        self.peak2energy.set("13.5970")
        self.peak2e_entry = ttk.Entry(self.frameCalibration, width=8, textvariable=self.peak2energy, bootstyle=entryColor).place(x=15, y=240)

        ttk.Label(self.frameCalibration, text="Channel").place(x=85, y=180)

        self.peak1channel = StringVar()
        self.peak1channel.set("1")
        self.peak1channel_entry = ttk.Entry(self.frameCalibration, width=8, textvariable=self.peak1channel, bootstyle=entryColor).place(x=80, y=210)

        self.peak2channel = StringVar()
        self.peak2channel.set("2")
        self.peak2channel_entry = ttk.Entry(self.frameCalibration, width=8, textvariable=self.peak2channel, bootstyle=entryColor).place(x=80, y=240)
        
        

        ttk.Button(self.frameCalibration, text="Mark", style="TButton", command=self.peak_1, width=5).place(x=160, y=210)
        ttk.Button(self.frameCalibration, text="Mark", style="TButton", command=self.peak_2, width=5).place(x=160, y=240)
        
        # Save Button to the list - later callibration
        #ttk.Button(self.frameCalibration, text="Save to list", command=self.saveToList).place(x=15, y=280)
       
        BMappingPlot=ttk.Button(self.frameCalibration, text="PLOT", command=self.replotMappedData,state=bstate)
        BMappingPlot.place(x=160, y=280)
        self.buttons['mappPlot']=BMappingPlot
        
        #ttk.Button(self.frameWavelets, text="Savitzky-Golay filter", style="TButton", command=lambda: SG_filter(self)).place(x=20, y=135)



        ### Frame 3 - notebook 1 (Add spectra) #######

        ############ Remove noise ###############
        
        cMethod=iniSet.getValue('filtering','method')
        
        self.filterMethod=tk.StringVar()
        self.filterMethod.set(cMethod)
        comboArgs={'width': 13, 'justify': 'center','values': self.filtr.methods}
        lecmFl=LabeledInput(self.frameFiltering,label='Method ',input_class=ttk.Combobox,input_args=comboArgs,
                     label_args=labelArgs, input_var=self.filterMethod)
        lecmFl.place(x=20,y=40)                
        lecmFl.input.bind("<<ComboboxSelected>>",self.comboChangedFiterMethods)
        
        bstate='disabled' if self.filterMethod.get()=="Ignore" else 'enabled'
        
        ksize=int(iniSet.getValue('filtering','ksize'))
        self.filterKernelSize=tk.IntVar()
        self.filterKernelSize.set(ksize)
        editArgs={'width': 15, 'justify' : 'center', 'text': self.filterKernelSize, 'state':bstate}
        labelArgs={'justify': 'right','width' :15}
        Ekersize=LabeledInput(self.frameFiltering,label='Kernel size',input_class=ttk.Entry,input_args=editArgs,label_args=labelArgs)
        Ekersize.place(x=20,y=70)
        self.edits['ekersize']=Ekersize.input
        
        BfilterPlot=ttk.Button(self.frameFiltering, text="PLOT", state=bstate,command=self.replotFilteredData)
        BfilterPlot.place(x=160, y=170)        
        self.buttons['BfilterPot']=BfilterPlot        
        
        bstate='enabled' if self.filterMethod.get()==self.filtr.methods[2] else 'disabled'
        
        polyorder=int(iniSet.getValue('filtering','polyorder'))
        self.filterPolyOrder=tk.IntVar()
        self.filterPolyOrder.set(polyorder)
        comboArgs={'width': 13, 'justify': 'center','values': [2,3,4,5], 'state':bstate}
        CmBfilterPoly=LabeledInput(self.frameFiltering,label='Polyroder',input_class=ttk.Combobox,input_args=comboArgs,
                     label_args=labelArgs,input_var=self.filterPolyOrder)
        CmBfilterPoly.place(x=20,y=100)
        self.comboBoxes['filterPoly']=CmBfilterPoly.input
        
        
        deriv=int(iniSet.getValue('filtering','derivative'))
        self.filterDeriv=tk.IntVar()
        self.filterDeriv.set(deriv)
        comboArgs={'width': 13, 'justify': 'center','values': [0,1,2],'state':bstate}
        CmBfilterDeriv=LabeledInput(self.frameFiltering,label='Derivative',input_class=ttk.Combobox,input_args=comboArgs,
                     label_args=labelArgs,input_var=self.filterDeriv)
        CmBfilterDeriv.place(x=20,y=130)
        self.comboBoxes['filterDeriv']=CmBfilterDeriv.input


        plotRawData=int(iniSet.getValue('filtering','plotRaw'))
        self.filterShowRawData=tk.IntVar()
        self.filterShowRawData.set(plotRawData)
        checkArgs={'width': 13, 'justify': 'center'}
        LabeledInput(self.frameFiltering,label="plot raw data",input_class=ttk.Checkbutton,input_var=self.filterShowRawData).place(x=20,y=170)
        
        

        ### Frame 4 - notebook 1 (Add spectra) #######
        ############ Remove background ###############
        
        wfamily=iniSet.getValue('wavelet','family')

        self.waveletFamily=tk.StringVar()        
        self.waveletFamily.set(wfamily)
        comboArgs={'width': 15,
               'justify': 'center',
               'values': ['Daubechies db',
                          'Symlets sym',
                          'Coiflets coif',
                          'Biorthogonal bior',
                          'Reverse biorthogonal rbio',
                          ]  }
        labelArgs={'justify': 'right','width' :15}
        lecmFm=LabeledInput(self.frameWavelets,label='Family ',
                          input_class=ttk.Combobox,input_args=comboArgs, 
                          input_var=self.waveletFamily,label_args=labelArgs)
        lecmFm.place(x=20,y=40)
        lecmFm.input.bind("<<ComboboxSelected>>",self.comboChangedWaveletFamily)
                
        wlcollection=pywt.wavelist(self.waveletFamily.get().split(' ')[-1])
        
        wname=iniSet.getValue('wavelet','name')
        self.waveletName=tk.StringVar()
        self.waveletName.set(wname)
        comboArgs={'width': 15,
               'justify': 'center',
               'values': wlcollection }
        
        lecmFm=LabeledInput(self.frameWavelets,label='name',
                          input_class=ttk.Combobox,input_args=comboArgs,
                          input_var=self.waveletName,label_args=labelArgs)
        lecmFm.place(x=20,y=70)
        self.comboBoxes['waveletName']=lecmFm.input
                        
        
        maxLevels=1 if self.spectrum.empty() else int(iniSet.getValue('wavelet','levels'))
                    
        self.waveletLevels=tk.IntVar()
        self.waveletLevels.set(maxLevels)        
        comboArgs={'width': 15,
               'justify': 'center',
               'values':  ' '.join(map(str,np.arange(1,maxLevels+1,1))) }
        
        CmBoxLevels=LabeledInput(self.frameWavelets,label='Levels',input_class=ttk.Combobox,
                     input_args=comboArgs,input_var=self.waveletLevels,label_args=labelArgs)
        CmBoxLevels.place(x=20,y=100)
        self.comboBoxes['waveletLevels']=CmBoxLevels.input
        
        
        ttk.Button(self.frameWavelets, text="PLOT", command=self.plotRemoveBackground).place(x=160, y=130)        
        
        ### Frame 5 - notebook 1 (Add spectra) #######
        ################ Save ########################

        ttk.Button(self.frame_5, text="Save spectrum", style="TButton", command=self.saveData, width=13).place(x=20, y=40)
        #ttk.Button(self.frame_5, text="Save list of spectra", style="TButton", command=lambda: add_to_list(self, self.spectra_name.get()), width=13).place(x=20, y=80)
        #ttk.Button(self.frame_5, text="Save config parameters", style="TButton", command=lambda: InitialParameters(self.X, self.Y,).save_config("config.json")).place(x=20, y=230)

        ### Notebook 2 ###############################
        ################ Multispectra analysis ########################

        self.frame2_1 = ttk.Frame(self.notebook_2) 
        self.frame2_1.pack(pady=5, padx=5, fill="both")

        self.multispectraOptions = ["Open spectra", "Data Reduction", "Remove background", "Save"]
        
        self.frameOpenSpectra = ttk.Frame(self.notebook_2)
        #self.frameCalibrateMany = ttk.Frame(self.notebook_2)
        self.frameDataReduceMany = ttk.Frame(self.notebook_2)
        #self.frameFilteringMany = ttk.Frame(self.notebook_2)
        self.frameWaveletsMany = ttk.Frame(self.notebook_2)
        self.frameSaveMany = ttk.Frame(self.notebook_2)

        for step in self.multispectraOptions:
            frameStepMany = ttk.Frame(self.notebook_2)
            frameStepMany.pack(pady=5, fill=X)
            
            if step == "Open spectra":
                frame = self.frameOpenSpectra
            if step == "Data Reduction":
                frame = self.frameDataReduceMany
            if step == "Remove background":
                frame = self.frameWaveletsMany
            if step == "Save":
                frame = self.frameSaveMany

            varStepMany = IntVar() # sprawdza ktory przycisk jest klikniety
            varStepMany.set(0)

            self.bt_2 = ttk.Checkbutton(frameStepMany, width=40, text=step, variable = varStepMany, style='Toolbutton')
            self.bt_2.config(command=lambda var=varStepMany, btn=self.bt_2, sub=frame, frame=frameStepMany: self.toggleM(var, btn, sub, frame))
            self.bt_2.pack(pady=5, padx=2)

        ####
        # Open folder
        #l = tk.LabelFrame(self.frameOpenSpectra, text="Input data processing", padx=30, pady=20)
        #l.place(x=2, y=5, width=243, height=400)

        ttk.Button(self.frameOpenSpectra, text="Open", style="TButton", command=lambda: openFolder(self), width=18).pack(pady=30, padx=2)

        #self.list_for_plot = tk.Listbox(self.frameOpenSpectra, selectmode=tk.MULTIPLE)
        #self.list_for_plot.pack(pady=15, padx = 5)
        
        self.tree = ttk.Treeview(self.frameOpenSpectra, columns=("ID", "TiltAngle", "Plane", "Axis", "Material"), show="headings")

        #self.tree = ttk.Treeview(self.frameOpenSpectra, columns=("TiltAngle", "Plane", "Axis", "Material"), show="headings")

        
        self.tree.column("ID",  width=2) #, stretch=False)
        self.tree.column("Material", width=1)
        self.tree.column("Axis", width=1)
        self.tree.column("Plane", width=1)
        self.tree.column("TiltAngle", width=1)
        self.tree.heading("ID", text="ID")
        self.tree.heading("Material", text="Material")
        self.tree.heading("Axis", text="Axis")
        self.tree.heading("Plane", text="Plane")
        self.tree.heading("TiltAngle", text="Tilt angle")
        self.tree.configure(selectmode="extended")
        self.tree.pack(expand=True, fill="both")
        #tree = ttk.Treeview(root, columns=("ID", "Nazwa", "Wartość"), show="headings", selectmode="extended")
        self.tree.bind("<ButtonRelease-1>", self.clickTree)
        #ttk.Button(self.frameOpenSpectra, text="Plot", style="TButton", command=self.plotList, width=13).pack(pady=15, padx=2)

        # Data Reduction
        # srange, sreduction
        #srangeM=srange 
        self.spectra_range_M=StringVar()
        #self.spectra_range_M.set(srangeM)
        editArgs_M={'width': 15, 'justify' : 'center', 'text': self.spectra_range_M}
        labelArgs_M={'justify': 'right','width' :15}
        LabeledInput(self.frameDataReduceMany,label='Range [from;to]',input_class=ttk.Entry,input_args=editArgs_M,label_args=labelArgs_M).place(x=15,y=20)

        sreductionM='none'

        self.spectra_reduction_M=tk.StringVar()
        self.spectra_reduction_M.set(sreductionM)
        comboArgs_M={'width': 13, 'justify': 'center','values': self.spectrum.reducTypes}
        LabeledInput(self.frameDataReduceMany,label='Reduction ',input_class=ttk.Combobox,input_args=comboArgs_M,label_args=labelArgs_M, input_var=self.spectra_reduction_M).place(x=15,y=60)

        ttk.Button(self.frameDataReduceMany, text="Reduce", command=self.replotM).place(x=15, y=100)

        ### Filtering ####
        
        # Background
        
        wfamilyM=iniSet.getValue('wavelet','family')

        self.waveletFamilyM=tk.StringVar()        
        self.waveletFamilyM.set(wfamilyM)
        comboArgsM={'width': 15,
               'justify': 'center',
               'values': ['Daubechies db',
                          'Symlets sym',
                          'Coiflets coif',
                          'Biorthogonal bior',
                          'Reverse biorthogonal rbio',
                          ]  }
        #labelArgs={'justify': 'right','width' :15}
        lecmFmM=LabeledInput(self.frameWaveletsMany,label='Family ',
                          input_class=ttk.Combobox,input_args=comboArgsM, 
                          input_var=self.waveletFamilyM,label_args=labelArgs)
        lecmFmM.place(x=20,y=40)
        lecmFmM.input.bind("<<ComboboxSelected>>",self.comboChangedWaveletFamilyM)
                
        wlcollectionM=pywt.wavelist(self.waveletFamilyM.get().split(' ')[-1])
        
        wnameM=iniSet.getValue('wavelet','name')
        self.waveletNameM=tk.StringVar()
        self.waveletNameM.set(wnameM)
        comboArgs={'width': 15,
               'justify': 'center',
               'values': wlcollectionM }
        lecmFmM=LabeledInput(self.frameWaveletsMany,label='Name',
                          input_class=ttk.Combobox,input_args=comboArgs,
                          input_var=self.waveletNameM,label_args=labelArgs)
        lecmFmM.place(x=20,y=70)
        self.comboBoxes['waveletName']=lecmFmM.input
                        
        maxLevelsM=1 if self.spectrum.empty() else int(iniSet.getValue('wavelet','levels'))
                    
        self.waveletLevelsM=tk.IntVar()
        self.waveletLevelsM.set(maxLevelsM)        
        comboArgs={'width': 15,
               'justify': 'center',
               'values':  ' '.join(map(str,np.arange(1,maxLevelsM+1,1))) }
        CmBoxLevelsM=LabeledInput(self.frameWaveletsMany,label='Levels',input_class=ttk.Combobox,
                     input_args=comboArgs,input_var=self.waveletLevelsM,label_args=labelArgs)
        CmBoxLevelsM.place(x=20,y=100)
        self.comboBoxes['waveletLevels']=CmBoxLevelsM.input
        
        ttk.Button(self.frameWaveletsMany, text="PLOT", command=self.plotRemoveBackgroundM).place(x=160, y=130)      

        # Save

        ttk.Button(self.frameSaveMany, text="Save", command=self.saveMany).place(x=30, y=30)

        ### Notebook 3 #########################
        ######### Counts #######################

        self.frame3_1 = ttk.Frame(self.notebook_3) 
        self.frame3_1.pack(pady=5, padx=5, fill="both")

        self.CountsOptions = ["PIXE/C", "RBS/C"]
        
        self.framePIXEC = ttk.Frame(self.notebook_3)
        self.frameRBSC = ttk.Frame(self.notebook_3)


        for step in self.CountsOptions:
            frameStepC = ttk.Frame(self.notebook_3)
            frameStepC.pack(pady=5, fill=X)
            
            if step == "PIXE/C":
                frame = self.framePIXEC
            if step == "RBS/C":
                frame = self.frameRBSC

            varStepC = IntVar() # sprawdza ktory przycisk jest klikniety
            varStepC.set(0)

            self.bt_3 = ttk.Checkbutton(frameStepC, width=40, text=step, variable = varStepC, style='Toolbutton')
            self.bt_3.config(command=lambda var=varStepC, btn=self.bt_3, sub=frame, frame=frameStepC: self.toggleC(var, btn, sub, frame))
            self.bt_3.pack(pady=5, padx=2)
        ### Show list and select spectra for counts

        ttk.Button(self.framePIXEC, text="Add spectrum", style="TButton", command=lambda: loadData(self), width=25).place(y=200, x=2)
        
        self.treeC = ttk.Treeview(self.framePIXEC, columns=("ID", "TiltAngle", "Plane", "Axis", "Material"), show="headings", height=10)
        self.treeC.heading("ID", text="ID")
        self.treeC.heading("Material", text="Material")
        self.treeC.heading("Axis", text="Axis")
        self.treeC.heading("Plane", text="Plane")
        self.treeC.heading("TiltAngle", text="Ang.")
        self.treeC.column("ID", width=1, stretch=True)
        self.treeC.column("Material", width=1)
        self.treeC.column("Axis", width=1)
        self.treeC.column("Plane", width=1)
        self.treeC.column("TiltAngle", width=1)
        #self.treeC.configure(selectmode="extended")
        self.treeC.pack(side="top", pady=5, fill="x")
        
        self.treeC.bind("<ButtonRelease-1>", self.clickTreeC)
        
        ### Select energy range of peaks
        
        ttk.Label(self.framePIXEC, text="Define energy range of selected peak",  width=40).place(x=20, y=250)
        
        ttk.Label(self.framePIXEC, text="E_min [keV]").place(x=20, y=290)
        self.emin_as = StringVar()
        self.emin_as_entry = ttk.Entry(self.framePIXEC, width=5, textvariable=self.emin_as)
        self.emin_as_entry.place(x=100, y=290)
        self.emin_as.set('4.4')
        ttk.Label(self.framePIXEC, text="E_max [keV]").place(x=20, y=330)
        self.emax_as = StringVar()
        self.emax_as_entry = ttk.Entry(self.framePIXEC, width=5, textvariable=self.emax_as)
        self.emax_as_entry.place(x=100, y=330)
        self.emax_as.set('5.6')

        ttk.Button(self.framePIXEC, text="Show on plot", style="TButton", command=lambda: showRange(self), width=13).place(x=37, y=370)
        ttk.Button(self.framePIXEC, text="Calculate", style="TButton", command=lambda: showPlotCounts(self), width=13).place(x=37, y=410)
        ttk.Button(self.framePIXEC, text="Calculate - R = 1", style="TButton", command=lambda: showPlotCountsN(self), width=13).place(x=37, y=450)
        #ttk.Button(self.framePIXEC, text="Show plot", style="TButton", command=lambda: show_plot_counts(self), width=13).place(x=37, y=450)
        ttk.Button(self.framePIXEC, text="Save to file", style="TButton", command=lambda: saveCount(self), width=13).place(x=37, y=490)
        
        ## RBS/C
        # Dane: wybierz plik, zdefiniuj kat, dawka
        # First add random
        # Then attach files with scan from differeat angles
        # Specify range of channels for counts
        # Give name to your data
        # Add list of counts to plot
        ttk.Label(self.frameRBSC, text="Load random spectrum").place(x=20, y=20)
        ttk.Button(self.frameRBSC, text="Open file", style="TButton", command=lambda: openRandom(self), width=25).place(x=20, y=50)
        ttk.Label(self.frameRBSC, text="Dose").place(x=20, y=90)
        self.doseRBS = StringVar()
        self.doseRBS_entry = ttk.Entry(self.frameRBSC, width=5, textvariable=self.doseRBS)
        self.doseRBS_entry.place(x=100, y=90)
        self.doseRBS.set('10')

        ttk.Label(self.frameRBSC, text="Load spectra").place(x=20, y=120)
        ttk.Button(self.frameRBSC, text="Open file", style="TButton", command=lambda: openRBSC(self), width=25).place(x=20, y=150)
        ttk.Label(self.frameRBSC, text="Angle").place(x=20, y=190)
        self.angleRBS = StringVar()
        self.angleRBS_entry = ttk.Entry(self.frameRBSC, width=5, textvariable=self.angleRBS)
        self.angleRBS_entry.place(x=100, y=190)
        self.angleRBS.set('0')
        ttk.Label(self.frameRBSC, text="Dose").place(x=20, y=220)
        self.doseRBSa = StringVar()
        self.doseRBSa_entry = ttk.Entry(self.frameRBSC, width=5, textvariable=self.doseRBSa)
        self.doseRBSa_entry.place(x=100, y=220)
        self.doseRBSa.set('10')

        ttk.Label(self.frameRBSC, text="Channel range").place(x=20, y=250)
        self.channelRBSmin = StringVar()
        self.channelRBSmin_entry = ttk.Entry(self.frameRBSC, width=4, textvariable=self.channelRBSmin)
        self.channelRBSmin_entry.place(x=50, y=280)
        self.channelRBSmin.set('1')
        self.channelRBSmax = StringVar()
        self.channelRBSmax_entry = ttk.Entry(self.frameRBSC, width=4, textvariable=self.channelRBSmax)
        self.channelRBSmax_entry.place(x=100, y=280)
        self.channelRBSmax.set('2000')
        ttk.Button(self.frameRBSC, text="Show", style="TButton", command=lambda: showRBSlim(self), width=10).place(x=180, y=280)
        ttk.Button(self.frameRBSC, text="Calculate", style="TButton", command=lambda: RBSas(self), width=10).place(x=20, y=320)
        
        ttk.Label(self.frameRBSC, text="Name").place(x=20, y=360)
        self.nameRBS = StringVar()
        self.nameRBS_entry = ttk.Entry(self.frameRBSC, width=10, textvariable=self.nameRBS)
        self.nameRBS_entry.place(x=100, y=360)
        ttk.Button(self.frameRBSC, text="Save", style="TButton", command=lambda: RBSsave(self), width=10).place(x=20, y=400)
        






        #########################
        ### Angular Scanning ####
        
        self.frame4_1 = ttk.Frame(self.notebook_4) 
        self.frame4_1.pack(pady=5, padx=5, fill="both")

        self.ASOptions = ["Load scans", "Plot options", "Save"]
        
        self.frameOpenAS = ttk.Frame(self.notebook_4)
        self.framePlotAS = ttk.Frame(self.notebook_4)
        self.frameSaveAS = ttk.Frame(self.notebook_4)


        for step in self.ASOptions:
            frameStepAS = ttk.Frame(self.notebook_4)
            frameStepAS.pack(pady=5, fill=X)
            
            if step == "Load scans":
                frame = self.frameOpenAS
            if step == "Plot options":
                frame = self.framePlotAS
            if step == "Save":
                frame = self.frameSaveAS

            varStepAS = IntVar() # sprawdza ktory przycisk jest klikniety
            varStepAS.set(0)

            self.bt_4 = ttk.Checkbutton(frameStepAS, width=40, text=step, variable = varStepAS, style='Toolbutton')
            self.bt_4.config(command=lambda var=varStepAS, btn=self.bt_4, sub=frame, frame=frameStepAS: self.toggleAS(var, btn, sub, frame))
            self.bt_4.pack(pady=5, padx=2)
        
        
        self.treeAS = ttk.Treeview(self.frameOpenAS, columns=("ID", "Name"), height=10, show="headings")
        self.treeAS.heading("ID", text="ID")
        self.treeAS.heading("Name", text="Name")
        self.treeAS.column("ID", width=1, stretch=True)
        self.treeAS.column("Name", width=1, stretch=True)
        self.treeAS.pack(side="top", pady=5, fill="x")
        self.treeAS.bind("<ButtonRelease-1>", self.clickTreeAS)




        ### Notebook 4 ###############################
        ################ INFO ########################

### DEF ######################################################

    def show_gui(*args):
        return AppGUI(self.root)

    def spectra_to_dict(spectra):
        return spectra.__dict__

    def saveToList(self):
        spectrum_ = [self.X, self.Y, self.spectrum.Info]
        self.list_of_spectra.append(spectrum_)
        messagebox.showinfo('Info', 'Done!')
        print(self.spectrum.Info)

    def saveList(self):
        file_name = self.file_list.get()
        with open(file_name, 'w', newline='') as csvfile:
        #fieldnames = ['Energy', 'Yield', 'Info']
            writer = csv.DictWriter(csvfile, fieldnames=['Energy [keV]', 'Yield', 'Info'])
            writer.writeheader()
            for spectrum_ in self.list_of_spectra:
                writer.writerow(self.X, self.Y, self.spectrum.Info)
        messagebox.showinfo('Info', 'Done!')


    def clickTree(self, event):
        self.selectedSpectrum = self.tree.selection()
        if self.selectedSpectrum in self.selectedSpectra:
            self.selectedSpectra.remove(self.selectedSpectrum)
        else:
            self.selectedSpectra.append(self.selectedSpectrum)
        self.showSpectra()

    def clickTreeC(self, event):
        self.selectedSpectrum = self.treeC.selection()
        if self.selectedSpectrum in self.selectedSpectraC:
            self.selectedSpectraC.remove(self.selectedSpectrum)
        else:
            self.selectedSpectraC.append(self.selectedSpectrum)
        print('d')
        plotSpectraCount(self)
        
    def clickTreeAS(self, event):
        self.selectedScan = self.treeAS.selection()
        if self.selectedScan in self.selectedScans:
            self.selectedScans.remove(self.selectedScan)
        else:
            self.selectedScans.append(self.selectedScan)
        showplotAS(self)

    def peak_1(self):
        self.cid = self.canvas_1.mpl_connect('button_press_event', self.onclick_1)
    def onclick_1(self, event):
        self.x_x = event.xdata
        self.y_y = event.ydata
        Y_u = []
        X_u = []
        Y_u.append(0)
        Y_u.append(self.spectrum.Yi.max())
        X_u.append(self.x_x)
        X_u.append(self.x_x)
        plt.plot(X_u, Y_u, linestyle='-', color='red')
        self.canvas_1.draw()
        self.peak1channel.set(self.x_x)
        self.canvas_1.mpl_disconnect(self.cid)
    
    def peak_2(self):
        self.cid = self.canvas_1.mpl_connect('button_press_event', self.onclick_2)
    def onclick_2(self, event):
        self.x_x = event.xdata
        self.y_y = event.ydata
        Y_u = []
        X_u = []
        Y_u.append(0)
        Y_u.append(self.spectrum.Yi.max())
        X_u.append(self.x_x)
        X_u.append(self.x_x)
        plt.plot(X_u, Y_u, linestyle='-', color='red')
        self.canvas_1.draw()
        self.peak2channel.set(self.x_x)
        self.canvas_1.mpl_disconnect(self.cid)

    def findPeaks(self):
        npoints=self.spectrum.Xi.shape[0]
        peaks, _ = find_peaks(self.spectrum.Yi, height=10, threshold=10, distance=50)
        P = []
        for peak in peaks:
            peak = peak + 0
            P.append(peak)
        P_array = np.array(P)
        plt.plot(P_array, self.spectrum.Yi[peaks], "*", color='red')
        self.canvas_1.draw()  

    def toggleM(self, show, toggle_button, frame_s, frame_n):
        if bool(show.get()):
            frame_s.pack(after=frame_n, fill="both", expand=1, padx=10, pady=5)  # Dodanie sub_frame pod frame
            if frame_s != self.frameOpenSpectra:
                self.frameOpenSpectra.pack_forget()
            if frame_s != self.frameDataReduceMany:
                self.frameDataReduceMany.pack_forget()
            if frame_s != self.frameWaveletsMany:
                self.frameWaveletsMany.pack_forget()
            if frame_s != self.frameSaveMany:
                self.frameSaveMany.pack_forget()
        else:
            frame_s.pack_forget()

    def toggleC(self, show, toggle_button, frame_s, frame_n):
        if bool(show.get()):
            frame_s.pack(after=frame_n, fill="both", expand=1, padx=10, pady=5)  # Dodanie sub_frame pod frame
            if frame_s != self.framePIXEC:
                self.framePIXEC.pack_forget()
            if frame_s != self.frameRBSC:
                self.frameRBSC.pack_forget()
        else:
            frame_s.pack_forget()

    def toggleAS(self, show, toggle_button, frame_s, frame_n):
        if bool(show.get()):
            frame_s.pack(after=frame_n, fill="both", expand=1, padx=10, pady=5)  # Dodanie sub_frame pod frame
            if frame_s != self.frameOpenAS:
                self.frameOpenAS.pack_forget()
            if frame_s != self.framePlotAS:
                self.framePlotAS.pack_forget()
            if frame_s != self.frameSaveAS:
                self.frameSaveAS.pack_forget()
        else:
            frame_s.pack_forget()

    def toggle(self, show, toggle_button, frame_s, frame_n):
        if bool(show.get()):
            frame_s.pack(after=frame_n, fill="both", expand=1, padx=10, pady=5)  # Dodanie sub_frame pod frame
            if frame_s != self.frameLoadData:
                self.frameLoadData.pack_forget()
            if frame_s != self.frameDataReduce:
                self.frameDataReduce.pack_forget()
            if frame_s != self.frameCalibration:
                self.frameCalibration.pack_forget()
            if frame_s != self.frameFiltering:
                self.frameFiltering.pack_forget()
            if frame_s != self.frameWavelets:
                self.frameWavelets.pack_forget()
            if frame_s != self.frame_5:
                self.frame_5.pack_forget()
        else:
            frame_s.pack_forget()

    
    def show_frame_spectra(self):
        self.clear_main_frame()
        self.notebook_1.pack(fill="both", side=LEFT) # expand=True)
        self.frame_6.pack(fill="both", expand=True, side=RIGHT, pady=5, padx=2)
   
    def show_frame_plot(self):
        self.clear_main_frame()
        self.notebook_2.pack(fill="both", side=LEFT) # expand=True)
        self.frame_6.pack(fill="both", expand=True, side=RIGHT, pady=5, padx=2)
    
    def show_frame_count(self):
        self.clear_main_frame()
        self.notebook_3.pack(fill="both", side=LEFT) # expand=True)
        self.frame_6.pack(fill="both", expand=True, side=RIGHT, pady=5, padx=2)

    def show_frame_angular_scan(self):
        self.clear_main_frame()
        self.notebook_4.pack(fill="both", side=LEFT) # expand=True)
        self.frame_6.pack(fill="both", expand=True, side=RIGHT, pady=5, padx=2)

    def show_info(self):
        self.clear_main_frame()
        self.notebook_5.pack(fill="both", side=LEFT) # expand=True)

    def clear_main_frame(self):
        for widget in self.notebook.winfo_children():
            widget.pack_forget()
            
            
    def Destroy(self,*arg):        
        #close = messagebox.askyesno("Exit?", "Are you sure you want to exit?")
        close = True
        if close:            
            iniSet.setValue('default','inpfile',self.open_spectra_file.get())
            
            iniSet.setValue('reduction','name',self.spectra_reduction.get())
            iniSet.setValue('reduction','fromto',self.spectra_range.get())
            
            iniSet.setValue('calibration','method',self.spectra_mapping.get())
            iniSet.setValue('calibration','ax+b',str(self.mapping_axb.get()))
            
            iniSet.setValue('filtering','method',self.filterMethod.get())
            iniSet.setValue('filtering','ksize',str(self.filterKernelSize.get()))
            iniSet.setValue('filtering','polyorder',str(self.filterPolyOrder.get()))
            iniSet.setValue('filtering','derivative',str(self.filterDeriv.get()))
            iniSet.setValue('filtering','plotRaw',str(self.filterShowRawData.get()))
            
            iniSet.setValue('wavelet','family',self.waveletFamily.get())
            iniSet.setValue('wavelet','name',self.waveletName.get())
            iniSet.setValue('wavelet','levels',str(self.waveletLevels.get()))
            
            self.root.destroy()
            
            
    def __del__(self):
        print('close')   
        
        #iniSet.setValue('calibration','ax+b',str(self.spectra_mapping.get()))
        #iniSet.setValue('calibration','ax+b','123;345')

            
            


#s=ttk.Style()
#print(s.theme_names())
iniSet.setDefaultIni()
#matplotlib.use('tkagg')

if __name__ == "__main__":
    root = tk.Tk()
    app = AppGUI(root)
    root.mainloop()
