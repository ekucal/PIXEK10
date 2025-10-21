import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import key_press_handler


from scipy.signal import find_peaks
from skued import gaussian, spectrum_colors
from scipy.signal import savgol_filter
from skued import baseline_dwt

from scipy import stats
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

import numpy as np
from math import sin, pi, exp, sqrt, log

import pywt
from modules.plot import plot_background, plot_background_list 
from modules.plot_list import show_plot_list_c, show_plot_list_c_b
from modules.data import Spectra, list_of_spectra, list_of_spectra_c

def SG_filter(self):
    dl = len(self.spectrum.channel[:])
    energy_array = np.array(self.spectrum.energy)
    counts_array = np.array(self.spectrum.counts)
    sg_min = int(self.sg_min.get())
    sg_max = int(self.sg_max.get())
    energy_array, counts_array = energy_array[sg_min:sg_max], counts_array[sg_min:sg_max]
    npoints=energy_array.shape[0]
    counts_filter = []
    signal = savgol_filter(counts_array, int(np.floor(npoints/12)),8)
    for c in range(0, sg_min):
        counts_filter.append(0)
    for c in range(sg_min, sg_max):
        s = c - sg_min
        counts_filter.append(signal[s])
    for c in range(sg_max, dl):
        counts_filter.append(0)
    self.spectrum.filter = counts_filter
    plt.plot(self.spectrum.energy, self.spectrum.filter, label='Signal with SG filter')
    plt.legend()
    self.canvas_1.draw()  

def my_function(x, a1, a2, a3, a4, a5):
    return np.exp(a1+(a2*x)+(a3*pow(x, 2))+(a4*pow(x, 3))+(a5*pow(x, 4)))

def background(self): #, sg, rb, wavelet_family_names, wava_names, level_number, maxiter, oM, a0, a1, a2, a3, a4, x_1, x_2, y_1, y_2):
    sg = int(self.SGf.get())
    if sg == 1: #yes
        spectrum_counts = self.spectrum.filter
    if sg == 0: #no
        spectrum_counts = self.spectrum.counts
    lev_number = int(self.level_number.get())
    miter = int(self.maxiter.get())
    wava_name = self.wavelet_names.get()
    dl=len(self.spectrum.energy[:])
    spectrum_counts_log = []
    for i in range(0, dl):
        if float(spectrum_counts[i]) > 0:
            spectrum_counts_log.append(log(float(spectrum_counts[i])))
        else:
            spectrum_counts_log.append(0)
    energy_arr = np.array(self.spectrum.energy)
    counts_arr = np.array(spectrum_counts_log)
    background = []
    counts_b = []
    rb = self.rb_entry.get()
    if rb == 'None':
        for i in range(0, dl):
            background.append(0)
            c_b = spectrum_counts[i]
            counts_b.append(c_b)
        self.spectrum.background = background
        self.spectrum.final_yield = counts_b
        label_rb = rb
    if rb == 'Wavelet':
        baseline = baseline_dwt(counts_arr, level = lev_number, max_iter = miter, wavelet = wava_name)
        for i in range(0, dl):
            b = exp(baseline[i])
            background.append(b)
            c_b = spectrum_counts[i] - b
            counts_b.append(c_b)
        self.spectrum.background = background
        self.spectrum.final_yield = counts_b
        label_rb = wava_name
    if rb == 'New':
        baseline = baseline_dwt(counts_arr, level = lev_number, max_iter = miter, wavelet = wava_name)
        for i in range(0, dl):
            b = exp(baseline[i])
            background.append(b)
            c_b = spectrum_counts[i] - b
            counts_b.append(c_b)
        self.spectrum.background = background
        self.spectrum.final_yield = counts_b
        label_rb = rb
    if rb == 'Polynomial':
        oM_in = int(self.oM.get())
        a0_in = float(self.a0.get())
        a1_in = float(self.a1.get())
        a2_in = float(self.a2.get())
        a3_in = float(self.a3.get())
        a4_in = float(self.a4.get())
        for i in range(0, dl):
            Rdata0 = np.array(spectrum_counts)
            maxima0 = argrelextrema(Rdata0, np.less, order=oM_in)
            lfa0 = list(maxima0[0])
            params, covariance = curve_fit(my_function, lfa0, Rdata0[lfa0], p0=[a0_in, a1_in, a2_in, a3_in, a4_in]) 
            a1_fit, a2_fit, a3_fit, a4_fit, a5_fit = params
            fb = my_function(i, a1_fit, a2_fit, a3_fit, a4_fit, a5_fit)
            if fb < 0:
                fb = 0
            background.append(fb)
            c_b = spectrum_counts[i] - fb
            counts_b.append(c_b)
        self.spectrum.background = background
        self.spectrum.final_yield = counts_b
        label_rb = rb
    if rb == 'Linear':
        y1_l = float(self.y_1.get()) 
        y2_l = float(self.y_2.get())
        x1_l = float(self.x_1.get())
        x2_l = float(self.x_2.get())
        for i in range(0, dl):
            a_linear = (y1_l-y2_l)/(x1_l-x2_l)
            b_linear = y1_l-(a_linear*x1_l)
            b = ((float(self.spectrum.energy[i])*a_linear)+b_linear)
            background.append(b)
            c_b = spectrum_counts[i] - b
            counts_b.append(c_b)
        self.spectrum.background = background
        self.spectrum.final_yield = counts_b
        label_rb = rb
    plot_background(self, label_rb) 

def background_list(self): #, sg, rb, wavelet_family_names, wava_names, level_number, maxiter, oM, a0, a1, a2, a3, a4, x_1, x_2, y_1, y_2):
    sg = int(self.SGf.get())
    lev_number = int(self.level_number.get())
    miter = int(self.maxiter.get())
    wava_name = self.wavelet_names.get()
    for sp in list_of_spectra_c:
        if sg == 1: #yes
            spectrum_counts = sp.filter
        if sg == 0: #no
            spectrum_counts = sp.counts
        dl=len(sp.energy[:])
        spectrum_counts_log = []
        for i in range(0, dl):
            if float(spectrum_counts[i]) > 0:
                spectrum_counts_log.append(log(float(spectrum_counts[i])))
            else:
                spectrum_counts_log.append(0)
        energy_arr = np.array(sp.energy)
        counts_arr = np.array(spectrum_counts_log)
        background = []
        counts_b = []
        rb = self.rb_entry.get()
        if rb == 'None':
            for i in range(0, dl):
                background.append(0)
                c_b = spectrum_counts[i]
                counts_b.append(c_b)
            sp.background = background
            sp.final_yield = counts_b
            label_rb = rb
        if rb == 'Wavelet':
            baseline = baseline_dwt(counts_arr, level = lev_number, max_iter = miter, wavelet = wava_name)
            for i in range(0, dl):
                b = exp(baseline[i])
                background.append(b)
                c_b = spectrum_counts[i] - b
                counts_b.append(c_b)
            sp.background = background
            sp.final_yield = counts_b
            label_rb = wava_name
        if rb == 'New':
            baseline = baseline_dwt(counts_arr, level = lev_number, max_iter = miter, wavelet = wava_name)
            for i in range(0, dl):
                b = exp(baseline[i])
                background.append(b)
                c_b = spectrum_counts[i] - b
                counts_b.append(c_b)
            sp.background = background
            sp.final_yield = counts_b
            label_rb = rb
        if rb == 'Polynomial':
            oM_in = int(self.oM.get())
            a0_in = float(self.a0.get())
            a1_in = float(self.a1.get())
            a2_in = float(self.a2.get())
            a3_in = float(self.a3.get())
            a4_in = float(self.a4.get())
            for i in range(0, dl):
                Rdata0 = np.array(spectrum_counts)
                maxima0 = argrelextrema(Rdata0, np.less, order=oM_in)
                lfa0 = list(maxima0[0])
                params, covariance = curve_fit(my_function, lfa0, Rdata0[lfa0], p0=[a0_in, a1_in, a2_in, a3_in, a4_in]) 
                a1_fit, a2_fit, a3_fit, a4_fit, a5_fit = params
                fb = my_function(i, a1_fit, a2_fit, a3_fit, a4_fit, a5_fit)
                if fb < 0:
                    fb = 0
                background.append(fb)
                c_b = spectrum_counts[i] - fb
                counts_b.append(c_b)
            sp.background = background
            sp.final_yield = counts_b
            label_rb = rb
        if rb == 'Linear':
            y1_l = float(self.y_1.get()) 
            y2_l = float(self.y_2.get())
            x1_l = float(self.x_1.get())
            x2_l = float(self.x_2.get())
            for i in range(0, dl):
                a_linear = (y1_l-y2_l)/(x1_l-x2_l)
                b_linear = y1_l-(a_linear*x1_l)
                b = ((float(sp.energy[i])*a_linear)+b_linear)
                background.append(b)
                c_b = spectrum_counts[i] - b
                counts_b.append(c_b)
            sp.background = background
            sp.final_yield = counts_b
            label_rb = rb
    show_plot_list_c_b(self)
    plot_background_list(self, self.root) 