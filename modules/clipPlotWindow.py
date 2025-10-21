

from tkinter import ttk, Tk, Frame,Menu
import tkinter as tk

from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
import matplotlib.pyplot as plt
from matplotlib.backend_tools import ToolBase
import matplotlib.image as mpimg
from PIL import Image


class ExtNavigationToolbar(NavigationToolbar2Tk):
    
    def __init__(self, canvas, parent=None,plt=None):
        super(ExtNavigationToolbar, self).__init__(canvas, parent)
        
        self.plt=plt        
        self.parent=parent
        

        


class clipPlotWindow(tk.Toplevel):

    
    __plt__=''
    __fig__=''
    __grid__=''
    img=None
    
    
    
    def __init__(self,parent,*args,**kwargs):
        super().__init__(parent,*args,**kwargs)
        
        wposition='400x300+'+str(parent.winfo_x()+600)+"+"+str(parent.winfo_y())
                
        self.geometry(wposition)
                        
        self.__fig__ = Figure(figsize = (4,3), dpi = 150) 
        fig=self.__fig__
                        
        self.__plt__=fig.add_subplot(111)        
        plt=self.__plt__
        
        fig.subplots_adjust(top=1,bottom=0,left=0,right=1,hspace=0,wspace=0)

        canvas = FigureCanvasTkAgg(fig, master = self)        
        #canvas.get_default_filename=lambda: self.imageDir.get()
        toolbar=ExtNavigationToolbar(canvas,self,self.__plt__)
        toolbar.update()
        canvas._tkcanvas.pack(fill=tk.BOTH,expand=1)
        
            
    def drawImage(self,fp):        
        with fp:
            img=mpimg.imread(fp,format='jpeg')
        
        self.__plt__.imshow(img)
        self.__plt__.axis('off')               
        self.__fig__.canvas.draw()
        
        
        

    