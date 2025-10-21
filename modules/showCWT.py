

from tkinter import ttk, Tk, Frame,Menu
import tkinter as tk
import tkinter.font as tkFont
from tkinter import messagebox
from tkinter import filedialog as fd
#from ttkthemes import themed_tk as tk

from modules.labeledInput import LabeledInput
import numpy as np
import pywt as pwt
from scipy import signal

from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
from matplotlib.backend_tools import ToolBase
from matplotlib import gridspec
import numpy as np
import pyperclip as pc
import pandas as pd

from mpl_toolkits.axes_grid1 import make_axes_locatable


cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma', 'cividis']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper'])]




class showCWT(tk.Toplevel):
    
    X,Y=None,None

    __grid__=''
    
    __menuCheckButton__=''
    dataInfo={}
   
    
    
    def __init__(self,parent,*args,**kwargs):
        super().__init__(parent,*args,**kwargs)
                                        
        
        self.__defFont__=tkFont.Font(family="FreeMono",size=12)
        self.__titFont__=tkFont.Font(family="Helvetica",size=14)
        bgcolor='bisque3'
        #----------------------------------------------------------------------
        #------- FRAME TOP ----------------
        self.pltFrame=tk.LabelFrame(self,text="",height=500)
        self.pltFrame.pack(fill='both',expand=True)
        #self.pltFrame.
        

        self.__fig__ = Figure(figsize = (4, 3), dpi = 150)         
                        
        canvas = FigureCanvasTkAgg(self.__fig__, master = self.pltFrame)         
        toolbar=NavigationToolbar2Tk(canvas)
        toolbar.update()
        canvas._tkcanvas.pack(fill=tk.BOTH,expand=1)
        
        
        #----------------------------------------------------------------------
        #------- FRAME BOTTOM ----------------
        stFrame=tk.LabelFrame(self,text="CWT settings",bg=bgcolor,height=150)
        stFrame.pack(fill='both',side=tk.BOTTOM,expand=True)
        stFrame.pack_propagate(False)
        stFrame.grid_propagate(False)
        
        #.......... families ..................
        self.WletFamily=tk.StringVar()
        self.WletFamily.set("mexh");
         
        cwavelist=pwt.wavelist(kind='continuous')
     
        inarg={'font': self.__defFont__,
                'width': 10,
                'justify': 'center',
                'values': cwavelist
                }
        lbarg={'font': self.__defFont__ , 
                'background' : bgcolor,
                'justify': 'left'
                }
        lecmFm=LabeledInput(stFrame,label='family ',
                           input_class=ttk.Combobox, input_args=inarg,
                           input_var=self.WletFamily,
                           label_args=lbarg)
        lecmFm.grid(row=0,column=0,columnspan=2,padx=10,pady=10)
        
        
        #.......... widths ..................
        self.WletWidths=tk.StringVar()
        self.WletWidths.set("1:0.25:8");
        
        
        inarg={'font': self.__defFont__,
                'width': 10,
                'justify': 'center',                
                }
        lbarg={'font': self.__defFont__ , 
                'background' : bgcolor,
                'justify': 'left'
                }
        
        leEdit=LabeledInput(stFrame,label='widths',
                            input_class=ttk.Entry,
                            input_args=inarg, input_var=self.WletWidths,
                            label_args=lbarg
                            )
        leEdit.grid(row=0,column=2,columnspan=2,padx=10,pady=10)
        
        
        #.......... thresh ..................
        self.WletThresh=tk.StringVar()
        self.WletThresh.set("0.003125");
        
        
        inarg={'font': self.__defFont__,
                'width': 10,
                'justify': 'center',                
                }
        lbarg={'font': self.__defFont__ , 
                'background' : bgcolor,
                'justify': 'left'
                }
        
        leEdit=LabeledInput(stFrame,label='threshold',
                            input_class=ttk.Entry,
                            input_args=inarg, input_var=self.WletThresh,
                            label_args=lbarg
                            )
        leEdit.grid(row=0,column=4,columnspan=2,padx=10,pady=10)
        
        
        brun=tk.Button(stFrame,text='replot',bg=bgcolor,command=self.plot)
        brun.grid(row=0,column=6,columnspan=1,padx=10,pady=10)
        
        #------------ colormaps
        #.......... families ..................
        self.colormaps=tk.StringVar()
        self.colormaps.set("Oranges");
         
        cmaps=['viridis', 'plasma', 'inferno', 'magma', 'cividis',
                'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                    'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                    'hot', 'afmhot', 'gist_heat', 'copper']
     
        inarg={'font': self.__defFont__,
                'width': 10,
                'justify': 'center',
                'values': cmaps
                }
        lbarg={'font': self.__defFont__ , 
                'background' : bgcolor,
                'justify': 'left'
                }
        leCmColors=LabeledInput(stFrame,label='colormaps ',
                           input_class=ttk.Combobox, input_args=inarg,
                           input_var=self.colormaps,
                           label_args=lbarg)
        leCmColors.grid(row=1,column=0,columnspan=2,padx=10,pady=10)
        #----------------------------------------------------------------------
        
        self.initMenu()        
        self.geometry("800x600")
        self.title('CWT')
        #self.iniDir='/home/fizyk/python/tloremover-git/inputData/'
        self.parent=parent
        

        
    
        


        
    def initMenu(self):               
        menubar = Menu(self)
        

        fileMenu = Menu(menubar)#font=("",14)
        #fileMenu.add_command(label="Open")
        fileMenu.add_command(label="Save",command=self.dataSave)
        fileMenu.add_command(label="Exit")
        
        optMenu =Menu(menubar)
        optMenu.add_command(label="data  to  clipboard",command=self.dataToClipboard)
        optMenu.add_command(label="data from clipboard",command=self.dataFromClipboard)
                        
        self.menuXAxisLog=tk.BooleanVar()
        self.menuXAxisLog.set(False)
        
        optMenu.add_checkbutton(label="X-axis log",onvalue=1,offvalue=0,
                                command=self.eventXAxisLog,variable=self.menuXAxisLog)
        
        self.menuYAxisLog=tk.BooleanVar()
        self.menuYAxisLog.set(False)
        
        optMenu.add_checkbutton(label="Y-axis log",onvalue=1,offvalue=0,
                                command=self.eventYAxisLog,variable=self.menuYAxisLog)  
        
        
        
        
        helpMenu=Menu(menubar)
        helpMenu.add_command(label="Authors/Credits")
        
        plotMenu=Menu(menubar)
        
        self.menuLegend=tk.BooleanVar()
        self.menuLegend.set(False)
        
        plotMenu.add_checkbutton(label="Legend",onvalue=0,offvalue=0,command=self.pltLegend,variable=self.menuLegend)
        plotMenu.add_checkbutton(label="Grid",onvalue=1,offvalue=0,command=self.pltGrid)
        
        
        
        ##self.menuCheckButton['legend']=mlegend
                
        menubar.add_cascade(label="File", menu=fileMenu)      
        menubar.add_cascade(label="Options",menu=optMenu)
        menubar.add_cascade(label="Plot",menu=plotMenu)
        menubar.add_cascade(label="Help",menu=helpMenu)
        
        self.config(menu=menubar)
        
        
        
        
    def dataSave(self):
        filetypes = ( ('results','*.cwt'),('All files', '*.*') )
        outputFileName=fd.asksaveasfilename(title='Save a file',
                              initialdir=self.iniDir,
                              filetypes=filetypes)
                    
        #self.data2D.savetxt(outputFileName)
        
            
            
       #-------------------------------------------------- 
                
# sudo apt install xclip        
# pip install paperclip

    def dataToClipboard(self):
        #df = pd.DataFrame(self.data2D)
        df.to_clipboard(sep='\t') 
        
    def dataFromClipboard(self):
        df=pd.read_clipboard().to_numpy()
        
        plt=self.__plt__
        plt.plot(df[:,0],df[:,1],'.k',markersize=0.5,label='data from clipboard')
        
        self.__fig__.canvas.draw()
    #--------------------------------------------------   
    
    def eventXAxisLog(self):
        if self.menuXAxisLog.get():
            self.__plt__.set_xscale('log',base=10)
        else:
            self.__plt__.set_xscale('linear')
            
        self.__fig__.canvas.draw()
            
            
            
    def eventYAxisLog(self):
        if self.menuYAxisLog.get():
            self.__plt__.set_yscale('log',base=10)
        else:
            self.__plt__.set_yscale('linear')        
            
        self.__fig__.canvas.draw()            
       


    def onExit(self):
        self.quit()         
    #--------------------------------------------------   


        
    def plot(self,plt_args=None,plt_kwargs=None): 
        
        if self.X is  None:
            return
        
        fig=self.__fig__          
        
        if fig is not None:
            fig.clear()        
        '''
        spec = gridspec.GridSpec(ncols=1, nrows=2, 
                         hspace=0.1, height_ratios=[3,1])
        '''
        
        ax = fig.add_gridspec(3,8)
        a0 = fig.add_subplot(ax[0:2, 0:7])                    
        a1 = fig.add_subplot(ax[2, 0:7])                    
        a2 = fig.add_subplot(ax[0:,7])                    
        #fig.subplots_adjust(right=0.95,sharex=True)

        plt_args=plt_args or {}
        plt_kwargs=plt_kwargs or {}
        
        wFamily=self.WletFamily.get()
        wWidths=self.WletWidths.get().split(':')
        wThresh=self.WletThresh.get()

        
        minWl,stWl,maxWl=float(wWidths[0]),float(wWidths[1]),float(wWidths[2])
        widths=2.0**np.arange(minWl,maxWl,stWl)
        X,Y=self.X, self.Y
        
        dX=(X.max()-X.min())/X.shape[0]
        #print('dX ',dX,'  ',X[1]-X[0])
        
        cwtmatr,freqs=pwt.cwt(Y,widths,wFamily)
        freqs/=dX
        #print('fr=',freqs)
        fmin,fmax=freqs[-1],freqs[0]

        thr=float(wThresh)        
        if thr<0:
            thr=0
            
        izero=np.where(cwtmatr<thr)
        cwtmatr[izero]=0

        
        #bnd=5
        
        pcm = a0.pcolormesh(X, freqs, (cwtmatr)**0.5,cmap=self.colormaps.get())
        #cbx=fig.add_subplot(a0)        
        fig.colorbar(pcm,ax=a2,shrink=0.6)        
        #divider = make_axes_locatable(a0)
        #cax = divider.append_axes('right', size="7%", pad=0.2,)
        #fig.colorbar(pcm,cax=cax)
        
        a0.set_yscale("log",base=2)
        #a0.set_xlabel('X')
        a0.set_ylabel('frequency [a.u.]')
        #a0.set_xlim([0,3500])
        a0.get_xaxis().set_ticks([])
        
        a1.plot(X,fmax*0.5*Y/Y.max()+fmin,'-k',linewidth=0.75)      
        a1.spines['top'].set_visible(False)
        a1.spines['right'].set_visible(False)
        a1.set_xlim([X[0],X[-1]])
          
        
        a2.axis('off')
                            
        fig.canvas.draw()
        
    #--------------------------------------------------
        
        
                
    def pltGrid(self):        
        self.__plt__.grid()
        self.__fig__.canvas.draw()
        
        
    def pltLegend(self):          
        self.__plt__.legend().set_visible(self.menuLegend.get())            
        self.__fig__.canvas.draw()
        

    def pltTitle(self,ptitle='title'):
        self.__plt__.set_title(ptitle)
        
    
    
    def graphOnClick(self,event):
        return
    '''
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
         ('double' if event.dblclick else 'single', event.button,
          event.x, event.y, event.xdata, event.ydata))
        '''
            

        
        
