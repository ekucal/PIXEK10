import numpy as np
import re

class SpectrumData:
    fileTypes = (
                 ('pixie files','*.spe'),
                 ('binary npy','*.npy'),
                 ('data files', '*.dat'),
                 ('diffdata files', '*.diff'),
                 ('All files', '*.*'))

    reducTypes=['none','sqrt','log10']
    plotScheme='.k'

    #input file name
    ifileName=''
    #output file name
    ofileName=''
    #spectrum data
    Xi,Yi = None,None
    Xc = None
    Xr, Yr = None, None
    Xf, Yf = None, None
    Xb, Yb = None, None
    Info = None
    # Spectrum Info:
    Material, CrystalOrientation, AxialChannel, TiltAngle, Plane = None, None, None, None, None
    Filter, Background = None, None
    Dose, Dose_e = None, None

    #select range:  from, to
    _from,_to=0.,-1.
    #rtypes
    rtype=reducTypes[0]

    ylabel='Yield'

    def clear(self):
        self.ifileName=''
        self.ofileName=''
        self.Xi,self.Yi, self.Info=None, None, None
        self.Xc = None
        self.Xr, self.Yr = None, None
        self.Xf, self.Yf = None, None
        self.Xb, self.Yb = None, None
        self.Y = None
        self.Material, self.CrystalOrientation, self.AxialChannel, self.TiltAngle, Plane = None, None, None, None, None
        self.Filter, self.Background = None, None
        self.Dose, self.Dose_e = None, None

    def loadSpectrum(self,fileName): # tu dodano przekazywanie Info do zapisu
        if fileName=='':
            print('no file name given')
            return False
        self.clear()
        data2D=None
        if re.search('\.dat$',fileName):
            data2D=np.loadtxt(fileName)
        else:
            if re.search('\.npy',fileName):
                data2D=np.load(fileName).T
            else:
                if re.search('\.spe',fileName):
                    fin=open(fileName,'r')
                    for _ in range(0,7): fin.readline()

                    f,t=fin.readline().split()
                    data2D=np.ndarray((int(t)-int(f),2),float)
                    fpos=0

                    for line in fin:
                        splitLine=line.split()
                        for v in splitLine:
                            data2D[fpos,0]=fpos
                            data2D[fpos,1]=float(v)
                            fpos+=1

                    fin.close()
                else:
                    if re.search('\.diff',fileName):
                        data2D=np.loadtxt(fileName)[:,1:5:2]
                    else:
                        print('ERROR: unknown extension')
        if data2D is None:
            print('ERROR: no data')
            return False
        if self.Dose or self.Dose_e is None:
            doseScale = 1
        else:
            doseScale = float(self.Dose)/float(self.Dose_e)
        self.Xi=data2D[:,0]
        self.Yi=(data2D[:,1])*doseScale
        self.ifileName=fileName
        return True

    def dataSize(self):
        Xiarray=np.array(self.Xi)
        return Xiarray.shape[0]

    def empty(self):
        return True if self.Xi is None else False

    











