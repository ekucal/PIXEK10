
import numpy as np
import re
from scipy.signal import savgol_filter


class Filter:


    methods=['Ignore','moving average','Savitzky-Golay']
    method=methods[0]
    plotScheme='-g'


    size=10
    polyOrder=3
    deriv=0

    def __call__(self, Y, kernelSize, method, Order, Deriv):
        self.method=method
        if method==self.methods[0]:
            return Y

        self.size=kernelSize

        if method==self.methods[1]:
            kernel=np.ones((kernelSize,),float)/kernelSize
            return np.convolve(kernel,Y,'same')

        if method==self.methods[2]:
            self.polyOrder=Order
            self.deriv=Deriv
            return savgol_filter(Y, kernelSize, Order,deriv=Deriv)

        print('ERROR: unknown filtering method')
        return Y

    def isIgnore(self,method):
        return True if method==self.methods[0] else False




