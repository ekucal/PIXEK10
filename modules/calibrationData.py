
import numpy as np


class XaxisMapping:
    methods=['Ignore','E=ax+b','Peaks']
    #method=methods[0]
    def __call__(self,x,a,b,method):
        if method==self.methods[0]:
            return x
        if method==self.methods[1]:
            return np.polyval(np.array([a,b]),x)
        if method==self.methods[2]:
            Xc = []
            for xi in x:
                Xc.append(a*xi+b)
            return Xc
    def isIgnore(self,method):
        return True if method==self.methods[0] else False





