from root_numpy import tree2array
from ROOT import *
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from matplotlib.colors import LogNorm
def GetHistArray2D(hist):
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    n2    = (n[0])*(n[1])
    d     = hist.GetArray()
    dX    = hist.GetXaxis().GetBinWidth(0)
    dY    = hist.GetYaxis().GetBinWidth(0)

    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax+dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    #Xaxis = (Xlow + Xhigh)/2

    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax+dY, n[0]-2)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-2)
    #Yaxis[1:-2] = (Ylow + Yhigh)/2
    #Yaxis[0]  = Ymin
    #Yaxis[-1] = Ymax

    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xlow,Ylow

def GetProjectionArray2D(hist):
    hist  = hist.ProjectionXY("","")    
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    n2    = (n[0])*(n[1])
    d     = hist.GetArray()
    dX    = hist.GetXaxis().GetBinWidth(0)
    dY    = hist.GetYaxis().GetBinWidth(0)

    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax+dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    #Xaxis = (Xlow + Xhigh)/2

    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax+dY, n[0]-2)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-2)
    #Yaxis[1:-2] = (Ylow + Yhigh)/2
    #Yaxis[0]  = Ymin
    #Yaxis[-1] = Ymax

    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xlow,Ylow


def GetCalibrationCurve(hist,Yaxis):
    NY, NX = hist.shape
    Max = []
    Std = []
    X   = []
    for n in range(NX):
        try:
            y= hist[:,n]
            maxval = np.max(y)
            maxid  = np.argmax(y)
            X.append(Yaxis[n])
            popt,pcov = curve_fit(gaus,Yaxis,y,p0=[maxval,Yaxis[maxid],10], maxfev = 10000)
            Max.append(popt[1])
            Std.append(popt[2])
        except:
            Max.append(0)
            Std.append(0)

    return np.array(Max),np.array(Std)
    
f          = TFile(sys.argv[1],"update")
calibGraph = f.Get("calWEPL/NewRange_XY")

##Load the phasespace data
tree   = f.Get("phase")
array  = tree2array(tree,branches=['y0','y1','z0','z1','wepl','E','dEEFilter','ThresholdFilter','MaxEnergyTransFilter','PriorFilter',])
#array  = array[np.where( (array['dEEFilter']==1))]# & (array['MaxEnergyTransFilter']==1) & (array['ThresholdFilter']==1))]
array  = array[np.where( (array['PriorFilter']==1))]
WET    = array['wepl']
y1     = array['y1']
z1     = array['z1']
y0     = array['y0']
z0     = array['z0']

Y      = (y1 + y0)/2
Z      = (z1 + z0)/2

bins,   y_edges,z_edges = np.histogram2d(Z,Y, weights = WET,bins=(300,300))
entries,y_edges,z_edges = np.histogram2d(Z,Y,bins=(300,300))

z_center = (z_edges[:-1] + z_edges[1:])/2
Final = bins/entries
Final[np.isnan(Final)]=0
plt.imshow(Final)
plt.show()
Profile = np.mean(Final[25:200],axis=0)
plt.plot(z_center,Profile)

## Load the perfect Bricks
calib = TFile(sys.argv[2])
profile =  calib.Get("Phantom/nBricks "+sys.argv[3]+" calibration phantom thickness")
Y = []
X = []
for i in range(profile.GetNbinsX()):
        X.append(profile.GetXaxis().GetBinCenter(i))
        Y.append(profile.GetBinContent(i))

plt.plot(X,Y)
plt.show()

f = interpolate.interp1d(z_center,Profile,fill_value="extrapolate")
plt.plot(X, Y-f(X))
#plt.plot(Y)
#plt.plot(f(z_center))
plt.ylim([-2,2])
plt.show()
