from ROOT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from root_numpy import tree2array
def PlotFilter(X,dEE):
    Ylow  = np.power(X,2)*dEE[0][0] + X*dEE[0][1] + dEE[0][2]
    Yhigh = np.power(X,2)*dEE[0][3] + X*dEE[0][4] + dEE[0][5]
    return Ylow, Yhigh

def GetHistArray2D(hist):
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    n2    = (n[0])*(n[1])
    dX    = hist.GetXaxis().GetBinWidth(0)
    dY    = hist.GetYaxis().GetBinWidth(0)
    d     = hist.GetArray()

    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax-dX, n[1])
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1])
    Xaxis = (Xlow + Xhigh)/2

    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax-dY, n[0]-1)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-1)
    Yaxis = (Ylow + Yhigh)/2                     

    
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xaxis, Yaxis, Xmin, Xmax


f                 = TFile(sys.argv[1])
header = f.Get("header")
MiddleWedge = 0
for event in header: MiddleWedge = event.MiddleWedge
#MiddleWedge+=1    
hist1             = f.Get("REhist/ProfileE_Tot_f")
hist1,Xaxis,Yaxis,Xmin, Xmax = GetHistArray2D(hist1)

print len(Xaxis), Xaxis[0], Xaxis[-1]

for i in range(5):
    phantom   = f.Get("Phantom/nBricks "+str(i)+" calibration phantom thickness")
    data      = phantom.ProjectionX("","")
    data      = np.array(data)[1:-1]
    n         = phantom.GetNbinsX()
    p_dX      = phantom.GetXaxis().GetBinWidth(0)
    p_Xmin    = phantom.GetXaxis().GetXmin()
    p_Xmax    = phantom.GetXaxis().GetXmax()
    p_Xlow    = np.linspace(p_Xmin, p_Xmax-p_dX, n)
    p_Xhigh   = np.linspace(p_Xmin+p_dX, p_Xmax, n)
    p_Xaxis   = (p_Xlow + p_Xhigh)/2
    #p_Xaxis  += 1
    plt.plot(p_Xaxis,data)
print hist1.shape
print Xmin, Xmax
plt.imshow(hist1,origin='lower', extent= [Xmin, Xmax, Yaxis[0], Yaxis[-1]],aspect='auto')
plt.plot([MiddleWedge,MiddleWedge],[0,1000],'r-')
plt.show()
