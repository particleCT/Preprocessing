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

    d     = hist.GetArray()
    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax-dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    Xaxis = (Xlow + Xhigh)/2                 
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xaxis


f = TFile(sys.argv[1])

for i in range(1,5):
    hist1  = f.Get("dEE/dE-EStage_"+str(i))
    header = f.Get("header")
    dEE = tree2array(header,branches=["dEElow_Stage"+str(i)+"_0", "dEElow_Stage"+str(i)+"_1", "dEElow_Stage"+str(i)+"_2",
                                      "dEEhigh_Stage"+str(i)+"_0", "dEEhigh_Stage"+str(i)+"_1", "dEEhigh_Stage"+str(i)+"_2"])
    hist1,Xaxis = GetHistArray2D(hist1)
    #plt.imshow(hist1)
    X = np.linspace(0,300)

    Ylow, Yhigh = PlotFilter(X,dEE)
    plt.imshow(hist1, origin='lower', norm=LogNorm(),extent=[Xaxis[0], Xaxis[-1], Xaxis[0], Xaxis[-1]])
    plt.title("Stage "+str(i))
    plt.plot(Ylow, X ,'k--',lw=2)
    plt.plot(Yhigh, X ,'r-',lw=2)
    plt.show()


plt.show()
