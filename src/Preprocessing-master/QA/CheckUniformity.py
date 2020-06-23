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
    Xlow  = np.linspace(Xmin, Xmax-dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    Xaxis = (Xlow + Xhigh)/2

    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax-dY, n[0]-2)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-2)
    Yaxis = (Ylow + Yhigh)/2                 
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xaxis,Yaxis

def GetProfileArray2D(hist):
    hist  = hist.ProjectionXY("")
    n     = [hist.GetNbinsY()+2,hist.GetNbinsX()+2]
    n2    = (n[0])*(n[1])
    dX    = hist.GetXaxis().GetBinWidth(0)
    dY    = hist.GetYaxis().GetBinWidth(0)
    
    d     = hist.GetArray()
    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax-dX, n[1]-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n[1]-2)
    Xaxis = (Xlow + Xhigh)/2                 
    Ymin  = hist.GetYaxis().GetXmin()
    Ymax  = hist.GetYaxis().GetXmax()
    Ylow  = np.linspace(Ymin, Ymax-dY, n[0]-2)
    Yhigh = np.linspace(Ymin+dY, Ymax, n[0]-2)
    Yaxis = (Ylow + Yhigh)/2                 
    d.SetSize(n2)
    histA = np.array(d)
    histA = np.reshape(histA,n)[1:-1,1:-1]
    return histA,Xaxis,Yaxis

f = TFile(sys.argv[1])

for stage in range(0,5):
    hist  = f.Get("E_map_front_TS_"+str(stage))
    hist,Xaxis,Yaxis = GetProfileArray2D(hist)
    print Xaxis
    idsX = np.where((Xaxis>-100) & (Xaxis<100))
    idsY = np.where((Yaxis>-30) & (Yaxis<30))
    
    ## X Direction
    MeanX = np.mean(hist[idsY[0]][:,idsX[0]],axis=0)
    plt.plot(Xaxis[idsX],np.mean(hist[idsY[0]][:,idsX[0]],axis=0),'bo',label="Stage "+str(stage))
    plt.plot([min(Xaxis[idsX]),max(Xaxis[idsX])] ,[np.mean(MeanX),np.mean(MeanX)],'k-')
    plt.plot([min(Xaxis[idsX]),max(Xaxis[idsX])] ,[np.mean(MeanX)+2,np.mean(MeanX)+2],'r--')
    plt.plot([min(Xaxis[idsX]),max(Xaxis[idsX])] ,[np.mean(MeanX)-2,np.mean(MeanX)-2],'r--')
    plt.title("T direction")
    plt.legend()
    plt.show()
    
    ## Y Direction
    plt.title("V direction")
    MeanY = np.mean(hist[idsY[0]][:,idsX[0]],axis=1)
    plt.plot(Yaxis[idsY],np.mean(hist[idsY[0]][:,idsX[0]],axis=1),'bo',label="Stage "+str(stage))
    plt.plot([min(Yaxis[idsX]),max(Yaxis[idsX])] ,[np.mean(MeanY),np.mean(MeanY)],'k-')
    plt.plot([min(Yaxis[idsX]),max(Yaxis[idsX])] ,[np.mean(MeanY)+2,np.mean(MeanY)+2],'r--')
    plt.plot([min(Yaxis[idsX]),max(Yaxis[idsX])] ,[np.mean(MeanY)-2,np.mean(MeanY)-2],'r--')
    plt.legend()
    plt.show()

