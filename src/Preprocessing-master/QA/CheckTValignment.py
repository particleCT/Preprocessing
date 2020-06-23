from ROOT import *
from root_numpy import tree2array
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import sys
import numpy as np
def GetHistArray(hist):
    n     = hist.GetNbinsX()+2
    dX    = hist.GetXaxis().GetBinWidth(0)
    d     = hist.GetArray()
    
    Xmin  = hist.GetXaxis().GetXmin()
    Xmax  = hist.GetXaxis().GetXmax()
    Xlow  = np.linspace(Xmin, Xmax-dX, n-2)
    Xhigh = np.linspace(Xmin+dX, Xmax, n-2)
    Xaxis = (Xlow + Xhigh)/2
    d.SetSize(n)
    histA = np.array(d)[1:-1]
    return histA,Xaxis


### T Position
f=TFile(sys.argv[1])
fig = plt.figure()
gs  = gridspec.GridSpec(ncols=1, nrows=4)
plt.title("T Position")
plt.subplots_adjust(hspace=0.001)
for histid in range(0,4):
    hist = f.Get("Alignement/TalignementPosT"+str(histid))
    hist,Xaxis = GetHistArray(hist)
    ax = fig.add_subplot(gs[histid,0])
    mean = np.average(Xaxis,weights=hist)
    meanid = np.argmin(np.abs(Xaxis-mean))
    label= "T%d Mean:%.3f"%tuple([histid,mean])
    ax.annotate(label, (0.01, 0.9), xycoords='axes fraction', va='center')
    ax.plot(Xaxis,hist)
    ax.plot([0,0],[np.min(hist),np.max(hist)],'r-', label='Center')
    ax.plot([Xaxis[meanid],Xaxis[meanid]],[np.min(hist),np.max(hist)],'k-',label='Mean')
plt.legend(frameon=False)
plt.show()

### T Direction
fig = plt.figure()
gs  = gridspec.GridSpec(ncols=1, nrows=2)
plt.title("T Direction")
plt.subplots_adjust(hspace=0.001)
for histid in range(0,2):
    if(histid==0): T0 = "T0";T1="T1"
    else : T0 = "T2"; T1="T3"
    hist = f.Get("Alignement/TalignementDir"+T0+T1)
    hist,Xaxis = GetHistArray(hist)
    ax = fig.add_subplot(gs[histid,0])
    mean = np.average(Xaxis,weights=hist)
    meanid = np.argmin(np.abs(Xaxis-mean))
    label= T1+"-"+T0+" Mean:%.3f"%tuple([mean])
    ax.annotate(label, (0.01, 0.9), xycoords='axes fraction', va='center')
    ax.plot(Xaxis,hist)
    ax.plot([0,0],[np.min(hist),np.max(hist)],'r-',label='Center')
    ax.plot([Xaxis[meanid],Xaxis[meanid]],[np.min(hist),np.max(hist)],'k-',label='Mean')
plt.legend(frameon=False)
plt.show()


### V Position
fig = plt.figure()
gs  = gridspec.GridSpec(ncols=1, nrows=4)
plt.title("V direction")
plt.subplots_adjust(hspace=0.001)
for histid in range(0,4):
    hist = f.Get("Alignement/ValignementPosV"+str(histid))
    hist,Xaxis = GetHistArray(hist)
    ax = fig.add_subplot(gs[histid,0])
    mean = np.average(Xaxis,weights=hist)
    meanid = np.argmin(np.abs(Xaxis-mean))
    label= "V%d Mean:%.3f"%tuple([histid,mean])
    ax.annotate(label, (0.1, 0.5), xycoords='axes fraction', va='center')
    ax.plot(Xaxis,hist)
    ax.plot([0,0],[np.min(hist),np.max(hist)],'r-',label='Center')
    ax.plot([Xaxis[meanid],Xaxis[meanid]],[np.min(hist),np.max(hist)],'k-',label='Mean')
plt.legend(frameon=False)
plt.show()

### T Direction
fig = plt.figure()
gs  = gridspec.GridSpec(ncols=1, nrows=2)
plt.title("V Direction")
plt.subplots_adjust(hspace=0.001)
for histid in range(0,2):
    if(histid==0): V0 = "V0";V1="V1"
    else : V0 = "V2"; V1="V3"
    hist = f.Get("Alignement/ValignementDir"+V0+V1)
    hist,Xaxis = GetHistArray(hist)
    ax = fig.add_subplot(gs[histid,0])
    mean = np.average(Xaxis,weights=hist)
    meanid = np.argmin(np.abs(Xaxis-mean))
    label= V1+"-"+V0+" Mean:%.3f"%tuple([mean])
    ax.annotate(label, (0.01, 0.9), xycoords='axes fraction', va='center')
    ax.plot(Xaxis,hist)
    ax.plot([0,0],[np.min(hist),np.max(hist)],'r-',label='Center')
    ax.plot([Xaxis[meanid],Xaxis[meanid]],[np.min(hist),np.max(hist)],'k-',label='Mean')
plt.legend(frameon=False)
plt.show()

