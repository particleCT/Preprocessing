import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
import sys
def GetHistArray(hist):
    n     = hist.GetNbinsX()+2
    d     = hist.GetArray()
    d.SetSize(n)
    histA = np.array(d)[1:-1]
    Xmin = hist.GetXaxis().GetXmin()
    Xmax = hist.GetXaxis().GetXmax()
    Xaxis = np.linspace(Xmin,Xmax,n)[1:-1]
    
    return histA,Xaxis



f = TFile(sys.argv[1])
print f.ls()
header = f.Get("header")
peds = []
for event in header:
    peds.append(event.peds_0)
    peds.append(event.peds_1)
    peds.append(event.peds_2)
    peds.append(event.peds_3)
    peds.append(event.peds_4)

directory = f.Get("Pedestals")
print directory.ls()
for subdirkey in directory.GetListOfKeys():
    subdir = f.Get("Pedestals/"+subdirkey.GetName())
    for files in subdir.GetListOfKeys():
        if(files.GetName().count("Pedestal") and files.GetName().count(sys.argv[2])):
            if(not files.GetName().count("Out")):
                if(not files.GetName().count("In")):
                    hist,Xaxis = GetHistArray(f.Get("Pedestals/"+subdirkey.GetName()+"/"+files.GetName()))            
                    plt.plot(Xaxis,hist/np.max(hist),label=files.GetName()+subdirkey.GetName())

                    #plt.plot(Xaxis,hist,label=files.GetName())
                   
ids = np.argmin(np.abs(Xaxis-peds[int(sys.argv[2])]))
plt.plot(Xaxis[ids],hist[ids]/np.max(hist),'ro')
plt.title(sys.argv[2])
plt.legend()
plt.show()

