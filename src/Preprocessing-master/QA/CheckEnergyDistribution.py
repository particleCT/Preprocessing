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

directory = f.Get("Pedestals")
for stage in range(0,5):
    for subdirkey in directory.GetListOfKeys():
        subdir = f.Get("Pedestals/"+subdirkey.GetName())
        for files in subdir.GetListOfKeys():        
            if(files.GetName() == "EnergyDistribution_"+str(stage)):
                hist,Xaxis = GetHistArray(f.Get("Pedestals/"+subdirkey.GetName()+"/"+files.GetName()))
                plt.plot(Xaxis,hist/np.max(hist),label=subdirkey.GetName())
        

    plt.title(stage)
    plt.legend()
    plt.show()

