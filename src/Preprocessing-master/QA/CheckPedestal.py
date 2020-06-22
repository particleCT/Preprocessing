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
print directory.ls()
for subdirkey in directory.GetListOfKeys():
    #if(subdirkey.GetName().count(sys.argv[2])):
    subdir = f.Get("Pedestals/"+subdirkey.GetName())
    for files in subdir.GetListOfKeys():
        if(files.GetName().count("Pedestal") and files.GetName().count(sys.argv[3])):
            if(not files.GetName().count("Out")):
                if(not files.GetName().count("In")):
                    hist,Xaxis = GetHistArray(f.Get("Pedestals/"+subdirkey.GetName()+"/"+files.GetName()))            
                    plt.plot(Xaxis,hist/np.max(hist),label=files.GetName()+subdirkey.GetName())
                    #plt.plot(Xaxis,hist,label=files.GetName())
                   

plt.title(sys.argv[2])
plt.legend()
plt.show()

